import hashlib
import json
import re

import pandas as pd
from biocypher import APIRequest, BioCypher
from .constants import REGISTRY_KEYS
from .mapping_utils import map_species_terms

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

def is_valid_peptide_sequence(seq: str) -> bool:    
    """Checks if a given sequence is a valid peptide sequence."""
    if isinstance(seq, str) and len(seq) > 2:
        return all([aa in AMINO_ACIDS for aa in seq])
    else:
        return False

def process_cdr3_sequence(seq: str, is_igh: bool = False) -> str | None:
    if seq is None:
        return None
    
    # Clean and normalize the sequence
    seq = str(seq).upper().strip().replace(" ", "").replace("\n", "")

    # Validate that the sequence contains only valid amino acids (optional: define valid AAs if needed)
    if not is_valid_peptide_sequence(seq):
        return None

    # Check if sequence has a valid CDR3 format
    starts_with_c = seq.startswith("C")
    ends_with_fw = seq.endswith("F") or (is_igh and seq.endswith("W"))
    
    if starts_with_c and ends_with_fw:
        return seq

    # Pad the sequence appropriately
    seq = seq.lstrip("C")  # remove leading C if already present
    if is_igh:
        seq = seq.rstrip("FW")  # remove existing F or W if present
        return f"C{seq}W"
    else:
        seq = seq.rstrip("F")
        return f"C{seq}F"

def process_epitope_sequence(seq: str | None) -> str | None:
    """Process a sequence string to remove non-peptidic characters and validate it."""
    if seq is None:
        return None
    seq = str(seq)
    result = seq.split("+")[0]  # split_epitope_sequence
    result = result.upper()
    result = "".join(result.split())

    return result

def normalise_gene_name(gene: str) -> str:
    if pd.isna(gene):
        return None
    gene = gene.strip()
    # Replace TCRA → TRA, TCRB → TRB, etc.
    gene = re.sub(r'^TCR([ABGD])', r'TR\1', gene)
    # Remove allele annotation like *01 or *01_F
    gene = re.sub(r'\*.*$', '', gene)
    return gene


def harmonize_sequences(table: pd.DataFrame) -> pd.DataFrame:
    """Preprocesses CDR3 sequences, epitope sequences, and gene names in a harmonized way."""
    # Clean CDR3 sequences (normalize junction_aas)
    for i in [1, 2]:
        cdr3_col = getattr(REGISTRY_KEYS, f"CHAIN_{i}_CDR3_KEY")
        type_col = getattr(REGISTRY_KEYS, f"CHAIN_{i}_TYPE_KEY")

        if cdr3_col in table.columns and type_col in table.columns:
            table[cdr3_col] = table.apply(
                lambda row: process_cdr3_sequence(
                    row[cdr3_col],
                    is_igh=(row[type_col] == "IGH")
                ),
                axis=1
            )

    # Clean epitope sequences
    if REGISTRY_KEYS.EPITOPE_KEY in table.columns:
        table[REGISTRY_KEYS.EPITOPE_KEY] = table[REGISTRY_KEYS.EPITOPE_KEY].apply(process_epitope_sequence)

    # Normalize V and J genes
    vj_genes_cols = [
        REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
        REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
        REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
        REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
    ]
    for col in vj_genes_cols:
        if col in table.columns:
            table[col] = table[col].apply(normalise_gene_name)


    # human_keys = [
    #     'HomoSapiens',
    #     'Homosapian',  # typo but likely refers to humans
    #     'Human',
    #     'http://purl.obolibrary.org/obo/NCBITaxon_9606',  # NCBI taxonomy ID for Homo sapiens
    #     'Homo sapiens (human)',
    # ]
    # # table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY] = table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY].replace(human_keys, 'HomoSapiens')
    # # table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY].replace(human_keys, 'HomoSapiens')
    # # table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY].replace(human_keys, 'HomoSapiens')

    # Apply mapping for Antigen species keys
    if REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY in table.columns:
        antigen_species_map = map_species_terms(table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY].dropna().unique().tolist())
        table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY] = table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY].map(antigen_species_map)

    # Clean antigen names
    # table[REGISTRY_KEYS.ANTIGEN_KEY] = table[REGISTRY_KEYS.ANTIGEN_KEY].apply(clean_antigen)


    return table

def get_iedb_ids_batch(bc: BioCypher, epitopes: list[str], chunk_size: int = 150) -> dict[str, int]:
    """Retrieve IEDB IDs for multiple epitopes using batched requests.

    First tries exact matches, then falls back to substring matches for unmatched epitopes.

    Args:
        bc: Biocypher instance for the donwnload
        epitopes: List of epitope sequences to query
        chunk_size: Size of chunks to break epitopes into (to avoid URL length limits)

    Returns:
        Dictionary mapping epitope sequences to their IEDB IDs (0 if not found)
    """
    base_url = "https://query-api.iedb.org/epitope_search"
    epitope_to_id = {}
    unmatched_epitopes = []
    epitopes = [e for e in epitopes if e is not None]

    # Step 1: Try exact matches first
    print("Mapping AA epitope sequences to IEDB IDs: exact matches...")

    for i in range(0, len(epitopes), chunk_size):
        chunk = epitopes[i : i + chunk_size]
        # print("Starting new batch...")

        epitope_matches = _get_epitope_data(bc, chunk, base_url, match_type="exact")

        # Map results to the dictionary
        for epitope in chunk:
            epitope_to_id[epitope] = f"no_iedb_iri:{epitope}"  # Default value if not found

        for match in epitope_matches:
            # Handle both possible API return formats for epitope sequences
            if match["structure_descriptions"]:
                epitope_seq = match["structure_descriptions"][0]
            elif "linear_sequence" in match:
                epitope_seq = match["linear_sequence"]
            else:
                continue

            if epitope_seq in epitope_to_id:  # Only update if it's one we requested
                epitope_to_id[epitope_seq] = f"iedb_iri:{match['structure_id']}"

    # Step 2: Collect epitopes without matches and try string matching
    unmatched_epitopes = [ep for ep, id_val in epitope_to_id.items() if id_val == f"no_iedb_iri:{ep}"]

    if unmatched_epitopes:
        print(
            f"Found {len(epitopes) - len(unmatched_epitopes)} exact IEDB ID matches.",
            f"Trying substring matches for {len(unmatched_epitopes)} remaining epitopes...",
        )
        chunk_size = chunk_size // 2

        for i in range(0, len(unmatched_epitopes), chunk_size):
            chunk = unmatched_epitopes[i : i + chunk_size]
            substring_matches = _get_epitope_data(bc, chunk, base_url, match_type="substring")

            for epitope in chunk:
                best_match = None
                min_length = float("inf")

                # Find matches where the epitope is a substring
                for match in substring_matches:
                    if epitope in match["linear_sequence"] and len(match["linear_sequence"]) < min_length:
                        best_match = match
                        min_length = len(match["linear_sequence"])

                # Update the dictionary if a match was found
                if best_match:
                    epitope_to_id[epitope] = f"iedb_iri:{best_match['structure_id']}"

    # Final statistics
    matched_count = sum(1 for id_val in epitope_to_id.values() if id_val.startswith("iedb_iri:"))
    print(
        f"Final results: {matched_count} of {len(epitopes)} epitopes matched to IEDB IDs ({matched_count / len(epitopes) * 100:.1f}%)"
    )
    return epitope_to_id


def _get_epitope_data(bc: BioCypher, epitopes: list[str], base_url: str, match_type: str = "exact") -> list[dict]:
    """Get epitope data.

    Args:
        bc: Biocypher instance for the donwnload
        epitopes: List of epitope sequences to query
        base_url: Base URL for the API endpoint
        match_type: Type of matching to perform ("exact" or "substring")

    Returns:
        List of epitope data dictionaries
    """
    if match_type == "exact":
        request_hash = hashlib.md5("_".join(sorted(epitopes)).encode()).hexdigest()
        request_name = f"iedb_exact_matches_{request_hash}"
        epitope_list = f"({','.join([f'{e}' for e in epitopes])})"
        check = f"linear_sequence=in.{epitope_list}"
        url = f"{base_url}?{check}&select=structure_id,structure_descriptions,linear_sequence&order=structure_id"
        print(f"Request URL: {url[:100]}..." if len(url) > 100 else f"Request URL: {url}")

    else:
        request_hash = hashlib.md5("_".join(sorted(epitopes)).encode()).hexdigest()
        request_name = f"iedb_substring_matches{request_hash}"
        conditions = [f"linear_sequence.ilike.*{e}*" for e in epitopes]
        check = f"or=({','.join(conditions)})"
        url = f"{base_url}?{check}&select=structure_id,structure_descriptions,linear_sequence&order=structure_id"

    try:
        iedb_request = APIRequest(
            name=request_name,
            url_s=[url],
            lifetime=30,
        )
        paths = bc.download(iedb_request)

        # Load the cached JSON file
        if paths and len(paths) > 0:
            with open(paths[0]) as f:
                return json.load(f)

    except Exception as e:
        print(f"API request failed: {e}")
        return []

# def clean_antigen(term):
#     return re.sub(r'\s*[\(\[].*?[\)\]]', '', term).strip()