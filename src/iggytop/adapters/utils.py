""" This module contains utility functions for harmonizing data for iggytop """

import gzip
import hashlib
import importlib.resources
import json
import os
import re
from datetime import datetime
from logging import getLogger
from typing import List

import pandas as pd
import requests
import scirpy as ir
import yaml
from biocypher import APIRequest, BioCypher
from scirpy.io._datastructures import AirrCell

from .constants import REGISTRY_KEYS
from .mapping_utils import map_antigen_names, map_species_terms

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def _set_up_config(output_format, cache_dir):
    # Load the base configuration
    with importlib.resources.open_text("iggytop.config", "biocypher_config.yaml") as file:
        config = yaml.safe_load(file)
    if output_format:
        # Check if the output format is allowed
        allowed_formats = ["airr", "networkx", "neo4j", "docker"]
        if output_format not in allowed_formats:
            raise ValueError(f"Invalid output_format: {output_format}. Allowed formats are: {allowed_formats}")
        # Modify the configuration
        if output_format == "neo4j" or output_format == "networkx":
            config["biocypher"]["offline"] = True

        config["biocypher"]["dbms"] = output_format
        if output_format == "docker":
            with importlib.resources.open_text("iggytop.config", "biocypher_docker_config.yaml") as d_file:
                config = yaml.safe_load(d_file)

    # Save the modified configuration to the cache directory
    modified_config_path = os.path.join(cache_dir, "biocypher_config.yaml")
    with open(modified_config_path, "w") as file:
        yaml.safe_dump(config, file)

    return modified_config_path


def _set_up_schema(cache_dir):
    # Load the schema configuration
    schema_config_path = os.path.join(cache_dir, "schema_config.yaml")
    with importlib.resources.open_text("iggytop.config", "schema_config.yaml") as schema_file:
        with open(schema_config_path, "w") as cache_schema_file:
            cache_schema_file.write(schema_file.read())

    return schema_config_path


def get_file_checksum(file_path: str) -> str | None:
    """Calculates the SHA256 checksum of a file."""
    if not os.path.exists(file_path):
        return None
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def _is_valid_peptide_sequence(seq: str) -> bool:
    """Checks if a given sequence is a valid peptide sequence."""
    if isinstance(seq, str) and len(seq) > 2:
        return all([aa in AMINO_ACIDS for aa in seq])
    else:
        return False


def _process_cdr3_sequence(seq: str, is_igh: bool = False) -> str | None:
    if seq is None:
        return None

    # Clean and normalize the sequence
    seq = str(seq).upper().strip().replace(" ", "").replace("\n", "")

    # Validate that the sequence contains only valid amino acids (optional: define valid AAs if needed)
    if not _is_valid_peptide_sequence(seq):
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


def _process_epitope_sequence(seq: str | None) -> str | None:
    """Remove flanking residues in epitope sequences."""
    if seq is None:
        return None
    seq = str(seq)
    result = seq.split("+")[0]  # split_epitope_sequence
    result = result.upper()
    result = "".join(result.split())

    return result


def _normalize_vdj_gene_name(gene: str) -> str:
    """Process VDJ-gene names to align with IMGT standards, skip alleles information"""
    if pd.isna(gene):
        return None
    gene = gene.strip()
    # Replace TCRA → TRA, TCRB → TRB, etc.
    gene = re.sub(r"^TCR([ABGD])", r"TR\1", gene)
    # Remove allele annotation like *01 or *01_F
    gene = re.sub(r"\*.*$", "", gene)

    return gene.strip()


def get_mhc_class(allele: str | None) -> str | None:
    """Find MHC class information from the MHC gene (allele) names.

    Args:
        allele: MHC allele name.

    Returns:
        MHC class (I or II) or None if not found.
    """
    if allele is None or not isinstance(allele, str):
        return None

    allele = allele.upper().strip()

    # If the allele is already exactly "I" or "II", return it
    if allele in ["I", "II", "1", "2"]:
        return "I" if allele in ["I", "1"] else "II"

    # Class II matches
    class_ii_indicators = [
        "HLA-DR",
        "HLA-DQ",
        "HLA-DP",
        "HLA-DM",
        "HLA-DO",
        "DRB",
        "DRA",
        "DQB",
        "DQA",
        "DPB",
        "DPA",
        "H-2A",
        "H-2E",
        "H2-A",
        "H2-E",
        "I-A",
        "I-E",
        "CLASS II",
        "MH2",
        "MHC2",
        "MHC-II",
        "CLASSII",
        "MHCII",
        "MHC II",
        "CD4",
        "DQ8",
        "DR3",
        "H-2G7",
    ]
    if any(x in allele for x in class_ii_indicators):
        return "II"

    # Class I matches
    class_i_indicators = [
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-E",
        "HLA-F",
        "HLA-G",
        "CD1",
        "MR1",
        "H-2K",
        "H-2D",
        "H-2L",
        "H-2B",
        "H-2Q",  # Added H-2B from list
        "H2-K",
        "H2-D",
        "H2-L",
        "H2-Q",
        "H2-B",
        "QA-",
        "QB-",
        "H2-M",
        "MAMU-A",
        "MAMU-B",
        "MAMU-I",
        "GAGA-BF",
        "CLASS I",
        "MH1",
        "MHC1",
        "MHC-I",
        "CLASSI",
        "MHCI",
        "MHC I",
        "CD8",
    ]
    if any(x in allele for x in class_i_indicators):
        return "I"

    # Generic hyphenated H-2 (Mouse) or HLA (Human) patterns if not already caught
    if any(x in allele for x in ["H-2", "H2", "HLA"]):
        return "I"

    return None


def get_tissue_source(tissue: str | None) -> str:
    """Standardize tissue source information while staying close to the original values.

    Args:
        tissue: Original tissue source from database.

    Returns:
        Standardized tissue name or original value.
    """
    if tissue is None or not isinstance(tissue, str):
        return "nan"
    else:
        return tissue.upper().strip()


def harmonize_sequences(bc, table: pd.DataFrame) -> pd.DataFrame:
    """
    Preprocesses CDR3 sequences, epitope sequences, and gene names in a harmonized way.
    The following steps are performed:
    1. Clean CDR3 sequences (normalizes junction_aas)
    2. Clean epitope sequences (remove flanking residues)
    3. Normalize VDJ-gene names to IMGT standards
    4. Add IEDB IRI and corresponding antigen information (species and antigen name) where missing
    5. Harmonize species terms for antigen species and receptor chain species

    """
    # Clean CDR3 sequences (normalize junction_aas)
    for i in [1, 2]:
        cdr3_col = getattr(REGISTRY_KEYS, f"CHAIN_{i}_CDR3_KEY")
        type_col = getattr(REGISTRY_KEYS, f"CHAIN_{i}_TYPE_KEY")

        if cdr3_col in table.columns and type_col in table.columns:
            table[cdr3_col] = table.apply(lambda row: _process_cdr3_sequence(row[cdr3_col], is_igh=(row[type_col] == "IGH")), axis=1)

    # Clean epitope sequences
    if REGISTRY_KEYS.EPITOPE_KEY in table.columns:
        table[REGISTRY_KEYS.EPITOPE_KEY] = table[REGISTRY_KEYS.EPITOPE_KEY].apply(_process_epitope_sequence)

    # Normalize V and J genes
    vj_genes_cols = [
        REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
        REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
        REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
        REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
    ]
    for col in vj_genes_cols:
        if col in table.columns:
            table[col] = table[col].apply(_normalize_vdj_gene_name)

    # Map epitope sequences to IEDB-IRI mapping + extract species names
    if REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY not in table.columns:
        valid_epitopes = table[REGISTRY_KEYS.EPITOPE_KEY].dropna().unique().tolist()
        if len(valid_epitopes) > 0:
            # Sent API request to get IEDB IRIss and antigen infirmation for epitopes
            epitope_map = get_iedb_ids_batch(bc, valid_epitopes)

        # Add column with IEDB IRIs corresponding to the epitope AA sequence
        iri_mapping = {epitope: data["iri"] for epitope, data in epitope_map.items()}
        table[REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY] = table[REGISTRY_KEYS.EPITOPE_KEY].replace(iri_mapping)

        # Fill missing antigen and antigen species pairs if at least one is missing using information from IEDB
        organism_mapping = {epitope: data["organism"] for epitope, data in epitope_map.items() if data["organism"] is not None}
        antigen_mapping = {epitope: data["antigen"] for epitope, data in epitope_map.items() if data["antigen"] is not None}

        # Check if columns exist, create them if they don't
        if REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY not in table.columns:
            table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY] = None
        if REGISTRY_KEYS.ANTIGEN_KEY not in table.columns:
            table[REGISTRY_KEYS.ANTIGEN_KEY] = None

        # Create boolean masks for missing values
        antigen_missing = table[REGISTRY_KEYS.ANTIGEN_KEY].isna()
        organism_missing = table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY].isna()
        at_least_one_missing = antigen_missing | organism_missing

        # Get epitopes where at least one value is missing
        epitopes_to_fill = table[REGISTRY_KEYS.EPITOPE_KEY][at_least_one_missing]

        # Fill both antigen and organism for rows where at least one is missing

        for idx in epitopes_to_fill.index:
            epitope = epitopes_to_fill.loc[idx]

            # If antigen or antigen species is missing, fill both antigen and antigen species with IEDB values
            if epitope in antigen_mapping:
                table.loc[idx, REGISTRY_KEYS.ANTIGEN_KEY] = antigen_mapping[epitope]
                table.loc[idx, REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY] = organism_mapping[epitope]

    # Harmonize/clean species terms for both, antigen species and receptor chain species, using rules defined in map_species_terms
    if REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY in table.columns:
        antigen_species_harmonized_map = map_species_terms(table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY].dropna().unique().tolist())
        table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY] = table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY].replace(antigen_species_harmonized_map)

    if REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY in table.columns:
        chain_1_species_map = map_species_terms(table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY].dropna().unique().tolist())
        table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY].replace(chain_1_species_map)

    if REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY in table.columns:
        chain_2_species_map = map_species_terms(table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY].dropna().unique().tolist())
        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY].replace(chain_2_species_map)

    # Clean/delete brackets from the antigen names
    antigen_names_clean = map_antigen_names(table[REGISTRY_KEYS.ANTIGEN_KEY].dropna().unique().tolist())
    table[REGISTRY_KEYS.ANTIGEN_KEY] = table[REGISTRY_KEYS.ANTIGEN_KEY].replace(antigen_names_clean)

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
    epitope_to_iedb = {}
    unmatched_epitopes = []
    epitopes = [e for e in epitopes if e is not None]

    # Step 1: Try exact matches first
    getLogger("biocypher").info("Mapping AA epitope sequences to IEDB IDs: exact matches...")

    for i in range(0, len(epitopes), chunk_size):
        chunk = epitopes[i : i + chunk_size]
        epitope_matches = _get_epitope_data(bc, chunk, base_url, match_type="exact")

        # Map results to the dictionary
        for epitope in chunk:
            epitope_to_iedb[epitope] = {
                "iri": f"seq:{epitope}",
                "antigen": None,
                "organism": None,
            }  # Default value if not found

        for match in epitope_matches:
            # Handle both possible API return formats for epitope sequences
            if match["structure_descriptions"]:
                epitope_seq = match["structure_descriptions"][0]
            elif "linear_sequence" in match:
                epitope_seq = match["linear_sequence"]
            else:
                continue

            if epitope_seq in epitope_to_iedb:  # Only update if it's one we requested
                antigens = match.get("curated_source_antigens")
                if antigens:
                    antigen = antigens[0].get("name")
                    organism = antigens[0].get("source_organism_name")
                else:
                    antigen = None
                    organism = None

                epitope_to_iedb[epitope_seq] = {
                    "iri": f"iedb:{match.get('structure_id')}",
                    "antigen": antigen,
                    "organism": organism,
                }

    # Step 2: Collect epitopes without matches and try string matching
    unmatched_epitopes = [ep for ep, info in epitope_to_iedb.items() if info["iri"] == f"seq:{ep}"]

    if unmatched_epitopes:
        getLogger("biocypher").info(
            f"Found {len(epitopes) - len(unmatched_epitopes)} exact IEDB ID matches. "
            f"Trying substring matches for {len(unmatched_epitopes)} remaining epitopes..."
        )
        chunk_size = chunk_size // 2

        for i in range(0, len(unmatched_epitopes), chunk_size):
            chunk = unmatched_epitopes[i : i + chunk_size]
            substring_matches = _get_epitope_data(bc, chunk, base_url, match_type="substring")

            for epitope_seq in chunk:
                best_match = None
                min_length = float("inf")

                # Find matches where the epitope is a substring
                for match in substring_matches:
                    if epitope_seq in match["linear_sequence"] and len(match["linear_sequence"]) < min_length:
                        best_match = match
                        min_length = len(match["linear_sequence"])

                # Update the dictionary if a match was found
                if best_match:
                    antigens = best_match.get("curated_source_antigens")
                    if antigens:
                        antigen = best_match.get("curated_source_antigens")[0].get("name")
                        organism = best_match.get("curated_source_antigens")[0].get("source_organism_name")
                    else:
                        antigen = None
                        organism = None

                    epitope_to_iedb[epitope_seq] = {
                        "iri": f"iedb:{best_match.get('structure_id')}",
                        "antigen": antigen,
                        "organism": organism,
                    }

    # Final statistics
    matched_count = sum(1 for ep, info in epitope_to_iedb.items() if info["iri"].startswith("iedb:"))
    total_eps = len(epitopes)
    getLogger("biocypher").info(
        f"Epitope mapping results: {matched_count} of {total_eps} epitopes " f"matched to IEDB IDs ({matched_count / total_eps * 100:.1f}%)"
    )
    return epitope_to_iedb


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
        url = f"{base_url}?{check}&select=structure_id,structure_descriptions,linear_sequence,curated_source_antigens&order=structure_id"
        getLogger("biocypher").info(f"Request URL: {url[:100]}..." if len(url) > 100 else f"Request URL: {url}")

    else:
        request_hash = hashlib.md5("_".join(sorted(epitopes)).encode()).hexdigest()
        request_name = f"iedb_substring_matches{request_hash}"
        conditions = [f"linear_sequence.ilike.*{e}*" for e in epitopes]
        check = f"or=({','.join(conditions)})"
        url = f"{base_url}?{check}&select=structure_id,structure_descriptions,linear_sequence,curated_source_antigens&order=structure_id"

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
        getLogger("biocypher").info(f"API request failed: {e}")
        return []


def get_pmids_batch(bc: BioCypher, reference_urls: list[int], chunk_size: int = 150) -> dict[int, str]:
    """Retrieve PubMed IDs for multiple IEDB reference IDs using batched requests.

    Args:
        bc: BioCypher instance for the download
        reference_urls: List of IEDB reference URLs with IDs to query
        chunk_size: Size of chunks to break reference IDs into (to avoid URL length limits)

    Returns:
        Dictionary mapping IEDB reference IDs to their PubMed IDs (None if not found)
    """
    reference_ids_dic = {url: re.findall(r"\d+", url)[-1] for url in reference_urls if re.findall(r"\d+", url)}
    reference_ids = reference_ids_dic.values()

    base_url = "https://query-api.iedb.org/reference_export"
    reference_to_pmid = {}
    reference_ids = [ref_id for ref_id in reference_ids if ref_id is not None]

    getLogger("biocypher").info(f"Mapping {len(reference_ids)} IEDB reference IDs to PubMed IDs...")

    for i in range(0, len(reference_ids), chunk_size):
        chunk = reference_ids[i : i + chunk_size]
        reference_data = _get_reference_data(bc, chunk, base_url)

        # Initialize all reference IDs in chunk with None
        for ref_id in chunk:
            reference_to_pmid[ref_id] = f"no_pmid_{ref_id}"

        # Update with found matches
        for match in reference_data:
            ref_id = match.get("reference_id")
            pmid = match.get("reference__pmid")

            if str(ref_id) in reference_to_pmid and pmid is not None:
                # Only update if it's one we requested
                reference_to_pmid[str(ref_id)] = str(pmid)

    # Final statistics
    matched_count = sum(1 for ref_id, pmid in reference_to_pmid.items() if pmid is not None)
    total_refs = len(reference_ids)
    getLogger("biocypher").info(
        f"PMID mapping results: {matched_count} of {total_refs} reference IDs "
        f"matched to PubMed IDs ({matched_count / total_refs * 100:.1f}%)"
    )

    urls_to_pmids = {url: reference_to_pmid[id_val] for url, id_val in reference_ids_dic.items() if id_val in reference_to_pmid}

    return urls_to_pmids


def _get_reference_data(bc: BioCypher, reference_ids: list[int], base_url: str) -> list[dict]:
    """Get reference data for PubMed ID mapping.

    Args:
        bc: Biocypher instance for the download
        reference_ids: List of IEDB reference IDs to query
        base_url: Base URL for the API endpoint

    Returns:
        List of reference data dictionaries
    """
    request_hash = hashlib.md5("_".join(sorted(map(str, reference_ids))).encode()).hexdigest()
    request_name = f"iedb_reference_pmids_{request_hash}"

    reference_list = f"({','.join(map(str, reference_ids))})"
    check = f"reference_id=in.{reference_list}"
    url = f"{base_url}?{check}&select=reference_id,reference__pmid"

    getLogger("biocypher").info(f"Request URL: {url[:100]}..." if len(url) > 100 else f"Request URL: {url}")

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
        getLogger("biocypher").info(f"API request failed: {e}")
        return []


def save_airr_cells_json(airrcells: List[AirrCell], directory: str, filename: str = None, metadata: dict = None) -> None:
    """
    Save a list of AirrCell objects to a compressed JSON file with auto-generated filename.

    Args:
        airrcells (List[AirrCell]): List of AirrCell objects to save.
        directory (str): Directory path where to save the JSON file (e.g., "../data").
        filename (str, optional): Filename without extension.
        metadata (dict, optional): Metadata to include in the JSON file.
    """
    serialized_data = []

    for cell in airrcells:
        cell_data = {
            **dict(cell),
            "chains": cell.chains,
        }
        serialized_data.append(cell_data)

    output_data = serialized_data
    if metadata:
        output_data = {"metadata": metadata, "cells": serialized_data}

    # Generate filename with current date
    current_date = datetime.now().strftime("%Y%m%d%H%M%S")  # Format: YYYYMMDDHHMMSS
    if filename:
        filename = f"{filename}.json.gz"
    else:
        filename = f"airr_cells_{current_date}.json.gz"

    # Create full filepath
    filepath = os.path.join(directory, filename)

    # Create directory if it doesn't exist
    os.makedirs(directory, exist_ok=True)

    # Save as compressed JSON
    with gzip.open(filepath, "wt", encoding="utf-8") as f:
        json.dump(output_data, f, indent=2, ensure_ascii=False)

    getLogger("biocypher").info(f"Compressed JSON saved to: {filepath}")


def get_previous_release_metadata(repo_name: str = "iggytop/iggytop") -> dict | None:
    """
    Fetches metadata from the latest GitHub release of iggytop.
    """
    api_url = f"https://api.github.com/repos/{repo_name}/releases/latest"
    headers = {"Accept": "application/vnd.github.v3+json", "User-Agent": "iggytop-metadata-fetcher"}
    github_token = os.getenv("GITHUB_TOKEN")
    if github_token:
        headers["Authorization"] = f"token {github_token}"

    try:
        response = requests.get(api_url, headers=headers, timeout=10)
        if response.status_code == 200:
            release_data = response.json()
            # Try to find a metadata asset or parse it from the release body/assets
            # For now, we'll look for a specific 'metadata.json' asset if it exists
            for asset in release_data.get("assets", []):
                if asset["name"] == "metadata.json":
                    meta_url = asset["browser_download_url"]
                    meta_res = requests.get(meta_url, timeout=10)
                    if meta_res.status_code == 200:
                        return meta_res.json()
    except Exception:
        pass
    return None


def save_airr_cells_csv(airr_cells: List, directory: str) -> None:
    """
    Convert list of AirrCell objects to CSV format and save as compressed file.

    Args:
        airr_cells: List of AirrCell objects
        directory: Directory path where to save the CSV file (e.g., "../data")
    """
    rows = []

    for cell in airr_cells:
        row = {}

        # Add cell-level attributes (no prefix needed)
        # AirrCell implements MutableMapping, so we can iterate through it
        for key, value in cell.items():
            row[key] = value

        # Process chains using the chains property
        tra_chain = None
        trb_chain = None
        other_chains = []

        for chain in cell.chains:
            locus = chain.get("locus", "").upper()

            if locus == "TRA":
                tra_chain = chain
            elif locus == "TRB":
                trb_chain = chain
            else:
                other_chains.append((locus, chain))

        # Add TRA chain with chain_1_ prefix
        if tra_chain:
            for key, value in tra_chain.items():
                prefixed_key = f"chain_1_{key}"
                row[prefixed_key] = value

        # Add TRB chain with chain_2_ prefix
        if trb_chain:
            for key, value in trb_chain.items():
                prefixed_key = f"chain_2_{key}"
                row[prefixed_key] = value

        # Add other chains with their locus as prefix
        for locus, chain in other_chains:
            for key, value in chain.items():
                prefixed_key = f"chain_{locus.lower()}_{key}"
                row[prefixed_key] = value

        rows.append(row)

    # Create DataFrame
    df = pd.DataFrame(rows)

    # Generate filename with current date
    current_date = datetime.now().strftime("%d%m%Y")  # Format: DDMMYYYY
    filename = f"airr_cells_tabular_{current_date}.csv.gz"

    # Create full filepath
    filepath = os.path.join(directory, filename)

    # Create directory if it doesn't exist
    os.makedirs(directory, exist_ok=True)

    # Save to compressed CSV
    df.to_csv(filepath, index=False, compression="gzip")

    getLogger("biocypher").info(f"Compressed CSV saved to: {filepath}")
    getLogger("biocypher").info(f"Shape: {df.shape}")


def aggregate_unique_joined(series, separator="|"):
    """
    Helper function to aggregate unique values into a joined string.
    Warns if string 'nan' are found.
    """
    values = set()
    for v in series:
        if pd.isna(v) or str(v).lower() == "nan":
            continue

        s_v = str(v).strip()
        if s_v:
            values.update([s_v])
    if len(values) == 0:
        values.update(["nan"])
    return separator.join(sorted(values))


def deduplicate_and_aggregate(adata, subset_cols, agg_cols, separator="|"):
    """
    Deduplicates AnnData based on subset_cols and aggregates values in agg_cols.
    Uses scirpy airr_context to access TCR-specific columns if needed.
    """
    with ir.get.airr_context(adata, ["v_call", "junction_aa"], chain=["VJ_1", "VDJ_1"]) as m:
        obs_df = m.obs.copy()

        # Verify columns exist
        for col in subset_cols + agg_cols:
            if col not in obs_df.columns:
                raise KeyError(f"Column '{col}' not found in AnnData observations.")

        # Identify first occurrences based on subset_cols
        # Use drop_duplicates to find indices of records to keep
        # This is more robust than ~duplicated() when dealing with AIRR context/different views
        keep_indices = obs_df.reset_index().drop_duplicates(subset=subset_cols).set_index("cell_id").index

        # Define aggregation map
        agg_map = {col: (lambda s, sep=separator: aggregate_unique_joined(s, sep)) for col in agg_cols}

        # Group and aggregate
        # observed=True is used to avoid 'Product space too large' errors
        aggs = obs_df.groupby(subset_cols, observed=True, dropna=False, sort=False).agg(agg_map).reset_index()

        # Map aggregated info back to the deduplicated AnnData
        first_entries_keys = obs_df.loc[keep_indices, subset_cols].reset_index()
        final_info = first_entries_keys.merge(aggs, on=subset_cols, how="left")

    # Slice original AnnData using the identified indices
    deduplicated_adata = adata[keep_indices, :].copy()

    # Update aggregated columns
    # Ensure index alignment
    final_info.set_index("cell_id", inplace=True)
    for col in agg_cols:
        deduplicated_adata.obs[col] = final_info.loc[deduplicated_adata.obs.index, col]

    return deduplicated_adata
