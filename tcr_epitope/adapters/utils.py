import hashlib
import json
from datetime import datetime

import scanpy as sc
import scirpy.pp as scp
from biocypher import APIRequest, BioCypher
from scirpy.io import AirrCell, from_airr_cells

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def validate_peptide_sequence(seq: str) -> bool:
    """Checks if a given sequence is a valid peptide sequence."""
    if isinstance(seq, str):
        return all([aa in AMINO_ACIDS for aa in seq])
    else:
        return False


def get_iedb_ids_batch(epitopes: list[str], chunk_size: int = 50) -> dict[str, int]:
    """Retrieve IEDB IDs for multiple epitopes using batched requests.

    First tries exact matches, then falls back to substring matches for unmatched epitopes.

    Args:
        epitopes: List of epitope sequences to query
        chunk_size: Size of chunks to break epitopes into (to avoid URL length limits)

    Returns:
        Dictionary mapping epitope sequences to their IEDB IDs (0 if not found)
    """
    base_url = "https://query-api.iedb.org/epitope_search"
    epitope_to_id = {}
    unmatched_epitopes = []

    # Step 1: Try exact matches first
    print("Mapping AA epitope sequences to IEDB IDs: exact matches...")
    for i in range(0, len(epitopes), chunk_size):
        chunk = epitopes[i : i + chunk_size]
        # print(f"Processing batch {i//chunk_size + 1}/{(len(epitopes)-1)//chunk_size + 1}")
        epitope_matches = _get_epitope_data(chunk, base_url, match_type="exact")

        # Map results to the dictionary
        for epitope in chunk:
            epitope_to_id[epitope] = 0  # Default value if not found

        for match in epitope_matches:
            # Handle both possible API return formats for epitope sequences
            if match["structure_descriptions"]:
                epitope_seq = match["structure_descriptions"][0]
            elif "linear_sequence" in match:
                epitope_seq = match["linear_sequence"]
            else:
                continue

            if epitope_seq in epitope_to_id:  # Only update if it's one we requested
                epitope_to_id[epitope_seq] = match["structure_id"]

    # Step 2: Collect epitopes without matches
    unmatched_epitopes = [ep for ep, id_val in epitope_to_id.items() if id_val == 0]

    if unmatched_epitopes:
        print(
            f"Found {len(epitopes) - len(unmatched_epitopes)} exact IEDB ID matches. Trying substring matches for {len(unmatched_epitopes)} remaining epitopes..."
        )
        chunk_size = chunk_size // 2

        for i in range(0, len(unmatched_epitopes), chunk_size):
            chunk = unmatched_epitopes[i : i + chunk_size]
            # print(f"Processing unmatched batch {i//chunk_size + 1}/{(len(unmatched_epitopes)-1)//chunk_size + 1}")
            substring_matches = _get_epitope_data(chunk, base_url, match_type="substring")

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
                    epitope_to_id[epitope] = best_match["structure_id"]

    # Final statistics
    matched_count = sum(1 for id_val in epitope_to_id.values() if id_val > 0)
    print(
        f"Final results: {matched_count} of {len(epitopes)} epitopes matched to IEDB IDs ({matched_count / len(epitopes) * 100:.1f}%)"
    )

    return epitope_to_id


def _get_epitope_data(epitopes: list[str], base_url: str, match_type: str = "exact") -> list[dict]:
    """Get epitope data.

    Args:
        epitopes: List of epitope sequences to query
        base_url: Base URL for the API endpoint
        match_type: Type of matching to perform ("exact" or "substring")

    Returns:
        List of epitope data dictionaries
    """
    bc = BioCypher()

    if match_type == "exact":
        request_hash = hashlib.md5("_".join(sorted(epitopes)).encode()).hexdigest()
        request_name = f"iedb_exact_{request_hash}"
        epitope_list = f"({','.join([f'{e}' for e in epitopes])})"
        check = f"linear_sequence.in.{epitope_list}"
        url = f"{base_url}?{check}&select=structure_id,structure_descriptions,linear_sequence&order=structure_id"

    else:
        request_hash = hashlib.md5("_".join(sorted(epitopes)).encode()).hexdigest()
        request_name = f"iedb_substring_{request_hash}"
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


def kg_pd_to_anndata(df):
    """Convert KG data in pandas DataFrame format to an AnnData object.

    Handles cases where one chain may be missing or marked as NONE.

    Parameters
    ----------
    df : dict of DataFrames
        Dictionary containing filtered DataFrames for nodes and edges

    Returns:
    -------
    AnnData
        An AnnData object with TCR data formatted as AIRR cells
    """
    # Filter to get only the relevant nodes and edges
    tra_nodes = df["tra sequence"]
    trb_nodes = df["trb sequence"]
    epitope_nodes = df["epitope"]
    tcr_edges = df["alpha sequence to beta sequence association"]
    tcr_epitopes_edges = df["t cell receptor sequence to epitope association"]

    # Create dictionaries for faster lookups
    tra_dict = {row["node_id"]: row.to_dict() for _, row in tra_nodes.iterrows()}
    trb_dict = {row["node_id"]: row.to_dict() for _, row in trb_nodes.iterrows()}
    epitope_dict = {
        row["node_id"]: row.to_dict() for _, row in epitope_nodes.iterrows()
    }  # Keep it to add metadata to AnnData later

    # Initialize tracking structures
    tcr_pairs = {}  # Maps pair IDs to TRA and TRB IDs
    receptor_to_epitopes = {}  # Maps receptor IDs to epitope IDs

    # Extract TRA-TRB pairing information directly from the DataFrame
    for _, edge in tcr_edges.iterrows():
        relationship_id = edge["relationship_id"]
        tra_id = edge["source_id"]
        trb_id = edge["target_id"]

        # Store the pair
        tcr_pairs[relationship_id] = {"tra_id": tra_id, "trb_id": trb_id, "epitopes": set()}

        # Initialize epitope sets for each chain if not already done
        if tra_id not in receptor_to_epitopes:
            receptor_to_epitopes[tra_id] = set()
        if trb_id not in receptor_to_epitopes:
            receptor_to_epitopes[trb_id] = set()

    # Extract TCR-epitope associations
    for _, edge in tcr_epitopes_edges.iterrows():
        source_id = edge["source_id"]  # TCR chain (e.g., tra:CIRSGSARQLTF)
        target_id = edge["target_id"]  # Epitope (e.g., epitope:LLFGYPVYV)

        # Add to individual receptor's epitopes
        if source_id.startswith("tra:") or source_id.startswith("trb:"):
            if source_id not in receptor_to_epitopes:
                receptor_to_epitopes[source_id] = set()
            receptor_to_epitopes[source_id].add(target_id)

            # Find and add to all pairs containing this receptor
            for pair_id, pair_data in tcr_pairs.items():
                if source_id == pair_data["tra_id"] or source_id == pair_data["trb_id"]:
                    pair_data["epitopes"].add(target_id)

    # Process all receptors and create AIRR cells
    airr_cells = []
    processed_receptors = set()

    # Process paired receptors first
    for pair_id, pair_data in tcr_pairs.items():
        tra_id = pair_data["tra_id"]
        trb_id = pair_data["trb_id"]
        epitope_ids = pair_data["epitopes"]

        # Check if we should skip placeholder values
        is_tra_valid = tra_id in tra_dict and not tra_id.endswith(":NONE")
        is_trb_valid = trb_id in trb_dict and not trb_id.endswith(":NONE")

        # Skip if both chains are invalid
        if not (is_tra_valid or is_trb_valid):
            continue

        # Create a paired cell
        cell = AirrCell(cell_id=pair_id)

        # Mark chains as processed
        if is_tra_valid:
            processed_receptors.add(tra_id)
        if is_trb_valid:
            processed_receptors.add(trb_id)

        # Add TRA chain if valid
        if is_tra_valid:
            tra_data = tra_dict[tra_id]
            alpha_chain = AirrCell.empty_chain_dict()
            alpha_chain.update({"locus": "TRA", "junction_aa": extract_sequence_from_id(tra_id), "productive": True})

            # Add optional fields if available
            for field, key in [
                ("chain_1_v_gene", "v_call"),
                ("v_gene", "v_call"),
                ("chain_1_j_gene", "j_call"),
                ("j_gene", "j_call"),
                ("chain_1_organism", "species"),
                ("organism", "species"),
                ("species", "species"),
            ]:
                if tra_data.get(field):
                    alpha_chain[key] = tra_data[field]

            cell.add_chain(alpha_chain)

        # Add TRB chain if valid
        if is_trb_valid:
            trb_data = trb_dict[trb_id]
            beta_chain = AirrCell.empty_chain_dict()
            beta_chain.update({"locus": "TRB", "junction_aa": extract_sequence_from_id(trb_id), "productive": True})

            # Add optional fields if available
            for field, key in [
                ("chain_2_v_gene", "v_call"),
                ("v_gene", "v_call"),
                ("chain_2_d_gene", "d_call"),
                ("d_gene", "d_call"),
                ("chain_2_j_gene", "j_call"),
                ("j_gene", "j_call"),
                ("chain_2_organism", "species"),
                ("organism", "species"),
                ("species", "species"),
            ]:
                if trb_data.get(field):
                    beta_chain[key] = trb_data[field]

            cell.add_chain(beta_chain)

        # Add epitope metadata
        add_epitope_metadata(cell, epitope_ids)

        # Add cell metadata
        cell["data_source"] = "KG_DataFrames"
        cell["is_paired"] = is_tra_valid and is_trb_valid

        # Add source info if available
        if is_tra_valid and "source" in tra_data:
            cell["source_database"] = tra_data["source"]
        elif is_trb_valid and "source" in trb_data:
            cell["source_database"] = trb_data["source"]

        airr_cells.append(cell)

    # Process unpaired receptors with epitope associations
    # ... (code for unpaired receptors, similar to your original implementation)

    # Convert to AnnData
    if not airr_cells:
        return sc.AnnData()

    adata = from_airr_cells(airr_cells)
    scp.index_chains(adata)

    # Add metadata
    adata.uns["KG"] = {
        "source": "DataFrame_KG",
        "date_created": datetime.now().isoformat(),
        "cell_count": len(airr_cells),
        "paired_count": sum(1 for cell in airr_cells if cell.get("is_paired", True)),
        "unpaired_count": sum(1 for cell in airr_cells if not cell.get("is_paired", True)),
    }

    print("AnnData object created")
    print(adata)

    return adata


def extract_sequence_from_id(receptor_id):
    """Extract CDR3 sequence from receptor ID like 'tra:CIRSGSARQLTF'."""
    if ":" in receptor_id:
        seq = receptor_id.split(":", 1)[1]
        # Return None or empty string for placeholder values
        if seq == "NONE":
            return ""
        return seq
    return receptor_id


def add_epitope_metadata(cell, epitope_ids):
    """Helper function to add epitope metadata to an AirrCell."""
    epitope_sequences = []

    for epitope_id in epitope_ids:
        if epitope_id.startswith("epitope:"):
            # Extract the epitope sequence from the ID
            epitope_seq = extract_sequence_from_id(epitope_id)
            epitope_sequences.append(epitope_seq)

    # Store epitope data in cell metadata
    if epitope_sequences:
        if len(epitope_sequences) == 1:
            cell["antigen.epitope"] = epitope_sequences[0]
        else:
            # Use a string representation to avoid array issues
            cell["antigen.epitope"] = ", ".join(epitope_sequences)

    return cell
