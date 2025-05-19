import hashlib
import json

from biocypher import APIRequest, BioCypher

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def validate_peptide_sequence(seq: str) -> bool:
    """Checks if a given sequence is a valid peptide sequence."""
    if isinstance(seq, str):
        return all([aa in AMINO_ACIDS for aa in seq])
    else:
        return False


def process_sequence(x):
    if x is None:
        return None

    x = str(x)
    result = x.split("+")[0]  # split_epitope_sequence
    result = result.upper()
    result = "".join(result.split())
    # TODO: add here sequence validation but avoid the non-peptidic epitopes in IEDB db!
    return result


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

    # Step 1: Try exact matches first
    print("Mapping AA epitope sequences to IEDB IDs: exact matches...")

    for i in range(0, len(epitopes), chunk_size):
        chunk = epitopes[i : i + chunk_size]
        # print("Starting new batch...")

        epitope_matches = _get_epitope_data(bc, chunk, base_url, match_type="exact")

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
            f"Found {len(epitopes) - len(unmatched_epitopes)} exact IEDB ID matches.",
            f"Trying substring matches for {len(unmatched_epitopes)} remaining epitopes..."
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
                    epitope_to_id[epitope] = best_match["structure_id"]

    # Final statistics
    matched_count = sum(1 for id_val in epitope_to_id.values() if id_val > 0)
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
    