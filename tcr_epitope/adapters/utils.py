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
        print(f"Processing batch {i // chunk_size + 1}/{(len(epitopes) - 1) // chunk_size + 1}")
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
            f"Found {len(epitopes) - len(unmatched_epitopes)} exact IEDB ID matches. Trying substring matches for {len(unmatched_epitopes)} remaining epitopes..."
        )
        chunk_size = chunk_size // 2

        for i in range(0, len(unmatched_epitopes), chunk_size):
            chunk = unmatched_epitopes[i : i + chunk_size]
            # print(f"Processing unmatched batch {i//chunk_size + 1}/{(len(unmatched_epitopes)-1)//chunk_size + 1}")
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


# TODO: Temporary stored here for testing purposes before introduced in Biocypher
from datetime import datetime

import scanpy as sc
import scirpy.pp as scp
from biocypher._create import BioCypherEdge, BioCypherNode, BioCypherRelAsNode
from biocypher._deduplicate import Deduplicator
from biocypher.output.in_memory._in_memory_kg import _InMemoryKG
from scirpy.io import AirrCell, from_airr_cells


class AnnDataKG(_InMemoryKG):
    def __init__(self, deduplicator=None):
        super().__init__()
        self.deduplicator = deduplicator or Deduplicator()
        # Store entities by type directly
        self.entities_by_type = {}

    def get_kg(self):
        """Convert directly to AnnData instead of going through DataFrames"""
        return self.to_anndata()

    def add_nodes(self, nodes):
        """Add BioCypher nodes, organizing them by type."""
        lists = self._separate_entity_types(nodes)
        self._add_to_entities_by_type(lists)

    def add_edges(self, edges):
        """Add BioCypher edges, organizing them by type."""
        lists = self._separate_entity_types(edges)
        self._add_to_entities_by_type(lists)

    def _add_to_entities_by_type(self, lists):
        """Add entities from lists to the entities_by_type dictionary."""
        for _type, _entities in lists.items():
            if _type not in self.entities_by_type:
                self.entities_by_type[_type] = []
            self.entities_by_type[_type].extend(_entities)

    def _separate_entity_types(self, entities):
        """Given mixed iterable of BioCypher objects, separate them into lists by
        type. Also deduplicates using the `Deduplicator` instance.
        """
        lists = {}
        count = 0
        total = len(entities) if hasattr(entities, "__len__") else 0

        # Set up progress reporting to show progress without cluttering output
        report_interval = max(1, min(1000, total // 5)) if total > 0 else 0

        for entity in entities:
            count += 1

            # Report progress at reasonable intervals
            if report_interval > 0 and (count % report_interval == 0 or count == total):
                print(f"Processing entities: {count}/{total} ({count / total * 100:.1f}%)")

            # Type checking
            if (
                not isinstance(entity, BioCypherNode)
                and not isinstance(entity, BioCypherEdge)
                and not isinstance(entity, BioCypherRelAsNode)
            ):
                raise TypeError(f"Expected a BioCypherNode / BioCypherEdge / BioCypherRelAsNode, got {type(entity)}.")

            # Deduplication check
            if isinstance(entity, BioCypherNode):
                seen = self.deduplicator.node_seen(entity)
            elif isinstance(entity, BioCypherEdge):
                seen = self.deduplicator.edge_seen(entity)
            elif isinstance(entity, BioCypherRelAsNode):
                seen = self.deduplicator.rel_as_node_seen(entity)

            if seen:
                continue

            if isinstance(entity, BioCypherRelAsNode):
                node = entity.get_node()
                source_edge = entity.get_source_edge()
                target_edge = entity.get_target_edge()

                _type = node.get_type()
                if _type not in lists:
                    lists[_type] = []
                lists[_type].append(node)

                _source_type = source_edge.get_type()
                if _source_type not in lists:
                    lists[_source_type] = []
                lists[_source_type].append(source_edge)

                _target_type = target_edge.get_type()
                if _target_type not in lists:
                    lists[_target_type] = []
                lists[_target_type].append(target_edge)
                continue

            # Regular node or edge
            _type = entity.get_type()

            if _type not in lists:
                lists[_type] = []
            lists[_type].append(entity)

        # Final progress report
        if total > 0:
            print(f"Processed {count} entities into {len(lists)} types")
        return lists

    def to_anndata(self, verbose=False):
        """Directly convert stored BioCypher entities to AnnData."""
        # Get required entities
        tra_nodes = self.entities_by_type.get("tra sequence", [])
        trb_nodes = self.entities_by_type.get("trb sequence", [])
        epitope_nodes = self.entities_by_type.get("epitope", [])
        tcr_edges = self.entities_by_type.get("alpha sequence to beta sequence association", [])
        tcr_epitope_edges = self.entities_by_type.get("t cell receptor sequence to epitope association", [])

        # Create efficient lookups directly from BioCypher objects
        tra_dict = {node.get_id(): node.get_properties() for node in tra_nodes}
        trb_dict = {node.get_id(): node.get_properties() for node in trb_nodes}
        epitope_dict = {node.get_id(): node.get_properties() for node in epitope_nodes}

        tcr_pairs = {}  # Maps relationship_id to pair data
        receptor_to_epitopes = {}  # Maps receptor IDs to epitope sets
        processed_receptors = set()

        # Process epitope associations
        print("Processing epitope associations")
        for edge in tcr_epitope_edges:
            edge_dict = edge.get_dict()
            source_id = edge_dict["source_id"]  # TCR chain
            target_id = edge_dict["target_id"]  # Epitope

            if source_id not in receptor_to_epitopes:
                receptor_to_epitopes[source_id] = set()
            receptor_to_epitopes[source_id].add(target_id)
            # Process paired TCRs first

        print("Processing paired TCRs")
        for edge in tcr_edges:
            edge_dict = edge.get_dict()
            relationship_id = edge_dict["relationship_id"]
            tra_id = edge_dict["source_id"]
            trb_id = edge_dict["target_id"]

            tcr_pairs[relationship_id] = {"tra_id": tra_id, "trb_id": trb_id, "epitopes": set()}

            processed_receptors.add(tra_id)
            processed_receptors.add(trb_id)

            # Initialize epitope tracking for receptors
            if tra_id in receptor_to_epitopes:
                tcr_pairs[relationship_id]["epitopes"].update(receptor_to_epitopes[tra_id])
            if trb_id in receptor_to_epitopes:
                tcr_pairs[relationship_id]["epitopes"].update(receptor_to_epitopes[trb_id])

        # Create AIRR cells
        airr_cells = []

        for pair_id, pair_data in tcr_pairs.items():
            tra_id = pair_data["tra_id"]
            trb_id = pair_data["trb_id"]
            epitope_ids = pair_data["epitopes"]

            # Create cell
            cell = AirrCell(cell_id=pair_id)

            # Add chains
            if tra_id in tra_dict:
                tra_data = tra_dict[tra_id]
                alpha_chain = AirrCell.empty_chain_dict()
                alpha_chain.update(
                    {
                        "locus": "TRA",
                        "junction_aa": extract_sequence_from_id(tra_id),
                        "v_call": tra_data["chain_1_v_gene"],
                        "j_call": tra_data["chain_1_j_gene"],
                        "consensus_count": 0,
                        "productive": True,
                    }
                )
                cell.add_chain(alpha_chain)

            if trb_id in trb_dict:
                trb_data = trb_dict[trb_id]
                beta_chain = AirrCell.empty_chain_dict()
                beta_chain.update(
                    {
                        "locus": "TRB",
                        "junction_aa": extract_sequence_from_id(trb_id),
                        "v_call": trb_data["chain_2_v_gene"],
                        "j_call": trb_data["chain_2_j_gene"],
                        "consensus_count": 0,
                        "productive": True,
                    }
                )
                cell.add_chain(beta_chain)

            # Add epitope metadata
            if epitope_ids:
                add_epitope_metadata(epitope_dict, cell, epitope_ids)

            # Add cell metadata
            cell["data_source"] = "BioCypher"
            cell["is_paired"] = True
            airr_cells.append(cell)

        print("Processing unpaired TCRs")
        unpaired_count = 0

        # Single pass through receptor_to_epitopes (more efficient)
        for receptor_id, epitope_ids in receptor_to_epitopes.items():
            # Skip if already in a pair or no epitopes
            if receptor_id in processed_receptors or not epitope_ids:
                continue

            # Create cell
            cell_id = f"unpaired_{receptor_id}"
            cell = AirrCell(cell_id=cell_id)

            # Add chains
            # TODO: isolate chain addition in a separate function
            if receptor_id.startswith("tra:"):
                tra_data = tra_dict[receptor_id]
                alpha_chain = AirrCell.empty_chain_dict()
                alpha_chain.update(
                    {
                        "locus": "TRA",
                        "junction_aa": extract_sequence_from_id(receptor_id),
                        "v_call": tra_data["chain_1_v_gene"],
                        "j_call": tra_data["chain_1_j_gene"],
                        "consensus_count": 0,
                        "productive": True,
                    }
                )
                cell.add_chain(alpha_chain)

            if receptor_id.startswith("trb:"):
                trb_data = trb_dict[trb_id]
                beta_chain = AirrCell.empty_chain_dict()
                beta_chain.update(
                    {
                        "locus": "TRB",
                        "junction_aa": extract_sequence_from_id(receptor_id),
                        "v_call": trb_data["chain_2_v_gene"],
                        "j_call": trb_data["chain_2_j_gene"],
                        "consensus_count": 0,
                        "productive": True,
                    }
                )
                cell.add_chain(beta_chain)

            # Add epitope metadata
            if epitope_ids:
                add_epitope_metadata(epitope_dict, cell, epitope_ids)

            # Add cell metadata
            cell["data_source"] = "BioCypher"
            cell["is_paired"] = False
            airr_cells.append(cell)

            unpaired_count += 1

        if verbose and unpaired_count > 0:
            print(f"Added {unpaired_count} unpaired receptors with epitope associations")

        # Convert to AnnData
        if not airr_cells:
            return sc.AnnData()

        adata = from_airr_cells(airr_cells)
        scp.index_chains(adata)

        # Add metadata
        adata.uns["DB"] = {
            "name": "BioCypher_KG",
            "date_created": datetime.now().isoformat(),
            "cell_count": len(airr_cells),
            "paired_count": sum(1 for cell in airr_cells if cell.get("is_paired", True)),
            "unpaired_count": sum(1 for cell in airr_cells if not cell.get("is_paired", False)),
        }

        if verbose:
            print(f"Created AnnData object with {len(airr_cells)} cells")
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


def add_epitope_metadata(epitope_dict, cell, epitope_ids):
    """Helper function to add epitope metadata to an AirrCell."""
    # TODO: now I only add the data about one epitope to the cell, figure how to store several
    for epitope_id in epitope_ids:
        if epitope_id in epitope_dict:
            if "properties" in epitope_dict[epitope_id]:
                epitope_data = epitope_dict[epitope_id]["properties"]
            else:
                epitope_data = epitope_dict[epitope_id]

            keys_to_remove = ["node_id", "node_label", "id", "preferred_id"]
            epitope_data = {key: value for key, value in epitope_data.items() if key not in keys_to_remove}

    for prop_key, prop_value in epitope_data.items():
        cell[prop_key] = prop_value

    return cell
