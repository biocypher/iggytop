from __future__ import annotations

import re
from abc import abstractmethod
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

from .constants import REGISTRY_KEYS

if TYPE_CHECKING:
    import pandas as pd
    from biocypher import BioCypher


class BaseAdapter:
    """Base class for all adapters. This class is responsible for downloading and reading the data from the source.
    It also provides methods for generating BioCypher nodes and edges from the data.
    """

    def __init__(self, bc: BioCypher, cache_dir: str | None = None, test: bool = False):
        cache_dir = cache_dir or TemporaryDirectory().name
        table_path = self.get_latest_release(bc)
        self.table = self.read_table(bc, table_path, test)

    @abstractmethod
    def get_latest_release(self, bc: BioCypher, cache_dir: str) -> str:
        pass

    @abstractmethod
    def read_table(self, table_path: str, test: bool = False) -> pd.DataFrame:
        pass

    @abstractmethod
    def get_nodes(self):
        pass

    @abstractmethod
    def get_edges(self):
        pass

    def _generate_nodes_from_table(
        self,
        subset_cols: list[str],
        unique_cols: list[str] | None = None,
        property_cols: list[str] | None = None,
    ):
        if not isinstance(subset_cols, list):
            subset_cols = [subset_cols]

        unique_cols = unique_cols or subset_cols
        if not isinstance(unique_cols, list):
            unique_cols = [unique_cols]

        property_cols = property_cols or list(set(subset_cols) - set(unique_cols))
        if not isinstance(property_cols, list):
            property_cols = [property_cols]

        subset_table = self.table[subset_cols].dropna(subset=unique_cols)

        for _, row in subset_table.iterrows():
            if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in subset_cols:
                _type = row[REGISTRY_KEYS.CHAIN_1_TYPE_KEY]
            elif REGISTRY_KEYS.CHAIN_2_TYPE_KEY in subset_cols:
                _type = row[REGISTRY_KEYS.CHAIN_2_TYPE_KEY]
            else:
                _type = "epitope"

            # _id = ":".join([_type.lower(), *row[unique_cols].to_list()])
            
            # For TCR chains, use sequence + V gene + J gene as the identifier
            if _type.lower() != "epitope":
                # Get V gene and J gene if available
                v_gene_key = REGISTRY_KEYS.CHAIN_1_V_GENE_KEY if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in subset_cols else REGISTRY_KEYS.CHAIN_2_V_GENE_KEY
                j_gene_key = REGISTRY_KEYS.CHAIN_1_J_GENE_KEY if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in subset_cols else REGISTRY_KEYS.CHAIN_2_J_GENE_KEY
                
                # Check if V and J genes are available in the row
                v_gene = row.get(v_gene_key)
                j_gene = row.get(j_gene_key)
                
                # Create an ID that includes V and J genes if available
                id_components = [_type.lower()]
                id_components.extend(row[unique_cols].to_list())
                if v_gene:
                    id_components.append(f"{v_gene}")
                # if j_gene:
                    # id_components.append(f"j_{j_gene}")
                
                _id = ":".join(id_components)
            else:
                # For epitopes and other types, keep the original ID format
                _id = ":".join([_type.lower(), *row[unique_cols].to_list()])
            
            _props = {re.sub("chain_\d_", "", k): row[k] for k in property_cols}
            # _props["junction_aa"] = row[unique_cols[0]] if unique_cols else None

            yield _id, _type.lower(), _props

    def _generate_edges_from_table(
        self,
        source_subset_cols: list[str],
        target_subset_cols: list[str],
        source_unique_cols: list[str] | None = None,
        target_unique_cols: list[str] | None = None,
    ):
        source_subset_cols = source_subset_cols or []
        if not isinstance(source_subset_cols, list):
            source_subset_cols = [source_subset_cols]

        source_unique_cols = source_unique_cols or source_subset_cols
        if not isinstance(source_unique_cols, list):
            source_unique_cols = [source_unique_cols]

        target_subset_cols = target_subset_cols or []
        if not isinstance(target_subset_cols, list):
            target_subset_cols = [target_subset_cols]

        target_unique_cols = target_unique_cols or target_subset_cols
        if not isinstance(target_unique_cols, list):
            target_unique_cols = [target_unique_cols]

        subset_table = (
            self.table[source_subset_cols + target_subset_cols]
            .drop_duplicates(subset=source_unique_cols + target_unique_cols)
            .dropna(subset=source_unique_cols + target_unique_cols)
        )

        for _, row in subset_table.iterrows():

            node_data = {}
            for i in ["source", "target"]:
                cols = locals()[f"{i}_subset_cols"]
                if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in cols:
                    node_type = row[REGISTRY_KEYS.CHAIN_1_TYPE_KEY]
                    v_gene_key = REGISTRY_KEYS.CHAIN_1_V_GENE_KEY
                elif REGISTRY_KEYS.CHAIN_2_TYPE_KEY in cols:
                    node_type = row[REGISTRY_KEYS.CHAIN_2_TYPE_KEY]
                    v_gene_key = REGISTRY_KEYS.CHAIN_2_V_GENE_KEY
                else:
                    node_type = "epitope"
                    v_gene_key = None

                id_components = [node_type.lower()]
                id_components.extend(row[locals()[f"{i}_unique_cols"]].tolist())

                if v_gene_key:
                    v_gene = row[v_gene_key]
                    if v_gene:
                        id_components.append(v_gene)

                node_data[i] = {
                    "id": ":".join(id_components),
                    "type": node_type
                }

            _source_id = node_data["source"]["id"]
            _target_id = node_data["target"]["id"]
            _source_type = node_data["source"]["type"]
            _target_type = node_data["target"]["type"]

            _id = f"{_source_id}-{_target_id}"
            _type = f"{_source_type.lower()}_to_{_target_type.lower()}"

            yield (_id, _source_id, _target_id, _type, {})