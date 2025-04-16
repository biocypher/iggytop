from __future__ import annotations

from abc import abstractmethod
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

from .constants import REGISTRY_KEYS

if TYPE_CHECKING:
    import pandas as pd
    from biocypher import BioCypher


class BaseAdapter:
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

        subset_table = self.table[subset_cols].drop_duplicates(subset=unique_cols).dropna(subset=unique_cols)
        for _, row in subset_table.iterrows():
            if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in subset_cols:
                _type = row[REGISTRY_KEYS.CHAIN_1_TYPE_KEY]
            elif REGISTRY_KEYS.CHAIN_2_TYPE_KEY in subset_cols:
                _type = row[REGISTRY_KEYS.CHAIN_2_TYPE_KEY]
            else:
                _type = "epitope"

            _id = ":".join([_type.lower()] + row[unique_cols].to_list())
            _props = {k: row[k] for k in property_cols}

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
            if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in source_subset_cols:
                _source_type = row[REGISTRY_KEYS.CHAIN_1_TYPE_KEY]
            elif REGISTRY_KEYS.CHAIN_2_TYPE_KEY in source_subset_cols:
                _source_type = row[REGISTRY_KEYS.CHAIN_2_TYPE_KEY]
            else:
                _source_type = "epitope"

            if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in target_subset_cols:
                _target_type = row[REGISTRY_KEYS.CHAIN_1_TYPE_KEY]
            elif REGISTRY_KEYS.CHAIN_2_TYPE_KEY in target_subset_cols:
                _target_type = row[REGISTRY_KEYS.CHAIN_2_TYPE_KEY]
            else:
                _target_type = "epitope"

            _source_id = ":".join([_source_type.lower()] + row[source_unique_cols].to_list())
            _target_id = ":".join([_target_type.lower()] + row[target_unique_cols].to_list())
            _id = "-".join([_source_id, _target_id])
            _type = "_".join([_source_type.lower(), "to", _target_type.lower()])

            yield (_id, _source_id, _target_id, _type, {})
