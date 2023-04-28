from abc import abstractmethod
from tempfile import TemporaryDirectory
from typing import List, Optional

import pandas as pd

from .constants import REGISTRY_KEYS


class BaseAdapter:

    def __init__(self, cache_dir: Optional[str] = None, test: bool = False):
        cache_dir = cache_dir or TemporaryDirectory().name
        table_path = self.get_latest_release(cache_dir)
        self.table = self.read_table(table_path, test=test)

    @abstractmethod
    def get_latest_release(self, save_dir: str) -> str:
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
        subset_cols: List[str], 
        unique_cols: Optional[List[str]] = None,
        property_cols: Optional[List[str]] = None,
    ):
        if not isinstance(subset_cols, list):
            subset_cols = [subset_cols]

        unique_cols = unique_cols or subset_cols
        if not isinstance(unique_cols, list):
            unique_cols = [unique_cols]

        property_cols = property_cols or list(set(subset_cols) - set(unique_cols))
        if not isinstance(property_cols, list):
            property_cols = [property_cols]

        subset_table = self.table[subset_cols].drop_duplicates(
            subset=unique_cols
        ).dropna(subset=unique_cols)
        for _, row in subset_table.iterrows():
            if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in subset_cols:
                _type = row[REGISTRY_KEYS.CHAIN_1_TYPE_KEY]
            elif REGISTRY_KEYS.CHAIN_2_TYPE_KEY in subset_cols:
                _type = row[REGISTRY_KEYS.CHAIN_2_TYPE_KEY]
            else:
                _type = "epitope"

            _id = "_".join([_type] + row[unique_cols].to_list())
            _props = {k: row[k] for k in property_cols}

            yield _id, _type, _props
    
    def _generate_edges_from_table(
        self,
        source_subset_cols: List[str], 
        target_subset_cols: List[str],
        source_unique_cols: Optional[List[str]] = None,
        target_unique_cols: Optional[List[str]] = None,
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

        subset_table = self.table[source_subset_cols + target_subset_cols].drop_duplicates(
            subset=source_unique_cols + target_unique_cols
        ).dropna(subset=source_unique_cols + target_unique_cols)

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

            _source_id = "_".join([_source_type] + row[source_unique_cols].to_list())
            _target_id = "_".join([_target_type] + row[target_unique_cols].to_list())
            _id = "-".join([_source_id, _target_id])
            _type = "_".join([_source_type, "to", _target_type])

            yield (_id, _source_id, _target_id, _type, {})
