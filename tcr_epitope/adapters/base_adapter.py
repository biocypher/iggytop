from abc import abstractmethod
from tempfile import TemporaryDirectory
from typing import Optional

import pandas as pd

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
