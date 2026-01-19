from __future__ import annotations

import re
import os
import pandas as pd
from pathlib import Path
from abc import abstractmethod
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING
from datetime import datetime
from typing import cast
from tqdm.auto import tqdm
from scirpy.io._datastructures import AirrCell
from scirpy.io._convert_anndata import from_airr_cells
from scirpy.pp import index_chains
from .constants import REGISTRY_KEYS

if TYPE_CHECKING:
    import pandas as pd
    from biocypher import BioCypher


class BaseAdapter:
    """
    Base class for all adapters.

    This class is responsible for the basic structure and function for Iggytop adapters.
    It initializes any adapter by calling the corresponding function for downloading and reading the data from the source.
    It also provides methods for generating BioCypher nodes and edges from the data.

    Attributes:
        table (pd.DataFrame): The data table read from the source.
    
    Args:
        bc (BioCypher): An instance of the BioCypher class.
        cache_dir (str | None, optional): Directory to cache data. Defaults to None.
        test (bool, optional): Whether to run in test mode. Defaults to False.
    """

    def __init__(self, bc: BioCypher, cache_dir: str | None = None, test: bool = False):
        """
        Initializes the BaseAdapter instance.

        Args:
            bc (BioCypher): An instance of the BioCypher class.
            cache_dir (str | None, optional): Directory to cache data. Defaults to None.
            test (bool, optional): Whether to run in test mode. Defaults to False.
        """
        self.cache_dir = cache_dir or TemporaryDirectory().name
        table_path = self.get_latest_release(bc)
        self.table = self.read_table(bc, table_path, test)
        self.airr_cells = None

        if not hasattr(self.__class__, "DB_NAME"):
            raise TypeError(f"Class {self.__class__.__name__} must define a 'DB_NAME' class attribute.")
        
    @abstractmethod
    def get_latest_release(self, bc: BioCypher) -> str:
        """
        Abstract method to get the latest release of the data.

        Args:
            bc (BioCypher): An instance of the BioCypher class.
            cache_dir (str): Directory to cache data.

        Returns:
            str: Path to the latest release.
        """
        pass

    @abstractmethod
    def read_table(self, bc: BioCypher, table_path: str, test: bool = False) -> pd.DataFrame:
        """
        Abstract method to read and harmonize the data table from the source.

        Args:
            table_path (str): Path to the data table.
            test (bool, optional): Whether to run in test mode. Defaults to False.

        Returns:
            pd.DataFrame: The data table.
        """
        pass

    @abstractmethod
    def get_nodes(self):
        """
        Abstract method to generate BioCypher nodes from the data. 
        This method is intended to use _generate_nodes_from_table with the right parameters for each edge type. This requires parameters depending on the adapter used.

        Returns:
            Iterable: An iterable of BioCypher nodes.
        """
        pass

    @abstractmethod
    def get_edges(self):
        """
        Abstract method to generate BioCypher edges from the data.
        This method is intended to call _generate_edges_from_table with the right parameters for each edge type. This requires parameters depending on the adapter used.

        Returns:
            Iterable: An iterable of BioCypher edges.
        """
        pass

    def create_airr_cells(self) -> bool:
        """
        Converts the table data to AIRR cell format. The list of AIRR cells is stored in self.airr_cells.

        Returns:
            bool: True if successful, False otherwise.

        """
        if self.airr_cells is not None:
            return True
        
        df = self.table
        self.airr_cells = []

        for idx, row in tqdm(df.iterrows(), total=df.shape[0], desc=f"Processing {self.DB_NAME} entries"):
            cell = AirrCell(cell_id=str(idx))
            if not pd.isnull(row[REGISTRY_KEYS.CHAIN_1_CDR3_KEY]):
                alpha_chain = AirrCell.empty_chain_dict()
                alpha_chain.update(
                    {
                        "locus": row[REGISTRY_KEYS.CHAIN_1_TYPE_KEY],
                        "junction_aa": row[REGISTRY_KEYS.CHAIN_1_CDR3_KEY],
                        "v_call": row[REGISTRY_KEYS.CHAIN_1_V_GENE_KEY],
                        "j_call": row[REGISTRY_KEYS.CHAIN_1_J_GENE_KEY],
                        "consensus_count": 0,
                        "productive": True,
                    }
                )
                cell.add_chain(alpha_chain)

            if not pd.isnull(row[REGISTRY_KEYS.CHAIN_2_CDR3_KEY]):
                beta_chain = AirrCell.empty_chain_dict()
                beta_chain.update(
                    {
                        "locus": row[REGISTRY_KEYS.CHAIN_2_TYPE_KEY],
                        "junction_aa": row[REGISTRY_KEYS.CHAIN_2_CDR3_KEY],
                        "v_call": row[REGISTRY_KEYS.CHAIN_2_V_GENE_KEY],
                        "j_call": row[REGISTRY_KEYS.CHAIN_2_J_GENE_KEY],
                        "consensus_count": 0,
                        "productive": True,
                    }
                )
                cell.add_chain(beta_chain)

            for f in REGISTRY_KEYS:
                if f not in row:
                    continue
                cell[f] = row[f]
            self.airr_cells.append(cell)
        return True

    def create_anndata(self) -> None:
        """
        Creates an Anndata object from the AIRR cell data and saves it to a file in the cache directory.
        """
        self.create_airr_cells()   
        adata = from_airr_cells(self.airr_cells)
        index_chains(adata)

        # Convert object columns to string to avoid serialization issues with h5py (e.g. for PMID)
        for col in adata.obs.columns:
            if adata.obs[col].dtype == object:
                adata.obs[col] = adata.obs[col].astype(str)

        adata.uns["DB"] = {"name": self.DB_NAME, "date_downloaded": datetime.now().isoformat()}
        anndata_path = Path(self.cache_dir) / f"{self.DB_NAME}_anndata.h5ad"
        anndata_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(cast(os.PathLike, anndata_path))
        print(f"Saved Anndata to {anndata_path}")
    

    def _generate_nodes_from_table(
        self,
        subset_cols: list[str],
        unique_cols: list[str] | None = None,
        property_cols: list[str] | None = None,
    ):
        """
        Generates BioCypher nodes from the data table.

        The unique_cols are used for selecting the rows which contain relevant information.
        They do NOT correspond to the unique identifier.
        To create the unique identifier, we use unique_cols + V gene (if available) for TCR chains.

        Args:
            subset_cols (list[str]): List of columns to subset the table.
            unique_cols (list[str] | None, optional): List of columns to check for uniqueness. Defaults to None.
            property_cols (list[str] | None, optional): List of columns to include as properties. Defaults to None.

        Yields:
            tuple: A tuple containing the node ID, node type, and properties.
        """
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
        """
        Generates BioCypher edges from the data table.

        The unique_cols are used for selecting the rows which contain relevant information.
        They do NOT correspond to the unique identifier.
        To create the unique identifier, we use unique_cols + V gene (if available) for TCR chains.

        Args:
            source_subset_cols (list[str]): List of columns for the source node.
            target_subset_cols (list[str]): List of columns for the target node.
            source_unique_cols (list[str] | None, optional): List of unique columns for the source node. Defaults to None.
            target_unique_cols (list[str] | None, optional): List of unique columns for the target node. Defaults to None.

        Yields:
            tuple: A tuple containing the edge ID, source ID, target ID, edge type, and properties.
        """
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