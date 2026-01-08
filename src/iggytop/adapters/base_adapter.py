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
        cache_dir = cache_dir or TemporaryDirectory().name
        table_path = self.get_latest_release(bc)
        self.table = self.read_table(bc, table_path, test)

    @abstractmethod
    def get_latest_release(self, bc: BioCypher, cache_dir: str) -> str:
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
    def read_table(self, table_path: str, test: bool = False) -> pd.DataFrame:
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

        subset_table = (
            self.table[list(set(subset_cols).union(property_cols or []))]
            .dropna(subset=unique_cols)
        )
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
        exclude_cols: list[str] | None = None,
        target_unique_cols: list[str] | None = None,
        property_cols: list[str] | None = None,
    ):
        """
        Generates BioCypher edges from the data table.

        The unique_cols are used for selecting the rows which contain relevant information.
        They do NOT correspond to the unique identifier.
        To create the unique identifier, we use unique_cols + V gene (if available) for TCR chains.

        Args:
            source_subset_cols (list[str]): List of columns for the source node. Includes data used of the identifier.
            target_subset_cols (list[str]): List of columns for the target node. Includes data used of the identifier.
            source_unique_cols (list[str] | None, optional): List of unique columns for the source node. Used to filter rows with non-missing values. Not necessarily identical to the identifier columns, eg V gene is included in the id but only if it is available.
            exclude_cols (list[str] | None, optional): List of columns which must be NaN. Can be used to aviod duplicating edges which are represented in other edge types.
            target_unique_cols (list[str] | None, optional): List of unique columns for the target node. Used to filter rows with non-missing values. Not necessarily identical to the identifier columns, eg V gene is included in the id but only if it is available.

        Yields:
            tuple: A tuple containing the edge ID, source ID, target ID, edge type, and properties.
        """
        source_subset_cols = source_subset_cols or []
        if not isinstance(source_subset_cols, list):
            source_subset_cols = [source_subset_cols]

        source_unique_cols = source_unique_cols or source_subset_cols
        if not isinstance(source_unique_cols, list):
            source_unique_cols = [source_unique_cols]

        if exclude_cols and not isinstance(exclude_cols, list):
            exclude_cols = [exclude_cols]

        target_subset_cols = target_subset_cols or []
        if not isinstance(target_subset_cols, list):
            target_subset_cols = [target_subset_cols]

        target_unique_cols = target_unique_cols or target_subset_cols
        if not isinstance(target_unique_cols, list):
            target_unique_cols = [target_unique_cols]

        if exclude_cols:
            subset_table = self.table[self.table[exclude_cols].isna().all(axis=1)]
        else:
            subset_table = self.table
        subset_table = (
            subset_table[list(set(source_subset_cols + target_subset_cols + (property_cols or [])))]
            .dropna(subset=source_unique_cols + target_unique_cols)
        )

        for _, row in subset_table.iterrows():

            node_data = {}
            for i in ["source", "target"]:
                cols = locals()[f"{i}_subset_cols"]


                if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in cols and REGISTRY_KEYS.CHAIN_2_TYPE_KEY in cols:
                    chain1 = row[REGISTRY_KEYS.CHAIN_1_TYPE_KEY]+":"+row[locals()[f"{i}_unique_cols"][0]] 
                    if  row[REGISTRY_KEYS.CHAIN_1_V_GENE_KEY]:
                        chain1 += ":" + row[REGISTRY_KEYS.CHAIN_1_V_GENE_KEY]
                    chain2 = row[REGISTRY_KEYS.CHAIN_2_TYPE_KEY]+":"+ row[locals()[f"{i}_unique_cols"][1]]
                    if  row[REGISTRY_KEYS.CHAIN_2_V_GENE_KEY]:
                        chain2 += ":" + row[REGISTRY_KEYS.CHAIN_2_V_GENE_KEY]
                        
                    node_data[i] = {
                        "id": f"pair:{chain1}-{chain2}",
                        "type": row[REGISTRY_KEYS.CHAIN_1_TYPE_KEY].lower()+"_"+row[REGISTRY_KEYS.CHAIN_2_TYPE_KEY].lower()+"_pair"
                    }
                else:                    
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
            if _type == "tra_to_trb" or _type == "igh_to_igl": #if the edge defines a pair
                _id = f"pair:{_source_id}-{_target_id}"
                _type =  row[REGISTRY_KEYS.CHAIN_1_TYPE_KEY].lower()+"_"+row[REGISTRY_KEYS.CHAIN_2_TYPE_KEY].lower()+"_pair"
            if property_cols:
                _props = {re.sub("chain_\d_", "", k): row[k] for k in property_cols}
            else:
                _props = {}
            # _props["junction_aa"] = row[unique_cols[0]] if unique_cols else None

            yield (_id, _source_id, _target_id, _type, _props)