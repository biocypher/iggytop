from __future__ import annotations

import os
import re
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Sequence, cast

import pandas as pd
from scirpy.io._convert_anndata import from_airr_cells
from scirpy.io._datastructures import AirrCell
from scirpy.pp import index_chains
from tqdm.auto import tqdm

from .constants import REGISTRY_KEYS
from .utils import get_file_checksum

if TYPE_CHECKING:
    from biocypher import BioCypher
import platformdirs


class BaseAdapter(ABC):
    """
    Base class for all adapters.

    This class is responsible for the basic structure and function for Iggytop adapters.
    It initializes any adapter by calling the corresponding function for downloading and reading the data from the source.
    It also provides methods for generating BioCypher nodes and edges from the data.

    Attributes:
        table (pd.DataFrame): The data table read from the source.
        DB_NAME (str): Name of the database. Must be defined in subclasses.
        available_receptors (list[str]): List of receptor types available in the database. Must be defined in subclasses.

    Args:
        bc: An instance of the BioCypher class.
        cache_dir: Directory to cache data. Defaults to None.
        receptors_to_include: Receptors to include. Defaults to ("TCR", "BCR").
        test: Whether to run in test mode. Defaults to False.
        filter_10x: Whether to filter out 10X Genomics datasets. Defaults to False.
    """

    # These should be defined in child classes
    DB_NAME: str
    available_receptors: list[str]

    def __init__(
        self,
        bc: BioCypher,
        cache_dir: str | None = None,
        receptors_to_include: Sequence[Literal["TCR", "BCR"]] | None = ("TCR", "BCR"),
        test: bool = False,
        filter_10x: bool = False,
    ):
        """
        Initializes the BaseAdapter instance.

        Args:
            bc: An instance of the BioCypher class.
            cache_dir: Directory to cache data. Defaults to None.
            receptors_to_include: Receptors to include. Defaults to ("TCR", "BCR").
            test: Whether to run in test mode. Defaults to False.
            filter_10x: Whether to filter out 10X Genomics datasets. Defaults to False.
        """

        if not hasattr(self.__class__, "DB_NAME"):
            raise TypeError(f"Class {self.__class__.__name__} must define a 'DB_NAME' class attribute.")
        if not hasattr(self.__class__, "available_receptors"):
            raise TypeError(f"Class {self.__class__.__name__} must define an 'available_receptors' class attribute.")
        self._bc = bc
        self._test = test
        self._filter_10x = filter_10x
        self._cache_dir = cache_dir
        self._receptors = list(set(receptors_to_include) & set(self.available_receptors))
        self._table: pd.DataFrame | None = None
        self._airr_cells: list[AirrCell] | None = None
        self._metadata = {
            "db_name": self.DB_NAME,
            "version": "latest",
            "download_date": datetime.now().isoformat(),
            "source_url": None,
            "checksum": None,
            "has_changed": False,
        }
        self._table_path = self.get_latest_release(bc)

        # Handle both single paths and tuples (e.g., IEDB returns TCR/BCR pair)
        paths_to_check = [self._table_path] if isinstance(self._table_path, (str, Path)) else list(self._table_path)
        checksums = []
        mtimes = []
        if paths_to_check:
            for p in paths_to_check:
                if p and os.path.exists(str(p)):
                    checksums.append(get_file_checksum(str(p)))
                    mtimes.append(os.path.getmtime(str(p)))

        if checksums:
            self._metadata["checksum"] = checksums[0] if len(checksums) == 1 else checksums

        if mtimes:
            # Use the latest mtime if multiple files (cached files keep their original mtime)
            latest_mtime = max(mtimes)
            self._metadata["download_date"] = datetime.fromtimestamp(latest_mtime).isoformat()

    def set_metadata(self, version: str = None, source_url: str = None, previous_version: str = None):
        """
        Sets the metadata for the adapter.

        Args:
            version: The version of the database. Defaults to None.
            source_url: The URL of the source. Defaults to None.
            previous_version: The version of the database in the previous release. Defaults to None.
        """
        if version:
            self._metadata["version"] = version
        if source_url:
            self._metadata["source_url"] = source_url
        if previous_version is not None and version is not None:
            self._metadata["has_changed"] = version != previous_version

    @property
    def metadata(self) -> dict[str, Any]:
        """
        Property to get the adapter metadata.

        Returns:
            The metadata dictionary.
        """
        return self._metadata

    @property
    def db_name(self) -> str:
        """
        Property to get the database name.

        Returns:
            The database name.
        """
        return self.DB_NAME

    @property
    def receptors(self) -> list[str]:
        """
        Property to get the available receptor types.

        Returns:
            List of receptor types.
        """
        return self.available_receptors

    @property
    def table(self) -> pd.DataFrame:
        """
        Property to get the data table. Reads the table if not already read.

        Returns:
            The data table.
        """
        if self._table is None:
            self._table = self.read_table(self._bc, self._table_path, self._receptors, self._test)
            if self._filter_10x:
                # Filter out 10X Genomics dataset as it has been criticized for poor confidence.
                self._table = self._table[
                    ~(self._table["PMID"] == "no_pmid_1036521")
                    & ~self._table["PMID"].astype(str).str.contains("https://www.10xgenomics.com", na=False)
                ]

            # Filter out entries without receptor information (both chains missing) or without epitope information
            self._table = self._table[
                ~(self._table[REGISTRY_KEYS.CHAIN_1_CDR3_KEY].isna() & self._table[REGISTRY_KEYS.CHAIN_2_CDR3_KEY].isna())
                & ~self._table[REGISTRY_KEYS.EPITOPE_KEY].isna()
            ]

        return self._table

    @property
    def cache_dir(self) -> str:
        """
        Property to get the cache directory.

        Returns:
            The cache directory.
        """
        if self._cache_dir is None:
            self._cache_dir = platformdirs.user_cache_dir("iggytop")
            os.makedirs(self._cache_dir, exist_ok=True)
        return Path(self._cache_dir)

    @property
    def airr_cells(self) -> list[AirrCell] | None:
        """
        Property to get the list of AIRR cells.

        Returns:
            The list of AIRR cells.
        """
        if self._airr_cells is None:
            self._airr_cells = []

            # Using itertuples() for better performance on large DataFrames
            for row in tqdm(self.table.itertuples(), total=self.table.shape[0], desc=f"Processing {self.DB_NAME} entries"):
                c1_cdr3 = getattr(row, REGISTRY_KEYS.CHAIN_1_CDR3_KEY, None)
                c2_cdr3 = getattr(row, REGISTRY_KEYS.CHAIN_2_CDR3_KEY, None)
                if (pd.isnull(c2_cdr3) and pd.isnull(c1_cdr3)) or pd.isnull(getattr(row, "epitope_sequence", None)):
                    continue  # skip entries without chains or without epitope

                idx = row.Index
                cell = AirrCell(cell_id=str(idx))

                if not pd.isnull(c1_cdr3):
                    alpha_chain = AirrCell.empty_chain_dict()
                    alpha_chain.update(
                        {
                            "locus": getattr(row, REGISTRY_KEYS.CHAIN_1_TYPE_KEY, None),
                            "junction_aa": c1_cdr3,
                            "v_call": getattr(row, REGISTRY_KEYS.CHAIN_1_V_GENE_KEY, None),
                            "j_call": getattr(row, REGISTRY_KEYS.CHAIN_1_J_GENE_KEY, None),
                            "consensus_count": 0,
                            "productive": True,
                        }
                    )
                    cell.add_chain(alpha_chain)

                if not pd.isnull(c2_cdr3):
                    beta_chain = AirrCell.empty_chain_dict()
                    beta_chain.update(
                        {
                            "locus": getattr(row, REGISTRY_KEYS.CHAIN_2_TYPE_KEY, None),
                            "junction_aa": c2_cdr3,
                            "v_call": getattr(row, REGISTRY_KEYS.CHAIN_2_V_GENE_KEY, None),
                            "j_call": getattr(row, REGISTRY_KEYS.CHAIN_2_J_GENE_KEY, None),
                            "consensus_count": 0,
                            "productive": True,
                        }
                    )
                    cell.add_chain(beta_chain)

                # Source organism data is on the chain level for now
                c1_org = getattr(row, REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY, None)
                c2_org = getattr(row, REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY, None)

                resolved_org = None
                if pd.isnull(c1_org):
                    resolved_org = c2_org
                elif pd.isnull(c2_org):
                    resolved_org = c1_org
                elif c1_org == c2_org:
                    resolved_org = c1_org
                else:
                    if c2_org in c1_org or c1_org in c2_org:
                        resolved_org = c1_org if len(str(c1_org)) > len(str(c2_org)) else c2_org
                    else:
                        resolved_org = f"{c1_org} X {c2_org}"  # Transgenic or ambiguous case, keep both organisms in the string for now
                if not pd.isnull(resolved_org):
                    cell["source_organism"] = resolved_org

                for f in REGISTRY_KEYS:
                    # f is a column name (value from REGISTRY_KEYS)
                    if "chain" in f:
                        continue

                    val = getattr(row, f, None)
                    if val is not None and not pd.isnull(val):
                        cell[f] = val
                self._airr_cells.append(cell)

        return self._airr_cells

    @abstractmethod
    def get_latest_release(self, bc: BioCypher) -> str | tuple[str, ...]:
        """
        Abstract method to get the latest release of the data.

        Args:
            bc: An instance of the BioCypher class.

        Returns:
            Path to the latest release file(s).
        """
        pass

    @abstractmethod
    def read_table(self, bc: BioCypher, table_path: str | tuple[str, ...], receptors: list[str], test: bool = False) -> pd.DataFrame:
        """
        Abstract method to read and harmonize the data table from the source.

        Args:
            bc: An instance of the BioCypher class.
            table_path: Path to the data table file(s).
            receptors: List of receptor types to include in the table.
            test: Whether to run in test mode. Defaults to False.

        Returns:
            The data table.
        """
        pass

    @abstractmethod
    def get_nodes(self):
        """
        Abstract method to generate BioCypher nodes from the data.

        This method is intended to use _generate_nodes_from_table with the right parameters for each edge type.
        This requires parameters depending on the adapter used.

        Yields:
            tuple: A BioCypher node (id, type, properties).
        """
        pass

    @abstractmethod
    def get_edges(self):
        """
        Abstract method to generate BioCypher edges from the data.

        This method is intended to call _generate_edges_from_table with the right parameters for each edge type.
        This requires parameters depending on the adapter used.

        Yields:
            tuple: A BioCypher edge (id, source, target, type, properties).
        """
        pass

    def create_anndata(self) -> None:
        """
        Creates an Anndata object from the AIRR cell data and saves it to a file in the cache directory.
        """
        adata = from_airr_cells(self.airr_cells)
        index_chains(adata)

        # Convert object columns to string to avoid serialization issues with h5py (e.g. for PMID)
        for col in adata.obs.columns:
            if adata.obs[col].dtype == object:
                adata.obs[col] = adata.obs[col].astype(str)

        adata.uns["DB"] = {"name": self.DB_NAME, "date_downloaded": datetime.now().isoformat()}

        anndata_path = Path(self.cache_dir) / f"{self.DB_NAME}_anndata.h5ad"
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
            subset_cols: List of columns to subset the table.
            unique_cols: List of columns to check for uniqueness. Defaults to None.
            property_cols: List of columns to include as properties. Defaults to None.

        Yields:
            A tuple containing the node ID, node type, and properties.
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

        # Using itertuples() for better performance
        for row in subset_table.itertuples(index=False):
            if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in subset_cols:
                _type = getattr(row, REGISTRY_KEYS.CHAIN_1_TYPE_KEY)
            elif REGISTRY_KEYS.CHAIN_2_TYPE_KEY in subset_cols:
                _type = getattr(row, REGISTRY_KEYS.CHAIN_2_TYPE_KEY)
            else:
                _type = "epitope"

            # For TCR chains, use sequence + V gene + J gene as the identifier
            if _type.lower() != "epitope":
                # Get V gene and J gene if available
                v_gene_key = (
                    REGISTRY_KEYS.CHAIN_1_V_GENE_KEY if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in subset_cols else REGISTRY_KEYS.CHAIN_2_V_GENE_KEY
                )

                # Check if V gene is available in the row
                v_gene = getattr(row, v_gene_key, None)

                # Create an ID that includes V and J genes if available
                id_components = [_type.lower()]
                id_components.extend([str(getattr(row, col)) for col in unique_cols])
                if v_gene:
                    id_components.append(str(v_gene))

                _id = ":".join(id_components)
            else:
                # For epitopes and other types, keep the original ID format
                _id = ":".join([_type.lower(), *[str(getattr(row, col)) for col in unique_cols]])

            _props = {re.sub(r"chain_\d_", "", k): getattr(row, k) for k in property_cols}

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
            source_subset_cols: List of columns for the source node.
            target_subset_cols: List of columns for the target node.
            source_unique_cols: List of unique columns for the source node. Defaults to None.
            target_unique_cols: List of unique columns for the target node. Defaults to None.

        Yields:
            A tuple containing the edge ID, source ID, target ID, edge type, and properties.
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

        subset_table = self.table[source_subset_cols + target_subset_cols].dropna(subset=source_unique_cols + target_unique_cols)

        # Using itertuples() for better performance
        for row in subset_table.itertuples(index=False):
            node_data = {}
            for i in ["source", "target"]:
                cols = locals()[f"{i}_subset_cols"]
                unique_cols_list = locals()[f"{i}_unique_cols"]

                if REGISTRY_KEYS.CHAIN_1_TYPE_KEY in cols:
                    node_type = getattr(row, REGISTRY_KEYS.CHAIN_1_TYPE_KEY)
                    v_gene_key = REGISTRY_KEYS.CHAIN_1_V_GENE_KEY
                elif REGISTRY_KEYS.CHAIN_2_TYPE_KEY in cols:
                    node_type = getattr(row, REGISTRY_KEYS.CHAIN_2_TYPE_KEY)
                    v_gene_key = REGISTRY_KEYS.CHAIN_2_V_GENE_KEY
                else:
                    node_type = "epitope"
                    v_gene_key = None

                id_components = [node_type.lower()]
                id_components.extend([str(getattr(row, col)) for col in unique_cols_list])

                if v_gene_key:
                    v_gene = getattr(row, v_gene_key, None)
                    if v_gene:
                        id_components.append(str(v_gene))

                node_data[i] = {"id": ":".join(id_components), "type": node_type}

            _source_id = node_data["source"]["id"]
            _target_id = node_data["target"]["id"]
            _source_type = node_data["source"]["type"]
            _target_type = node_data["target"]["type"]

            _id = f"{_source_id}-{_target_id}"
            _type = f"{_source_type.lower()}_to_{_target_type.lower()}"

            yield (_id, _source_id, _target_id, _type, {})
