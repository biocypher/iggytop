from __future__ import annotations

import importlib.resources
import os
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Sequence, cast

import anndata as ad
import ontoweaver
import pandas as pd
import yaml
from scirpy.io._convert_anndata import from_airr_cells
from scirpy.io._datastructures import AirrCell
from scirpy.pp import index_chains
from tqdm.auto import tqdm

from . import ontoweaver_transformers  # noqa: F401 (registers the custom id transformers)
from .constants import REGISTRY_KEYS
from .utils import get_file_checksum

if TYPE_CHECKING:
    from biocypher import BioCypher
import platformdirs


def _load_ontoweaver_mapping(filename: str) -> dict:
    with importlib.resources.open_text("iggytop.config", filename) as f:
        return yaml.safe_load(f)


# One mapping per hub in the graph's hub-and-spoke hierarchy:
#   binding -> {receptor_complex, pmhc, database, PMID}
#   receptor_complex -> {chain_1, chain_2}
#   chain_N -> {v_gene, j_gene}
#   pmhc -> {epitope, mhc}
#   epitope -> {antigen}
# OntoWeaver maps one row-subject with radiating edges per pass, so each hub needs its own pass;
# the results of all passes are reconciled together in `_get_ontoweaver_kg`.
_ONTOWEAVER_MAPPINGS = [
    _load_ontoweaver_mapping(f"ontoweaver_mapping_{name}.yaml")
    for name in (
        "chain1",
        "chain2",
        "receptor_complex",
        "epitope",
        "antigen",
        "mhc",
        "pmhc",
        "binding",
        "database",
        "pmid",
    )
]


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
    ):
        """
        Initializes the BaseAdapter instance.

        Args:
            bc: An instance of the BioCypher class.
            cache_dir: Directory to cache data. Defaults to None.
            receptors_to_include: Receptors to include. Defaults to ("TCR", "BCR").
            test: Whether to run in test mode. Defaults to False.
        """

        if not hasattr(self.__class__, "DB_NAME"):
            raise TypeError(f"Class {self.__class__.__name__} must define a 'DB_NAME' class attribute.")
        if not hasattr(self.__class__, "available_receptors"):
            raise TypeError(f"Class {self.__class__.__name__} must define an 'available_receptors' class attribute.")
        self._bc = bc
        self._test = test
        self._cache_dir = cache_dir
        self._receptors = list(set(receptors_to_include) & set(self.available_receptors))
        self._table: pd.DataFrame | None = None
        self._airr_cells: list[AirrCell] | None = None
        self._ontoweaver_kg: tuple[list[tuple], list[tuple]] | None = None
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

    def set_metadata(
        self,
        version: str = None,
        source_url: str = None,
        previous_version: str = None,
        previous_checksum: str | list | None = None,
    ):
        """
        Sets the metadata for the adapter.

        Args:
            version: The version of the database. Defaults to None.
            source_url: The URL of the source. Defaults to None.
            previous_version: The version of the database in the previous release. Defaults to None.
            previous_checksum: The checksum(s) of the source file(s) in the previous release. Defaults to None.
        """
        if version:
            self._metadata["version"] = version
        if source_url:
            self._metadata["source_url"] = source_url

        if previous_checksum is not None:
            self._metadata["has_changed"] = self._metadata.get("checksum") != previous_checksum
        elif previous_version is not None:
            self._metadata["has_changed"] = self._metadata["version"] != previous_version

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

            # Filter out entries without receptor information (both chains missing) or without epitope information
            self._table = self._table[
                ~(self._table[REGISTRY_KEYS.CHAIN_1_CDR3_KEY].isna() & self._table[REGISTRY_KEYS.CHAIN_2_CDR3_KEY].isna())
                & ~self._table[REGISTRY_KEYS.EPITOPE_KEY].isna()
            ]

        return self._table

    def _get_ontoweaver_kg(self) -> tuple[list[tuple], list[tuple]]:
        """
        Generates BioCypher nodes and edges from the harmonized table via OntoWeaver.

        Runs the shared per-hub mappings in `_ONTOWEAVER_MAPPINGS` (driven purely by REGISTRY_KEYS
        column names, so they apply unmodified to any adapter's table) against a copy of the table
        augmented with a few precomputed columns (the `complete` flags, which need boolean
        combinations no single-column mapping can express, and the source database's identity,
        which lives in `self.metadata` rather than as a table column). The passes are reconciled
        together afterwards, so nodes/edges produced by more than one pass (e.g. a chain_2 node,
        created as a target by the receptor_complex mapping and as the subject of the chain_2
        mapping) are merged into one.

        Returns:
            A tuple of (nodes, edges), each a list of BioCypher tuples.
        """
        if self._ontoweaver_kg is None:
            table = self.table.copy()

            # Add any missing columns that OntoWeaver expects to see, but which may not be present in the table
            for _col in (REGISTRY_KEYS.MHC_GENE_2_KEY, REGISTRY_KEYS.TISSUE_KEY):
                if _col not in table.columns:
                    table[_col] = None

            chain_1_complete = (
                table[REGISTRY_KEYS.CHAIN_1_CDR3_KEY].notna()
                & table[REGISTRY_KEYS.CHAIN_1_V_GENE_KEY].notna()
                & table[REGISTRY_KEYS.CHAIN_1_J_GENE_KEY].notna()
            )
            chain_2_complete = (
                table[REGISTRY_KEYS.CHAIN_2_CDR3_KEY].notna()
                & table[REGISTRY_KEYS.CHAIN_2_V_GENE_KEY].notna()
                & table[REGISTRY_KEYS.CHAIN_2_J_GENE_KEY].notna()
            )
            receptor_complex_complete = chain_1_complete & chain_2_complete

            mhc_cols = [
                c for c in (REGISTRY_KEYS.MHC_CLASS_KEY, REGISTRY_KEYS.MHC_GENE_1_KEY, REGISTRY_KEYS.MHC_GENE_2_KEY) if c in table.columns
            ]
            mhc_present = table[mhc_cols].notna().any(axis=1) if mhc_cols else pd.Series(False, index=table.index)
            binding_complete = receptor_complex_complete & mhc_present

            # BioCypher's neo4j-admin-import writer declares these columns `:boolean` in the CSV
            # header but serializes values via plain str(), giving Python's capitalized "True"/
            # "False" — neo4j-admin only recognizes the lowercase literal "true" as true, so every
            # row was silently importing as false. Stringify to the exact literal it expects.
            table["chain_1_complete"] = chain_1_complete.map({True: "true", False: "false"})
            table["chain_2_complete"] = chain_2_complete.map({True: "true", False: "false"})
            table["receptor_complex_complete"] = receptor_complex_complete.map({True: "true", False: "false"})
            table["mhc_present"] = mhc_present.map({True: "true", False: "false"})
            table["binding_complete"] = binding_complete.map({True: "true", False: "false"})

            table["_db_name"] = self.DB_NAME
            table["_db_version"] = self.metadata.get("version") or "latest"

            nodes, edges = ontoweaver.extract(
                [(table, mapping) for mapping in _ONTOWEAVER_MAPPINGS],
                affix="none",
            )
            fnodes, fedges = ontoweaver.fusion.reconciliate(
                [n.as_tuple() for n in nodes],
                [e.as_tuple() for e in edges],
                reconciliate_sep=",",
                progress_bar=True,
            )
            self._ontoweaver_kg = (fnodes, fedges)
        return self._ontoweaver_kg

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
                            "junction_aa": getattr(row, REGISTRY_KEYS.CHAIN_1_JUNCTION_AA_KEY, None),
                            "cdr3_aa": c1_cdr3,
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
                            "junction_aa": getattr(row, REGISTRY_KEYS.CHAIN_2_JUNCTION_AA_KEY, None),
                            "cdr3_aa": c2_cdr3,
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
                    cell[REGISTRY_KEYS.SOURCE_ORGANISM_KEY] = resolved_org

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

        Adapters are expected to implement this via `_get_ontoweaver_kg`.

        Yields:
            tuple: A BioCypher node (id, type, properties).
        """
        pass

    @abstractmethod
    def get_edges(self):
        """
        Abstract method to generate BioCypher edges from the data.

        Adapters are expected to implement this via `_get_ontoweaver_kg`.

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
                adata.obs[col] = adata.obs[col].astype("string")

        adata.uns["DB"] = {"name": self.DB_NAME, "date_downloaded": datetime.now().isoformat()}

        anndata_path = Path(self.cache_dir) / f"{self.DB_NAME}_anndata.h5ad"
        # Opt in to writing pd.arrays.StringArray in obs/var.
        ad.settings.allow_write_nullable_strings = True
        adata.write_h5ad(cast(os.PathLike, anndata_path), compression="gzip")
        print(f"Saved Anndata to {anndata_path}")
