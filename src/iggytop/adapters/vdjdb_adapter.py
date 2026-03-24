import json
import os
from pathlib import Path

import pandas as pd
import requests
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import get_mhc_class, get_tissue_source, harmonize_sequences


class VDJDBAdapter(BaseAdapter):
    """
    BioCypher adapter for the `VDJDB database <https://vdjdb.cdr3.net/>`_.

    This adapter handles the downloading, reading, and processing of the VDJDB database.
    """

    REPO_NAME = "antigenomics/vdjdb-db"
    """GitHub repository name for the VDJDB database."""

    DB_DIR = "vdjdb_latest"
    """Directory name for the downloaded database."""

    DB_FNAME = "vdjdb.txt"
    """File name of the database."""

    DB_NAME = "VDJDB"
    """Name of the database."""

    available_receptors = ["TCR"]
    """Receptor types available in VDJDB."""

    def get_latest_release(self, bc: BioCypher) -> str:
        """
        Retrieves the latest release of the VDJDB database from GitHub.

        Args:
            bc: An instance of the BioCypher class.

        Returns:
            The file path of the downloaded database.

        Raises:
            FileNotFoundError: If the database file cannot be found after downloading.
        """
        # Use GitHub REST API directly to avoid PyGithub authentication issues
        api_url = f"https://api.github.com/repos/{self.REPO_NAME}/releases/latest"
        headers = {"Accept": "application/vnd.github.v3+json", "User-Agent": "iggytop-adapter"}
        github_token = os.getenv("GITHUB_TOKEN")

        # Try with token first if available
        if github_token:
            headers["Authorization"] = f"token {github_token}"

        response = requests.get(api_url, headers=headers, timeout=30)

        # If we get 401 (unauthorized), the token might be invalid/expired
        # For public repos, we can retry without authentication
        if response.status_code == 401:
            if github_token:
                # Token is invalid, try without it for public repo access
                clean_headers = {"Accept": "application/vnd.github.v3+json", "User-Agent": "iggytop-adapter"}
                response = requests.get(api_url, headers=clean_headers, timeout=30)

        response.raise_for_status()

        release_data = response.json()
        version_tag = release_data.get("tag_name", "latest")
        if not release_data.get("assets"):
            raise FileNotFoundError(f"No assets found in latest release for {self.REPO_NAME}")

        db_url = release_data["assets"][0]["browser_download_url"]

        self.set_metadata(version=version_tag, source_url=db_url)

        vdjdb_resource = FileDownload(
            name=self.DB_DIR,
            url_s=db_url,
            lifetime=30,
            is_dir=False,
        )

        vdjdb_paths = bc.download(vdjdb_resource)
        db_path = None
        db_dir = Path(vdjdb_paths[0]).parent
        for root, _dirs, files in os.walk(db_dir):
            for file in files:
                if file == self.DB_FNAME:
                    db_path = os.path.join(root, file)

        if not db_path or not os.path.exists(db_path):
            raise FileNotFoundError(f"Failed to download VDJDB database from {db_url}")

        return db_path

    def read_table(self, bc: BioCypher, table_path: str, receptors: list[str], test: bool = False) -> pd.DataFrame:
        """
        Reads and processes the VDJdb table from the downloaded database file.

        Args:
            bc: An instance of the BioCypher class.
            table_path: Path to the table file.
            receptors: List of receptor types to include in the table (Ignored as only TCR is available).
            test: If `True`, loads only a subset of the data for testing (default is False).

        Returns:
            A DataFrame containing the processed table data.

        Raises:
            FileNotFoundError: If the table file cannot be found.
        """
        table = pd.read_csv(table_path, sep="\t")
        if test:
            table = table.sample(frac=0.01, random_state=42)

        table = table.replace(["", "nan"], None).where(pd.notnull, None)

        table = self._transform_paired_data_efficient(table)

        rename_cols = {
            "cdr3_chain_1": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "v.segm_chain_1": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "j.segm_chain_1": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            "cdr3_chain_2": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "v.segm_chain_2": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            "j.segm_chain_2": REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            "species": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            "antigen.epitope": REGISTRY_KEYS.EPITOPE_KEY,
            "antigen.gene": REGISTRY_KEYS.ANTIGEN_KEY,
            "antigen.species": REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            "reference.id": REGISTRY_KEYS.PUBLICATION_KEY,
            "mhc.class": REGISTRY_KEYS.MHC_CLASS_KEY,
            "mhc.a": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "mhc.b": REGISTRY_KEYS.MHC_GENE_2_KEY,
            "meta.tissue": REGISTRY_KEYS.TISSUE_KEY,
        }
        table["meta.tissue"] = table["meta"].apply(
            lambda x: json.loads(x).get("tissue") if isinstance(x, str) else (x.get("tissue") if isinstance(x, dict) else None)
        )
        table = table.rename(columns=rename_cols)
        table = table[list(rename_cols.values())]

        # Remove 'PMID:' prefix from reference IDs
        table[REGISTRY_KEYS.PUBLICATION_KEY] = (
            table[REGISTRY_KEYS.PUBLICATION_KEY].astype(str).str.replace(r"^PMID:", "", regex=True).replace("None", None)
        )

        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY]
        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY

        table[REGISTRY_KEYS.MHC_CLASS_KEY] = table[REGISTRY_KEYS.MHC_CLASS_KEY].apply(get_mhc_class)

        table[REGISTRY_KEYS.TISSUE_KEY] = table[REGISTRY_KEYS.TISSUE_KEY].apply(get_tissue_source)

        # Preprocesses CDR3 sequences, epitope sequences, and gene names
        table_preprocessed = harmonize_sequences(bc, table)

        return table_preprocessed

    def _transform_paired_data_efficient(self, df):
        """
        Efficient transformation that handles ALL cases correctly.
        This is required to properly pair TRA and TRB chains based on complex.id. (They are on different rows in the raw database)
        Args:
            df: The input DataFrame containing the VDJdb data.

        Returns:
            A DataFrame with the transformed paired data.
        """

        # 1. Separate unpaired (complex.id == 0)
        unpaired = df[df["complex.id"] == 0].copy()

        # 2. Find ACTUALLY paired data (duplicated non-zero complex.ids)
        paired_mask = (df["complex.id"] != 0) & (df["complex.id"].duplicated(keep=False))
        paired = df[paired_mask].copy()

        # 3. Find incomplete pairs (non-zero, non-duplicated complex.ids)
        incomplete = df[(df["complex.id"] != 0) & (~paired_mask)].copy()

        result_parts = []

        # Process complete pairs
        if len(paired) > 0:
            # Check which complex.ids have both TRA and TRB
            chain_counts = paired.groupby("complex.id")["gene"].apply(set)
            complete_mask = chain_counts.apply(lambda x: {"TRA", "TRB"}.issubset(x))
            complete_complexes = chain_counts[complete_mask].index

            if len(complete_complexes) > 0:
                complete_data = paired[paired["complex.id"].isin(complete_complexes)]
                tra_complete = complete_data[complete_data["gene"] == "TRA"]
                trb_complete = complete_data[complete_data["gene"] == "TRB"]

                merge_cols = [
                    "complex.id",
                    "antigen.epitope",
                    "antigen.gene",
                    "antigen.species",
                    "mhc.class",
                    "mhc.a",
                    "mhc.b",
                    "reference.id",
                ]

                paired_result = tra_complete.merge(
                    trb_complete[merge_cols + ["cdr3", "v.segm", "j.segm"]],
                    on=merge_cols,
                    suffixes=("_chain_1", "_chain_2"),
                )
                result_parts.append(paired_result)

        # Process all single chains (unpaired + incomplete pairs)
        all_singles = pd.concat([unpaired, incomplete], ignore_index=True) if len(incomplete) > 0 else unpaired

        if len(all_singles) > 0:
            single_tra = self._process_single_chain(all_singles[all_singles["gene"] == "TRA"], "tra")
            single_trb = self._process_single_chain(all_singles[all_singles["gene"] == "TRB"], "trb")
            result_parts.extend([single_tra, single_trb])

        # Combine all
        result_parts = [df for df in result_parts if len(df) > 0]
        return pd.concat(result_parts, ignore_index=True) if result_parts else df

    def _process_single_chain(self, df, chain_type):
        """
        Process single chain data (TRA or TRB only).

        Args:
            df: The input DataFrame containing the single chain data.
            chain_type: The type of chain, either "tra" or "trb".

        Returns:
            A DataFrame with the processed single chain data.
        """
        if len(df) == 0:
            return df

        result = df.copy()
        if chain_type == "tra":
            result["cdr3_chain_1"] = result["cdr3"]
            result["v.segm_chain_1"] = result["v.segm"]
            result["j.segm_chain_1"] = result["j.segm"]
            result["cdr3_chain_2"] = None
            result["v.segm_chain_2"] = None
            result["j.segm_chain_2"] = None
        else:  # trb
            result["cdr3_chain_2"] = result["cdr3"]
            result["v.segm_chain_2"] = result["v.segm"]
            result["j.segm_chain_2"] = result["j.segm"]
            result["cdr3_chain_1"] = None
            result["v.segm_chain_1"] = None
            result["j.segm_chain_1"] = None

        return result.drop(columns=["cdr3", "v.segm", "j.segm"])

    def get_nodes(self):
        """
        Generates nodes for the VDJdb data.

        This method yields node data for chain 1, chain 2, and epitopes.
        """

        # chain 1
        yield from self._generate_nodes_from_table(
            subset_cols=[
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            ],
            unique_cols=[
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            ],
            property_cols=[
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            ],
        )

        # chain 2
        yield from self._generate_nodes_from_table(
            subset_cols=[
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY,
            ],
            unique_cols=[
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            ],
            property_cols=[
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY,
            ],
        )

        # epitope
        yield from self._generate_nodes_from_table(
            subset_cols=[
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.MHC_GENE_2_KEY,
                REGISTRY_KEYS.TISSUE_KEY,
                REGISTRY_KEYS.PUBLICATION_KEY,
            ],
            unique_cols=[
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
            ],
            property_cols=[
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.MHC_GENE_2_KEY,
                REGISTRY_KEYS.TISSUE_KEY,
                REGISTRY_KEYS.PUBLICATION_KEY,
            ],
        )

    def get_edges(self):
        """
        Generates edges for the VDJdb data.

        This method yields edge data for chain 1 to chain 2, chain 1 to epitope, and chain 2 to epitope.
        """

        # chain 1 to chain 2
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            ],
            [
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            ],
            source_unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
        )

        # chain 1 to epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            ],
            [REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY],
            source_unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
        )

        # chain 2 to epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            ],
            [REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY],
            source_unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
        )
