import pandas as pd
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import harmonize_sequences, normalize_table_strings


class NEOTCRAdapter(BaseAdapter):
    """BioCypher adapter for the `NeoTCR <http://neotcrdb.bioxai.cn/home>`_ dataset."""

    RAW_URL = "https://github.com/lyotvincent/NeoTCR/raw/main/data/NeoTCR%20data-20221220.xlsx"
    """URL to download the NeoTCR database."""
    FILE_NAME = "NeoTCR_data-20221220.xlsx"
    """File name of the NeoTCR database."""
    DB_VERSION = "2022-12-20"
    """Version of the NeoTCR database, taken from the date embedded in `RAW_URL`."""
    DB_NAME = "NEOTCR"
    """Name of the database."""
    DB_DIR = "neotcr_latest"
    """Directory name for the downloaded database."""
    available_receptors = ["TCR"]
    """Receptor types available in NeoTCR."""

    def get_latest_release(self, bc: BioCypher) -> str:
        """
        Retrieves the latest release of the NeoTCR database.

        Args:
            bc: An instance of the BioCypher class.

        Returns:
            Path to the latest release file.
        """
        self.set_metadata(version=self.DB_VERSION, source_url=self.RAW_URL)
        # Download NEOTCR
        neotcr_resource = FileDownload(
            name=self.DB_DIR,
            url_s=self.RAW_URL,
            lifetime=30,
            is_dir=False,
        )

        neotcr_path = bc.download(neotcr_resource)

        if not neotcr_path:
            raise FileNotFoundError(f"Failed to download NEOTCR database from {self.RAW_URL}")

        return neotcr_path[0]

    def read_table(self, bc: BioCypher, table_path: str, receptors: list[str], test: bool = False) -> pd.DataFrame:
        """
        Reads and processes the NeoTCR table from the downloaded database file.

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
        table = pd.read_excel(table_path)

        if test:
            table = table.sample(frac=0.05, random_state=42)

        table = normalize_table_strings(table)

        # Rename and harmonize columns
        rename_cols = {
            "TRA_CDR3": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "TRAV": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "TRAJ": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            "TRB_CDR3": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "TRBV": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            "TRBJ": REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            "Neoepitope": REGISTRY_KEYS.EPITOPE_KEY,
            "Antigen": REGISTRY_KEYS.ANTIGEN_KEY,
            "HLA Allele": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "Source": REGISTRY_KEYS.TISSUE_KEY,
            "PubMed ID": REGISTRY_KEYS.PUBLICATION_KEY,
        }

        table = table.rename(columns=rename_cols)
        table = table.replace("n.a.", None)

        # Add organism (human) and TCR types
        table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY] = "Homo sapiens"
        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = "Homo sapiens"
        table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY] = "Homo sapiens"

        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY

        # For the rows with multiple epitopes, separate them into multiple rows
        table[REGISTRY_KEYS.EPITOPE_KEY] = table[REGISTRY_KEYS.EPITOPE_KEY].apply(
            lambda x: x.split(",") if x is not None and "," in x else x
        )
        table = table.explode(REGISTRY_KEYS.EPITOPE_KEY).reset_index(drop=True)

        # Trim Pubmed IDs
        table[REGISTRY_KEYS.PUBLICATION_KEY] = table[REGISTRY_KEYS.PUBLICATION_KEY].astype(str).str.replace("PMID:", "").str.strip()

        table_preprocessed = harmonize_sequences(bc, table)

        return table_preprocessed

    def get_nodes(self):
        """Yield BioCypher nodes generated via OntoWeaver."""
        nodes, _ = self._get_ontoweaver_kg()
        yield from nodes

    def get_edges(self):
        """Yield BioCypher edges generated via OntoWeaver."""
        _, edges = self._get_ontoweaver_kg()
        yield from edges
