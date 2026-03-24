import pandas as pd
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import get_mhc_class, harmonize_sequences


class TCR3DAdapter(BaseAdapter):
    """BioCypher adapter for the `TCR3D <https://tcr3d.ibbr.umd.edu>`_ database."""

    DB_URL = "https://tcr3d.ibbr.umd.edu/static/download/tcr_complexes_data.tsv"
    """URL to download the TCR3D database."""
    DB_DIR = "tcr3d_latest"
    """Directory name for the downloaded database."""
    DB_NAME = "TCR3D"
    """Name of the database."""
    available_receptors = ["TCR"]
    """Receptor types available in TCR3D."""

    def get_latest_release(self, bc: BioCypher) -> str:
        """
        Retrieves the latest release of the TCR3D database.

        Args:
            bc: An instance of the BioCypher class.

        Returns:
            Path to the latest release file.
        """
        self.set_metadata(source_url=self.DB_URL)
        tcr3d_resource = FileDownload(
            name=self.DB_DIR,
            url_s=self.DB_URL,
            lifetime=30,
            is_dir=False,
        )

        tcr3d_path = bc.download(tcr3d_resource)

        if not tcr3d_path:
            raise FileNotFoundError(f"Failed to download TCR3D database from {self.DB_URL}")

        return tcr3d_path[0]

    def read_table(self, bc: BioCypher, table_path: str, receptors: list[str], test: bool = False) -> pd.DataFrame:
        """
        Reads and processes the TCR3D table from the downloaded database file.

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

        # Replace missing values
        table = table.replace(["", "nan", "n.a.", "null"], None).where(pd.notnull, None)

        rename_cols = {
            "TCR_complex": REGISTRY_KEYS.MHC_CLASS_KEY,
            "CDR3_alpha": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "TRAV_gene": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "CDR3_beta": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "TRBV_gene": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            "Epitope": REGISTRY_KEYS.EPITOPE_KEY,
            "MHC_allele": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "TCR_organism": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            "Pubmed": REGISTRY_KEYS.PUBLICATION_KEY,
        }

        table = table.rename(columns=rename_cols)
        table = table[list(rename_cols.values())]

        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY
        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY]

        table[REGISTRY_KEYS.CHAIN_1_J_GENE_KEY] = None
        table[REGISTRY_KEYS.CHAIN_2_J_GENE_KEY] = None

        # For the rows with multiple epitopes, separate them into multiple rows
        table[REGISTRY_KEYS.EPITOPE_KEY] = table[REGISTRY_KEYS.EPITOPE_KEY].apply(
            lambda x: x.split(",") if x is not None and "," in x else x
        )
        table = table.explode(REGISTRY_KEYS.EPITOPE_KEY).reset_index(drop=True)

        table[REGISTRY_KEYS.MHC_CLASS_KEY] = table[REGISTRY_KEYS.MHC_CLASS_KEY].apply(get_mhc_class)

        # Create a column placeholder for the antigen species
        table_preprocessed = harmonize_sequences(bc, table)

        return table_preprocessed

    def get_nodes(self):
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
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
                REGISTRY_KEYS.PUBLICATION_KEY,
            ],
            unique_cols=[
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
            ],
            property_cols=[
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
                REGISTRY_KEYS.PUBLICATION_KEY,
            ],
        )

    def get_edges(self):
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
