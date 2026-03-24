import os

import pandas as pd
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import get_mhc_class, get_tissue_source, harmonize_sequences


class TRAITAdapter(BaseAdapter):
    """A Comprehensive Database for T-cell Receptor-Antigen Interactions `TRAIT <https://pgx.zju.edu.cn/traitdb/>`_."""

    DB_URL = "https://pgx.zju.edu.cn/download.trait/Interactive_TCR-pMHC_Pairs.zip_20250312.zip"
    """URL to download the TRAIT database."""

    DB_DIR = "trait_latest"
    """Directory name for the downloaded database."""

    DB_NAME = "TRAIT"
    """Name of the database."""

    available_receptors = ["TCR"]
    """Receptor types available in TRAIT."""

    def get_latest_release(self, bc: BioCypher) -> str:
        """
        Retrieves the latest release of the TRAIT database.

        Args:
            bc: An instance of the BioCypher class.

        Returns:
            Path to the latest release file.
        """
        self.set_metadata(source_url=self.DB_URL)
        trait_resource = FileDownload(
            name=self.DB_DIR,
            url_s=self.DB_URL,
            lifetime=30,
            # is_dir=False,
        )

        trait_paths = bc.download(trait_resource)

        if not trait_paths:
            raise FileNotFoundError(f"Failed to download TRAIT database from {self.DB_URL}")

        # Flatten to actual file(s), not just the top-level dir
        final_files = []
        for path in trait_paths:
            if os.path.isdir(path):
                # walk recursively to get files inside
                for root, _, files in os.walk(path):
                    for file in files:
                        final_files.append(os.path.join(root, file))
            else:
                if not path.endswith(".zip"):
                    final_files.append(path)

        return final_files[0]

    def read_table(self, bc: BioCypher, table_path: str, receptors: list[str], test: bool = False) -> pd.DataFrame:
        """
        Reads and processes the TRAIT table from the downloaded database file.

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
            table = table.sample(frac=0.01, random_state=42)
        # Replace NaN and empty strings with None
        table = table.replace(["", "nan"], None).where(pd.notnull, None)

        rename_cols = {
            "CDR3α": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "CDR3β": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "Epitope": REGISTRY_KEYS.EPITOPE_KEY,
            "Epitope_gene": REGISTRY_KEYS.ANTIGEN_KEY,
            "Epitope_species": REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            "MHC_class": REGISTRY_KEYS.MHC_CLASS_KEY,
            "MHC_A": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "MHC_B": REGISTRY_KEYS.MHC_GENE_2_KEY,
            "TRAV": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "TRAJ": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            "TRBV": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            "TRBJ": REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            "Species": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            "PubMed.ID": REGISTRY_KEYS.PUBLICATION_KEY,
            "Tissue": REGISTRY_KEYS.TISSUE_KEY,
        }

        table = table.rename(columns=rename_cols)
        table = table[list(rename_cols.values())]

        # Remove 'PMID:' prefix from reference IDs
        table[REGISTRY_KEYS.PUBLICATION_KEY] = (
            table[REGISTRY_KEYS.PUBLICATION_KEY].astype(str).str.replace(r"^PMID:", "", regex=True).replace("None", None)
        )

        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY
        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY]

        table_preprocessed = harmonize_sequences(bc, table)
        table_preprocessed[REGISTRY_KEYS.MHC_CLASS_KEY] = table_preprocessed[REGISTRY_KEYS.MHC_GENE_1_KEY].apply(get_mhc_class)

        table_preprocessed[REGISTRY_KEYS.TISSUE_KEY] = table_preprocessed[REGISTRY_KEYS.TISSUE_KEY].apply(get_tissue_source)

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
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.MHC_GENE_2_KEY,
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
