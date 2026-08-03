import logging
import os
from pathlib import Path

import pandas as pd
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import get_pmids_batch, harmonize_sequences, normalize_table_strings

logger = logging.getLogger(__name__)


class CEDARAdapter(BaseAdapter):
    """BioCypher adapter for the Cancer Epitope Database and Analysis Resource `CEDAR <https://cedar.iedb.org/>`_."""

    DB_URL = "https://cedar.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip"
    """URL to download the CEDAR database."""
    DB_DIR = "cedar_latest"
    """Directory name for the downloaded database."""
    DB_NAME = "CEDAR"
    """Name of the database."""
    TCR_FNAME = "tcr_full_v3.csv"
    """File name of the TCR data in CEDAR."""
    BCR_FNAME = "bcr_full_v3.csv"
    """File name of the BCR data in CEDAR."""
    available_receptors = ["TCR", "BCR"]
    """Receptor types available in CEDAR."""

    def get_latest_release(self, bc: BioCypher) -> tuple[str, str]:
        """
        Retrieves the latest release of the CEDAR database.

        Args:
            bc: An instance of the BioCypher class.

        Returns:
            Paths to the latest TCR and BCR release files.
        """
        self.set_metadata(version="v3", source_url=self.DB_URL)
        # Download CEDAR
        cedar_resource = FileDownload(
            name=self.DB_DIR,
            url_s=self.DB_URL,
            lifetime=30,
            is_dir=False,
        )

        cedar_paths = bc.download(cedar_resource)
        db_dir = Path(cedar_paths[0]).parent
        for root, _dirs, files in os.walk(db_dir):
            for file in files:
                if file == self.TCR_FNAME:
                    tcr_path = os.path.join(root, file)
                elif file == self.BCR_FNAME:
                    bcr_path = os.path.join(root, file)

        if not tcr_path or not os.path.exists(tcr_path) or not bcr_path or not os.path.exists(bcr_path):
            raise FileNotFoundError(f"Failed to download CEDAR database from {self.DB_URL}")

        return tcr_path, bcr_path

    def read_table(
        self,
        bc: BioCypher,
        table_path: tuple[str, str],
        receptors: list[str],
        test: bool = False,
    ) -> pd.DataFrame:
        """
        Reads and processes the CEDAR table from the downloaded database file.

        Args:
            bc: An instance of the BioCypher class.
            table_path: Paths to the TCR and BCR table files.
            receptors: List of receptor types to include in the table.
            test: If `True`, loads only a subset of the data for testing (default is False).

        Returns:
            A DataFrame containing the processed table data.

        Raises:
            FileNotFoundError: If the table file cannot be found.
        """
        tcr_table_path, bcr_table_path = table_path

        tcr_table = None
        bcr_table = None

        if "TCR" in receptors:
            # We use dtype=str to avoid DtypeWarning and optimize memory usage for large files.
            tcr_table = pd.read_csv(tcr_table_path, header=[0, 1], dtype=str, escapechar="\\")
            tcr_table.columns = tcr_table.columns.map(" ".join)
            tcr_table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
            tcr_table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY

        if "BCR" in receptors:
            bcr_table = pd.read_csv(bcr_table_path, header=[0, 1], dtype=str, escapechar="\\")
            bcr_table.columns = bcr_table.columns.map(" ".join)
            bcr_table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.IGH_KEY
            bcr_table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.IGL_KEY

        # Combine tables that were loaded
        tables_to_concat = [t for t in [tcr_table, bcr_table] if t is not None]
        if not tables_to_concat:
            return pd.DataFrame()
        table = pd.concat(tables_to_concat, ignore_index=True)

        if test:
            table = table.sample(frac=0.01, random_state=42)
        table = normalize_table_strings(table)

        # Same logic as in IEDB
        rename_cols = {
            "Epitope Name": REGISTRY_KEYS.EPITOPE_KEY,
            "Epitope CEDAR IRI": REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
            "Epitope Source Molecule": REGISTRY_KEYS.ANTIGEN_KEY,
            "Epitope Source Organism": REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            "Assay MHC Allele Names": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "Chain 1 CDR3 Calculated": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "Chain 2 CDR3 Calculated": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "Chain 1 Junction Calculated": REGISTRY_KEYS.CHAIN_1_JUNCTION_AA_KEY,
            "Chain 2 Junction Calculated": REGISTRY_KEYS.CHAIN_2_JUNCTION_AA_KEY,
            "Chain 1 Calculated V Gene": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "Chain 1 Calculated J Gene": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            "Chain 2 Calculated V Gene": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            "Chain 2 Calculated J Gene": REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            "Chain 1 Organism IRI": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            "Chain 2 Organism IRI": REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY,
            REGISTRY_KEYS.CHAIN_1_TYPE_KEY: REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
            REGISTRY_KEYS.CHAIN_2_TYPE_KEY: REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
            "Reference CEDAR IRI": REGISTRY_KEYS.PUBLICATION_KEY,
        }

        table = table.rename(columns=rename_cols)
        table = table[list(rename_cols.values())]

        # CEDAR's own export occasionally leaves a stray literal `"` in this field (e.g.
        # `Chlorobium chlorochromatii"`), which would otherwise corrupt BioCypher's
        # neo4j-admin-import CSV output.
        table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY] = table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY].str.replace('"', "", regex=False)

        # Extract iedb ID from the url
        table[REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY] = "iedb:" + table[REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY].str.extract(r"/epitope/(\d+)$")[0]

        # Preprocesses CDR3 sequences, epitope sequences, and gene names
        table_preprocessed = harmonize_sequences(bc, table)

        ref_urls = table_preprocessed[REGISTRY_KEYS.PUBLICATION_KEY].dropna().unique().tolist()
        ref_map = get_pmids_batch(bc, ref_urls)
        table_preprocessed[REGISTRY_KEYS.PUBLICATION_KEY] = table_preprocessed[REGISTRY_KEYS.PUBLICATION_KEY].replace(ref_map)

        return table_preprocessed

    def get_nodes(self):
        """Yield BioCypher nodes generated via OntoWeaver."""
        nodes, _ = self._get_ontoweaver_kg()
        yield from nodes

    def get_edges(self):
        """Yield BioCypher edges generated via OntoWeaver."""
        _, edges = self._get_ontoweaver_kg()
        yield from edges
