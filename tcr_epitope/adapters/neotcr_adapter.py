import os
import tempfile

import pandas as pd
import requests
from biocypher import BioCypher

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import harmonize_sequences


class NeoTCRAdapter(BaseAdapter):
    """BioCypher adapter for the NeoTCR dataset."""

    RAW_URL = "https://github.com/lyotvincent/NeoTCR/raw/main/data/NeoTCR%20data-20221220.xlsx"
    FILE_NAME = "NeoTCR_data-20221220.xlsx"

    def get_latest_release(self, bc: BioCypher) -> str:
        response = requests.get(self.RAW_URL)
        if response.status_code != 200:
            raise ConnectionError(f"Failed to download NeoTCR file: {self.RAW_URL}")

        tmp_dir = tempfile.mkdtemp()
        file_path = os.path.join(tmp_dir, self.FILE_NAME)

        with open(file_path, "wb") as f:
            f.write(response.content)

        return file_path

    def read_table(self, bc: BioCypher, table_path: str, test: bool = False) -> pd.DataFrame:
        table = pd.read_excel(table_path)

        if test:
            table = table.sample(frac=0.05, random_state=42)

        table = table.replace(["", "nan"], None).where(pd.notnull, None)

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
        table[REGISTRY_KEYS.PUBLICATION_KEY] = (
            table[REGISTRY_KEYS.PUBLICATION_KEY].astype(str).str.replace("PMID:", "").str.strip()
        )

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
            unique_cols=[REGISTRY_KEYS.CHAIN_1_CDR3_KEY],
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
            unique_cols=[REGISTRY_KEYS.CHAIN_2_CDR3_KEY],
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
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
                REGISTRY_KEYS.PUBLICATION_KEY,
            ],
            unique_cols=[REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY],
            property_cols=[
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
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
