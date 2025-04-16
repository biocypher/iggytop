import os
from pathlib import Path

import pandas as pd
from biocypher import BioCypher, FileDownload
from github import Github

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import get_iedb_ids_batch


class VDJDBAdapter(BaseAdapter):
    """BioCypher adapter for the VDJdb database (https://vdjdb.cdr3.net/).

    Parameters
    ----------
    bc
        BioCypher instance for DB download.
    test
        If `True`, only a subset of the data will be loaded for testing purposes.
    """

    REPO_NAME = "antigenomics/vdjdb-db"
    DB_DIR = "vdjdb_latest"
    DB_FNAME = "vdjdb.txt"

    def get_latest_release(self, bc: BioCypher) -> str:
        repo = Github().get_repo(self.REPO_NAME)
        db_url = repo.get_latest_release().get_assets()[0].browser_download_url
        # db_url = "https://github.com/antigenomics/vdjdb-db/releases/download/pyvdjdb-2025-02-21/vdjdb-2025-02-21.zip"

        vdjdb_resource = FileDownload(
            name=self.DB_DIR,
            url_s=db_url,
            lifetime=30,
            is_dir=False,
        )

        vdjdb_paths = bc.download(vdjdb_resource)

        db_dir = Path(vdjdb_paths[0]).parent
        for root, dirs, files in os.walk(db_dir):
            for file in files:
                if file == self.DB_FNAME:
                    db_path = os.path.join(root, file)

        if not db_path or not os.path.exists(db_path):
            raise FileNotFoundError(f"Failed to download VDJdb database from {db_url}")

        return db_path

    def read_table(self, bc: BioCypher, table_path: str, test: bool = False) -> pd.DataFrame:
        table = pd.read_csv(table_path, sep="\t")
        if test:
            table = table.sample(frac=0.1, random_state=42)
        table = table.where(pd.notnull(table), None)  # replace NaN with None

        rename_cols = {
            "gene": REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
            "cdr3": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "antigen.epitope": REGISTRY_KEYS.EPITOPE_KEY,
            "antigen.gene": REGISTRY_KEYS.ANTIGEN_KEY,
            "antigen.species": REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            "v.segm": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "j.segm": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            "species": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            "mhc.class": REGISTRY_KEYS.MHC_CLASS_KEY,
            "mhc.a": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "mhc.b": REGISTRY_KEYS.MHC_GENE_2_KEY,
        }
        table = table.rename(columns=rename_cols)
        table = table[list(rename_cols.values())]

        # validate peptide sequences
        sequence_cols = [
            REGISTRY_KEYS.EPITOPE_KEY,
            REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
        ]

        for col in sequence_cols:
            table[col] = table[col].apply(str)
            table[col] = table[col].apply(lambda x: x.upper())
            table[col] = table[col].apply(lambda x: "".join(x.split()))

        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY].apply(lambda x: x.lower())
        table[REGISTRY_KEYS.CHAIN_2_V_GENE_KEY] = table[REGISTRY_KEYS.CHAIN_1_V_GENE_KEY]
        table[REGISTRY_KEYS.CHAIN_2_J_GENE_KEY] = table[REGISTRY_KEYS.CHAIN_1_J_GENE_KEY]

        # Map epitope sequences to IEDB IDs
        valid_epitopes = table[REGISTRY_KEYS.EPITOPE_KEY].dropna().drop_duplicates().tolist()
        if len(valid_epitopes) > 0:
            epitope_map = get_iedb_ids_batch(bc, valid_epitopes)
        # Apply the mapping to create the IEDB ID column
        table[REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY] = table[REGISTRY_KEYS.EPITOPE_KEY].map(epitope_map)

        return table

    def get_nodes(self):
        # chain 1
        yield from self._generate_nodes_from_table(
            subset_cols=[
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            ],
            unique_cols=[
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            ],
            property_cols=[
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
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
            ],
            unique_cols=[
                REGISTRY_KEYS.EPITOPE_KEY,
            ],
            property_cols=[
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.MHC_GENE_2_KEY,
            ],
        )

    def get_edges(self):
        # chain 1 - epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            ],
            REGISTRY_KEYS.EPITOPE_KEY,
            source_unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
        )
