import logging
import os
from pathlib import Path
from typing import List, Literal, Optional

import pandas as pd
import pooch

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import validate_peptide_sequence

logger = logging.getLogger(__name__)


class IEDBAdapter(BaseAdapter):
    """
    BioCypher adapter for the Immune Epitope Database (IEDB)[https://www.iedb.org/].
    
    Parameters
    ----------
    cache_dir
        The directory to store the downloaded IEDB data in. If `None`, a temporary
        directory will be created.
    test
        If `True`, only a subset of the data will be loaded for testing purposes.
    """
    
    DB_URL = "https://www.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip"
    DB_DIR = "iedb_latest"
    TCR_FNAME = "tcr_full_v3.csv"
    BCR_FNAME = "bcr_full_v3.csv"

    def get_latest_release(self, save_dir: str) -> str:
        path = Path(pooch.retrieve(
            self.DB_URL,
            None,
            fname=self.DB_DIR,
            path=save_dir,
            processor=pooch.Unzip(),
        )[0]).parent

        return os.path.join(path, self.TCR_FNAME), os.path.join(path, self.BCR_FNAME)
    
    def read_table(self, table_path: str, test: bool = False) -> pd.DataFrame:
        tcr_table_path, bcr_table_path = table_path

        tcr_table = pd.read_csv(tcr_table_path, header=[0,1])
        tcr_table.columns = tcr_table.columns.map(' '.join)
        tcr_table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        tcr_table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY

        bcr_table = pd.read_csv(bcr_table_path, header=[0,1])
        bcr_table.columns = bcr_table.columns.map(' '.join)
        bcr_table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.IGH_KEY
        bcr_table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.IGL_KEY

        table = pd.concat([tcr_table, bcr_table], ignore_index=True)
        table = table.where(pd.notnull(table), None)  # replace NaN with None
        if test:
            table = table.sample(frac=0.1)

        rename_cols = {
            "Epitope Name": REGISTRY_KEYS.EPITOPE_KEY,
            "Epitope Source Molecule": REGISTRY_KEYS.ANTIGEN_KEY,
            "Epitope Source Organism": REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            "Chain 1 CDR1 Calculated": REGISTRY_KEYS.CHAIN_1_CDR1_KEY,
            "Chain 1 CDR2 Calculated": REGISTRY_KEYS.CHAIN_1_CDR2_KEY,
            "Chain 1 CDR3 Calculated": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "Chain 2 CDR1 Calculated": REGISTRY_KEYS.CHAIN_2_CDR1_KEY,
            "Chain 2 CDR2 Calculated": REGISTRY_KEYS.CHAIN_2_CDR2_KEY,
            "Chain 2 CDR3 Calculated": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "Chain 1 Calculated V Gene": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "Chain 1 Calculated J Gene": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            "Chain 2 Calculated V Gene": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            "Chain 2 Calculated J Gene": REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            "Chain 1 Organism IRI": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            "Chain 2 Organism IRI": REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY,
            REGISTRY_KEYS.CHAIN_1_TYPE_KEY: REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
            REGISTRY_KEYS.CHAIN_2_TYPE_KEY: REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
        }
        table = table.rename(columns=rename_cols)
        table = table[list(rename_cols.values())]

        # validate peptide sequences
        sequence_cols = [
            REGISTRY_KEYS.EPITOPE_KEY,
            REGISTRY_KEYS.CHAIN_1_CDR1_KEY,
            REGISTRY_KEYS.CHAIN_1_CDR2_KEY,
            REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            REGISTRY_KEYS.CHAIN_2_CDR1_KEY,
            REGISTRY_KEYS.CHAIN_2_CDR2_KEY,
            REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
        ]
        required_valid = [
            REGISTRY_KEYS.EPITOPE_KEY,
            REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
        ]

        # some sequences are in the format `sequence + additional info`
        def split_epitope_sequence(x: str) -> str:
            return x.split("+")[0]

        for col in sequence_cols:
            table[col] = table[col].apply(str)
            table[col] = table[col].apply(split_epitope_sequence)
            table[col] = table[col].apply(lambda x: x.upper())
            table[col] = table[col].apply(lambda x: "".join(x.split()))
            table[f"{col}_valid"] = table[col].apply(validate_peptide_sequence)

        for col in required_valid:
            table = table[table[f"{col}_valid"]]
            table = table[table[col] != "NAN"]

        return table
    
    def get_nodes(self):
        # chain 1
        yield from self._generate_nodes_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR1_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR2_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            ],
            unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
        )

        # chain 2
        yield from self._generate_nodes_from_table(
            subset_cols=[
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR1_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR2_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY,
            ],
            unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
        )

        # epitope
        yield from self._generate_nodes_from_table(
            subset_cols=[
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            ],
            unique_cols=REGISTRY_KEYS.EPITOPE_KEY,
        )

    def get_edges(self):        
        # chain 1 - chain 2
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            ],
            [
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            ],
            source_unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
        )

        # chain 1 - epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            ],
            REGISTRY_KEYS.EPITOPE_KEY,
            source_unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.EPITOPE_KEY,
        )

        # chain 2 - epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            ],
            REGISTRY_KEYS.EPITOPE_KEY,
            source_unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.EPITOPE_KEY,
        )
