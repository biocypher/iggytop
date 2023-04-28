import logging
import os
from pathlib import Path
from typing import List, Tuple

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
        tcr_table[REGISTRY_KEYS.TYPE_KEY] = "TCR"
        bcr_table = pd.read_csv(bcr_table_path, header=[0,1])
        bcr_table.columns = bcr_table.columns.map(' '.join)
        bcr_table[REGISTRY_KEYS.TYPE_KEY] = "BCR"

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
            "type": "type",
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
    
    def _generate_from_table(
        self, subset_cols: List[str], unique_cols: List[str]
    ) -> List[Tuple]:
        subset_table = self.table[subset_cols].drop_duplicates(
            subset=unique_cols
        )
        for _, row in subset_table.iterrows():
            
    
    def get_nodes(self):
        yield from self._generate_from_table(...)


        # chain 1 (alpha or heavy)
        chain_1_subset_cols = [
            REGISTRY_KEYS.TYPE_KEY,
            REGISTRY_KEYS.CHAIN_1_CDR1_KEY,
            REGISTRY_KEYS.CHAIN_1_CDR2_KEY,
            REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
        ]
        chain_1_df = self.table[chain_1_subset_cols].drop_duplicates(
            subset=[REGISTRY_KEYS.CHAIN_1_CDR3_KEY]
        )
        for _, row in chain_1_df.iterrows():
            _type = "TRA" if row["type"] == "TCR" else "IGH"
            _id = "_".join([_type, row["cdr3_alpha"]])
            _props = {
                "v_call": row["v_alpha"],
                "j_call": row["j_alpha"],
                "cdr1": row["cdr1_alpha"],
                "cdr2": row["cdr2_alpha"],
                "organism": row["organism_alpha"],
            }

            yield (_id, _type, _props)

        # chain 2 (beta or light)
        chain_2_subset_cols = [
            "type",
            "cdr1_beta",
            "cdr2_beta",
            "cdr3_beta",
            "v_beta",
            "j_beta",
            "organism_beta",
        ]
        chain_2_df = self.table[chain_2_subset_cols].drop_duplicates(
            subset=["cdr3_beta"]
        )
        for _, row in chain_2_df.iterrows():
            _type = "TRB" if row["type"] == "TCR" else "IGL"
            _id = "_".join([_type, row["cdr3_beta"]])
            _props = {
                "v_call": row["v_beta"],
                "j_call": row["j_beta"],
                "cdr1": row["cdr1_beta"],
                "cdr2": row["cdr2_beta"],
                "organism": row["organism_beta"],
            }

            yield (_id, _type, _props)

        # epitopes
        epitope_subset_cols = [
            "epitope_sequence",
            "antigen",
            "antigen_organism",
        ]
        epitope_df = self.table[epitope_subset_cols].drop_duplicates(
            subset=["epitope_sequence"]
        )
        for _, row in epitope_df.iterrows():
            _type = "epitope"
            _id = "_".join([_type, row["epitope_sequence"]])
            _props = {
                "antigen": row["antigen"],
                "organism": row["antigen_organism"],
            }

            yield (_id, _type, _props)

    def get_edges(self):        
        # alpha-beta
        alpha_beta_subset_cols = [
            "type",
            "cdr3_alpha",
            "cdr3_beta",
        ]
        alpha_beta_df = self.table[alpha_beta_subset_cols].drop_duplicates(
            subset=["cdr3_alpha", "cdr3_beta"]
        )
        for _, row in alpha_beta_df.iterrows():
            _from_type = "TRA" if row["type"] == "TCR" else "IGH"
            _from = "_".join([_from_type, row["cdr3_alpha"]])
            _to_type = "TRB" if row["type"] == "TCR" else "IGL"
            _to = "_".join([_to_type, row["cdr3_beta"]])
            _id = "-".join([_from, _to])
            _type = "_".join([_from_type, "to", "_to_type"])

            yield (_id, _from, _to, _type, {})

        # alpha-epitope
        alpha_epitope_subset_cols = [
            "type",
            "cdr3_alpha",
            "epitope_sequence",
        ]
        alpha_epitope_df = self.table[alpha_epitope_subset_cols].drop_duplicates(
            subset=["cdr3_alpha", "epitope_sequence"]
        )
        for _, row in alpha_epitope_df.iterrows():
            _from_type = "TRA" if row["type"] == "TCR" else "IGH"
            _from = "_".join([_from_type, row["cdr3_alpha"]])
            _to_type = "epitope"
            _to = "_".join([_type, row["epitope_sequence"]])
            _id = "-".join([_from, _to])
            _type = "_".join([_from_type, "to", "_to_type"])

            yield (_id, _from, _to, _type, {})
        
        # beta-epitope
        beta_epitope_subset_cols = [
            "type",
            "cdr3_beta",
            "epitope_sequence",
        ]
        beta_epitope_df = self.table[beta_epitope_subset_cols].drop_duplicates(
            subset=["cdr3_beta", "epitope_sequence"]
        )
        for _, row in beta_epitope_df.iterrows():
            _from_type = "TRB" if row["type"] == "TCR" else "IGL"
            _from = "_".join([_from_type, row["cdr3_beta"]])
            _to_type = "epitope"
            _to = "_".join([_type, row["epitope_sequence"]])
            _id = "-".join([_from, _to])
            _type = "_".join([_from_type, "to", "_to_type"])
            
            yield (_id, _from, _to, _type, {})
        



        

    def get_edges(self):
        for row in self.tcr["alpha_beta_edges"].itertuples():
            _from = "_".join(["TRA", row[1]])
            _to = "_".join(["TRB", row[2]])
            _type = "TRA_To_TRB"
            _props = {}

            yield (_from, _to, _type, _props)

        for row in self.tcr["alpha_epitope_edges"].itertuples():
            _from = "_".join(["TRA", row[1]])
            _to = row[2]
            _type = "TCR_Sequence_To_Epitope"
            _props = {}

            yield (_from, _to, _type, _props)
        
        for row in self.tcr["beta_epitope_edges"].itertuples():
            _from = "_".join(["TRB", row[1]])
            _to = row[2]
            _type = "TCR_Sequence_To_Epitope"
            _props = {}

            yield (_from, _to, _type, _props)
        
        for row in self.bcr["heavy_light_edges"].itertuples():
            _from = "_".join(["IGH", row[1]])
            _to = "_".join(["IGL", row[2]])
            _type = "IGH_To_IGL"
            _props = {}

            yield (_from, _to, _type, _props)
        
        for row in self.bcr["heavy_epitope_edges"].itertuples():
            _from = "_".join(["IGH", row[1]])
            _to = row[2]
            _type = "BCR_Sequence_To_Epitope"
            _props = {}

            yield (_from, _to, _type, _props)
        
        for row in self.bcr["light_epitope_edges"].itertuples():
            _from = "_".join(["IGL", row[1]])
            _to = row[2]
            _type = "BCR_Sequence_To_Epitope"
            _props = {}

            yield (_from, _to, _type, _props)
