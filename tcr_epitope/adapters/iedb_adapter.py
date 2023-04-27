import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional

import pandas as pd
import pooch

class IEDBAdapter:
    DB_URL = "https://www.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip"
    DB_DIR = "iedb_latest"
    TCR_FNAME = "tcr_full_v3.csv"
    BCR_FNAME = "bcr_full_v3.csv"

    def __init__(self, cache_dir: Optional[str] = None, test: bool = False):
        cache_dir = cache_dir or TemporaryDirectory().name
        tcr_path, bcr_path = self.download_latest_release(cache_dir)

        tcr_table = pd.read_csv(tcr_path, header=[0,1], index_col=0)
        tcr_table.columns = tcr_table.columns.map(' '.join)
        if test:
            tcr_table = tcr_table.sample(frac=0.1)
        self.tcr_table = tcr_table
        
        alpha = tcr_table[["Chain 1 CDR3 Calculated"]].drop_duplicates().dropna()
        beta = tcr_table[["Chain 2 CDR3 Calculated"]].drop_duplicates().dropna()
        epitopes = tcr_table[["Epitope Name"]].drop_duplicates().dropna()

        alpha_beta_edges = tcr_table[["Chain 1 CDR3 Calculated", "Chain 2 CDR3 Calculated"]].drop_duplicates().dropna()
        alpha_epitope_edges = tcr_table[["Chain 1 CDR3 Calculated", "Epitope Name"]].drop_duplicates().dropna()
        beta_epitope_edges = tcr_table[["Chain 2 CDR3 Calculated", "Epitope Name"]].drop_duplicates().dropna()

        self.tcr = {
            "alpha": alpha,
            "beta": beta,
            "epitopes": epitopes,
            "alpha_beta_edges": alpha_beta_edges,
            "alpha_epitope_edges": alpha_epitope_edges,
            "beta_epitope_edges": beta_epitope_edges,
        }

        bcr_table = pd.read_csv(bcr_path, header=[0,1], index_col=0)
        bcr_table.columns = bcr_table.columns.map(' '.join)
        if test:
            bcr_table = bcr_table.sample(frac=0.1)
        self.bcr_table = bcr_table

        heavy = bcr_table[["Chain 1 CDR3 Calculated"]].drop_duplicates().dropna()
        light = bcr_table[["Chain 2 CDR3 Calculated"]].drop_duplicates().dropna()
        epitopes = bcr_table[["Epitope Name"]].drop_duplicates().dropna()

        heavy_light_edges = bcr_table[["Chain 1 CDR3 Calculated", "Chain 2 CDR3 Calculated"]].drop_duplicates().dropna()
        heavy_epitope_edges = bcr_table[["Chain 1 CDR3 Calculated", "Epitope Name"]].drop_duplicates().dropna()
        light_epitope_edges = bcr_table[["Chain 2 CDR3 Calculated", "Epitope Name"]].drop_duplicates().dropna()

        self.bcr = {
            "heavy": heavy,
            "light": light,
            "epitopes": epitopes,
            "heavy_light_edges": heavy_light_edges,
            "heavy_epitope_edges": heavy_epitope_edges,
            "light_epitope_edges": light_epitope_edges,
        }

    def download_latest_release(self, save_dir: str) -> str:
        path = Path(pooch.retrieve(
            self.DB_URL,
            None,
            fname=self.DB_DIR,
            path=save_dir,
            processor=pooch.Unzip(),
        )[0]).parent

        return os.path.join(path, self.TCR_FNAME), os.path.join(path, self.BCR_FNAME)
        
    def get_nodes(self):
        for row in self.tcr["alpha"].itertuples():
            _id = "_".join(["TRA", row[1]])
            _type = "TRA"
            _props = {
                'v_call' : self.tcr_table.loc[row.index, 'Chain 1 Calculated V Gene'],
                'j_call' : self.tcr_table['Chain 1 Calculated J Gene'],
                'CDR1' : self.tcr_table['Chain 1 CDR1 Calculated'],
                'CDR2' : self.tcr_table['Chain 1 CDR2 Calculated'],
                'species' : self.tcr_table['Chain 1 Organism IRI']
            }

            yield (_id, _type, _props)

        for row in self.tcr["beta"].itertuples():
            _id = "_".join(["TRB", row[1]])
            _type = "TRB"
            _props = {
                'v_call' : self.tcr_table['Chain 2 Calculated V Gene'],
                'j_call' : self.tcr_table['Chain 2 Calculated J Gene'],
                'CDR1' : self.tcr_table['Chain 2 CDR1 Calculated'],
                'CDR2' : self.tcr_table['Chain 2 CDR2 Calculated'],
                'species' : self.tcr_table['Chain 2 Organism IRI']
            }

            yield (_id, _type, _props)

        for row in self.tcr["epitopes"].itertuples():
            _id = "_".join(["Epitope", row[1]])
            _type = "Epitope"
            _props = {
                'protein' : self.tcr_table['Epitope Source Molecule'],
                'species' : self.tcr_table['Epitope Source Organism'],
            }

            yield (_id, _type, _props)

        for row in self.bcr["heavy"].itertuples():
            _id = "_".join(["IGH", row[1]])
            _type = "IGH"
            _props = {}

            yield (_id, _type, _props)

        for row in self.bcr["light"].itertuples():
            _id = "_".join(["IGL", row[1]])
            _type = "IGL"
            _props = {}

            yield (_id, _type, _props)

        for row in self.bcr["epitopes"].itertuples():
            _id = "_".join(["Epitope", row[1]])
            _type = "Epitope"
            _props = {}

            yield (_id, _type, _props)

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
