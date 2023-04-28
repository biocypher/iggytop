import os
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
import pooch
from github import Github

from typing import Optional


class VDJDBAdapter:

    REPO_NAME = "antigenomics/vdjdb-db"
    DB_DIR = "vdjdb_latest"
    DB_FNAME = "vdjdb.txt"

    def __init__(self, cache_dir: Optional[str] = None, test: bool = False):
        cache_dir = cache_dir or TemporaryDirectory().name
        db_path = self.download_latest_release(cache_dir)

        table = pd.read_csv(db_path, sep="\t")
        if test:
            table = table.sample(frac=0.1, random_state=123456)
        self.tcr_table = table
        
        cdr3 = table[["gene", "cdr3"]].drop_duplicates().dropna()
        self.cdr3_alpha = cdr3[cdr3["gene"] == "TRA"][["cdr3"]]
        self.cdr3_beta = cdr3[cdr3["gene"] == "TRB"][["cdr3"]]

        self.epitopes = table[["antigen.epitope"]].drop_duplicates().dropna()

        cdr3_epitope_edges = table[["gene", "cdr3", "antigen.epitope"]].drop_duplicates().dropna()
        self.alpha_epitope_edges = cdr3_epitope_edges[cdr3_epitope_edges["gene"] == "TRA"][["cdr3", "antigen.epitope"]]
        self.beta_epitope_edges = cdr3_epitope_edges[cdr3_epitope_edges["gene"] == "TRB"][["cdr3", "antigen.epitope"]]

    def download_latest_release(self, save_dir: str):
        repo = Github().get_repo(self.REPO_NAME)
        db_url = repo.get_latest_release().get_assets()[0].browser_download_url

        path = Path(pooch.retrieve(
            db_url,
            None,
            fname=self.DB_DIR,
            path=save_dir,
            processor=pooch.Unzip(),
        )[0]).parent

        return os.path.join(path, self.DB_FNAME)

    def get_nodes(self):
        for row in self.cdr3_alpha.itertuples():
            _id = "_".join(["TRA", row[1]])
            _type = 'TRA'
            _props = {
                'v_call' : self.tcr_table.loc[row.Index, 'v.segm'],
                'j_call' : self.tcr_table.loc[row.Index, 'j.segm'],
                'species' : self.tcr_table.loc[row.Index, 'species']
            }
            
            yield (_id, _type, _props)

        for row in self.cdr3_beta.itertuples():
            _id = "_".join(["TRB", row[1]])
            _type = 'TRB'
            _props = {
                'v_call' : self.tcr_table.loc[row.Index, 'v.segm'],
                'j_call' : self.tcr_table.loc[row.Index, 'j.segm'],
                'species' : self.tcr_table.loc[row.Index, 'species']
            }
            
            yield (_id, _type, _props)

        for row in self.epitopes.itertuples():
            _id = "_".join(["Epitope", row[1]])
            _type = 'Epitope'
            _props = {
                'protein' : self.tcr_table.loc[row.Index, 'antigen.gene'],
                'MHC_class' : self.tcr_table.loc[row.Index, 'mhc.class'],
                'MHC_gene_1' : self.tcr_table.loc[row.Index, 'mhc.a'],
                'MHC_gene_2' : self.tcr_table.loc[row.Index, 'mhc.b'],
                'species' : self.tcr_table.loc[row.Index, 'antigen.species'],
            }
            
            yield (_id, _type, _props)

    def get_edges(self):
        for row in self.alpha_epitope_edges.itertuples():
            _from = "_".join(["TRA", row[1]])
            _to = row[2]
            _type = 'TCR_Sequence_To_Epitope'
            _props = {}
            
            yield (None, _from, _to, _type, _props)

        for row in self.beta_epitope_edges.itertuples():
            _from = "_".join(["TRB", row[1]])
            _to = row[2]
            _type = 'TCR_Sequence_To_Epitope'
            _props = {}
            
            yield (None, _from, _to, _type, _props)
