import os
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
import pooch
from github import Github

class VDJDBAdapter:

    REPO_NAME = "antigenomics/vdjdb-db"
    DB_DIR = "vdjdb_latest"
    DB_FNAME = "vdjdb_full.txt"

    def __init__(self):
        save_dir = TemporaryDirectory()
        db_path = self.download_latest_release(save_dir)
        table = pd.read_csv(db_path, sep="\t")

        cdr3 = table[["Gene", "CDR3"]].drop_duplicates().dropna()
        self.cdr3_alpha = cdr3[cdr3["Gene"] == "TRA"]["CDR3"]
        self.cdr3_beta = cdr3[cdr3["Gene"] == "TRB"]["CDR3"]

        self.epitopes = table[["Epitope"]].drop_duplicates().dropna()

        cdr3_epitope_edges = table[["Gene", "CDR3", "Epitope"]].drop_duplicates().dropna()
        self.alpha_epitope_edges = cdr3_epitope_edges[cdr3_epitope_edges["Gene"] == "TRA"][["CDR3", "Epitope"]]
        self.beta_epitope_edges = cdr3_epitope_edges[cdr3_epitope_edges["Gene"] == "TRB"][["CDR3", "Epitope"]]

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
        for row in self.cdr3_alpha.iteritems():
            _id = "_".join(["TRA", row[1]])
            _type = 'TRA'
            _props = {}
            
            yield (_id, _type, _props)

        for row in self.cdr3_beta.iteritems():
            _id = "_".join(["TRB", row[1]])
            _type = 'TRB'
            _props = {}
            
            yield (_id, _type, _props)

        for row in self.epitopes.itertuples():
            _id = row[1]
            _type = 'Epitope'
            _props = {}
            
            yield (_id, _type, _props)

    def get_edges(self):
        for row in self.alpha_epitope_edges.itertuples():
            _from = "_".join(["TRA", row.CDR3])
            _to = row.Epitope
            _type = 'TCR_Sequence_To_Epitope'
            _props = {}
            
            yield (None, _from, _to, _type, _props)

        for row in self.beta_epitope_edges.itertuples():
            _from = "_".join(["TRB", row.CDR3])
            _to = row.Epitope
            _type = 'TCR_Sequence_To_Epitope'
            _props = {}
            
            yield (None, _from, _to, _type, _props)
