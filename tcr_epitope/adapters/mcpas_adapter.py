import os

import pandas as pd


class MCPASAdapter:

    DB_PATH = "data/mcpas_full.csv"

    def __init__(self, test: bool = False):
        if not os.path.exists(self.DB_PATH):
            raise FileNotFoundError(
                "MCPAS database not found. Please download from "
                "http://friedmanlab.weizmann.ac.il/McPAS-TCR/ and save as "
                "`mcpas_full.csv` in the `data/` directory."
            )
        table = pd.read_csv(self.DB_PATH, encoding="unicode_escape")
        if test:
            table = table.sample(frac=0.1)

        self.cdr3_alpha = table[["CDR3.alpha.aa"]].drop_duplicates().dropna()
        self.cdr3_beta = table[["CDR3.beta.aa"]].drop_duplicates().dropna()
        self.epitopes = table[["Epitope.peptide"]].drop_duplicates().dropna()

        self.alpha_beta_edges = table[["CDR3.alpha.aa", "CDR3.beta.aa"]].drop_duplicates().dropna()
        self.alpha_epitope_edges = table[["CDR3.alpha.aa", "Epitope.peptide"]].drop_duplicates().dropna()
        self.beta_epitope_edges = table[["CDR3.beta.aa", "Epitope.peptide"]].drop_duplicates().dropna()

    def get_nodes(self):
        for row in self.cdr3_alpha.itertuples():
            _id = "_".join(["TRA", row[1]])
            _type = "TRA"
            _props = {}

            yield (_id, _type, _props)

        for row in self.cdr3_beta.itertuples():
            _id = "_".join(["TRB", row[1]])
            _type = "TRB"
            _props = {}

            yield (_id, _type, _props)

        for row in self.epitopes.itertuples():
            _id = row[1]
            _type = "Epitope"
            _props = {}

            yield (_id, _type, _props)

    def get_edges(self):
        for row in self.alpha_beta_edges.itertuples():
            _from = "_".join(["TRA", row[1]])
            _to = "_".join(["TRB", row[2]])
            _type = "TRA_To_TRB"
            _props = {}

            yield (_from, _to, _type, _props)

        for row in self.alpha_epitope_edges.itertuples():
            _from = "_".join(["TRA", row[1]])
            _to = row[2]
            _type = "TCR_Sequence_To_Epitope"
            _props = {}

            yield (_from, _to, _type, _props)

        for row in self.beta_epitope_edges.itertuples():
            _from = "_".join(["TRB", row[1]])
            _to = row[2]
            _type = "TCR_Sequence_To_Epitope"
            _props = {}

            yield (_from, _to, _type, _props)
        