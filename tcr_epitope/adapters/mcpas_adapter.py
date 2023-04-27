import pandas as pd


class MCPASAdapter:

    URL = "http://friedmanlab.weizmann.ac.il/McPAS-TCR/session/59ca1e2ae928627a163fef1ed83ae3c2/download/downloadDB?w="

    def __init__(self):
        table = pd.read_csv("data/mcpas_test.csv")

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
        