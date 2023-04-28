import os

import pandas as pd


class MCPASAdapter:
    """
    BioCypher adapter for the manually-curated catalogue of pathology-associated T cell 
    receptor sequences (McPAS-TCR)[http://friedmanlab.weizmann.ac.il/McPAS-TCR/].

    In order to use this adapter, please download the latest version of the database
    from http://friedmanlab.weizmann.ac.il/McPAS-TCR/ and save it as `mcpas_full.csv`
    in the `data/` directory.
    
    Parameters
    ----------
    cache_dir
        The directory to store the downloaded IEDB data in. If `None`, a temporary
        directory will be created.
    test
        If `True`, only a subset of the data will be loaded for testing purposes.
    """

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
        self.tcr_table = table

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
            _props = {
                'v_call' : self.tcr_table.loc[row.Index, 'TRAV'],
                'j_call' : self.tcr_table.loc[row.Index, 'TRAJ'],
                'species' : self.tcr_table.loc[row.Index, 'Species']
            }

            yield (_id, _type, _props)

        for row in self.cdr3_beta.itertuples():
            _id = "_".join(["TRB", row[1]])
            _type = "TRB"
            _props = {
                'v_call' : self.tcr_table.loc[row.Index, 'TRBV'],
                'j_call' : self.tcr_table.loc[row.Index, 'TRBJ'],
                'species' : self.tcr_table.loc[row.Index, 'Species']
            }

            yield (_id, _type, _props)

        for row in self.epitopes.itertuples():
            _id = "_".join(["Epitope", row[1]])
            _type = "Epitope"
            _props = {
                'protein' : self.tcr_table.loc[row.Index, 'Antigen.protein'],
                'MHC_gene_1' : self.tcr_table.loc[row.Index, 'MHC'],
                'species' : self.tcr_table.loc[row.Index, 'Pathology'],
            }

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
        