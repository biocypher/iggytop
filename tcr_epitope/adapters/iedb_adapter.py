# https://www.iedb.org/
# Receptor	Receptor	Receptor	Receptor	Reference	Epitope	Epitope	Epitope	Epitope	Assay	Assay	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 1	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2	Chain 2
# Group IRI	IEDB Receptor ID	Reference Name	Type	IEDB IRI	IEDB IRI	Name	Source Molecule	Source Organism	Type	IEDB IDs	Type	Organism IRI	Nucleotide Sequence	Curated V Gene	Calculated V Gene	Curated D Gene	Calculated D Gene	Curated J Gene	Calculated J Gene	Protein Sequence	Protein IRI	CDR3 Curated	CDR3 Calculated	CDR3 Start Curated	CDR3 End Curated	CDR3 Start Calculated	CDR3 End Calculated	CDR1 Curated	CDR1 Calculated	CDR1 Start Curated	CDR1 End Curated	CDR1 Start Calculated	CDR1 End Calculated	CDR2 Curated	CDR2 Calculated	CDR2 Start Curated	CDR2 End Curated	CDR2 Start Calculated	CDR2 End Calculated	Type	Organism IRI	Nucleotide Sequence	Curated V Gene	Calculated V Gene	Curated D Gene	Calculated D Gene	Curated J Gene	Calculated J Gene	Protein Sequence	Protein IRI	CDR3 Curated	CDR3 Calculated	CDR3 Start Curated	CDR3 End Curated	CDR3 Start Calculated	CDR3 End Calculated	CDR1 Curated	CDR1 Calculated	CDR1 Start Curated	CDR1 End Curated	CDR1 Start Calculated	CDR1 End Calculated	CDR2 Curated	CDR2 Calculated	CDR2 Start Curated	CDR2 End Curated	CDR2 Start Calculated	CDR2 End Calculated
# http://www.iedb.org/receptor/47	57	KK50.4	alphabeta	http://www.iedb.org/reference/1004539	http://www.iedb.org/epitope/69921	VMAPRTLIL	HLA class I histocompatibility antigen, Cw-3 alpha chain precursor	Homo sapiens (human)	T cell	1548960, 1583178	alpha	http://purl.obolibrary.org/obo/NCBITaxon_9606		TRAV26-1*01	TRAV26-1*01			TRAJ37*01	TRAJ37*01	KTTQPPSMDCAEGRAANLPCNHSTISGNEYVYWYRQIHSQGPQYIIHGLKNNETNEMASLIITEDRKSSTLILPHATLRDTAVYYCIVVRSSNTGKLIFGQGTTLQVKPDIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPS	https://www.ncbi.nlm.nih.gov/protein/2ESV_D	IVVRSSNTGKLI	IVVRSSNTGKLI			87	98		TISGNEY			24	30		GLKNN			48	52	beta	http://purl.obolibrary.org/obo/NCBITaxon_9606		TRBV14*01	TRBV14*01			TRBJ2-3*01	TRBJ2-3*01	GVTQFPSHSVIEKGQTVTLRCDPISGHDNLYWYRRVMGKEIKFLLHFVKESKQDESGMPNNRFLAERTGGTYSTLKVQPAELEDSGVYFCASSQDRDTQYFGPGTRLTVLEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD	https://www.ncbi.nlm.nih.gov/protein/2ESV_D	ASSQDRDTQY	ASSQDRDTQY			91	100		SGHDN			25	29		FVKESK			47	52
# http://www.iedb.org/receptor/47	57	KK50.4	alphabeta	http://www.iedb.org/reference/1004539	http://www.iedb.org/epitope/69921	VMAPRTLIL	HLA class I histocompatibility antigen, Cw-3 alpha chain precursor	Homo sapiens (human)	T cell	1583178	alpha	http://purl.obolibrary.org/obo/NCBITaxon_9606			TRAV26-1*01				TRAJ37*01	KTTQPPSMDCAEGRAANLPCNHSTISGNEYVYWYRQIHSQGPQYIIHGLKNNETNEMASLIITEDRKSSTLILPHATLRDTAVYYCIVVRSSNTGKLIFGQGTTLQVKPDIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPS	https://www.ncbi.nlm.nih.gov/protein/2ESV_D		IVVRSSNTGKLI			87	98		TISGNEY			24	30		GLKNN			48	52	beta	http://purl.obolibrary.org/obo/NCBITaxon_9606			TRBV14*01				TRBJ2-3*01	GVTQFPSHSVIEKGQTVTLRCDPISGHDNLYWYRRVMGKEIKFLLHFVKESKQDESGMPNNRFLAERTGGTYSTLKVQPAELEDSGVYFCASSQDRDTQYFGPGTRLTVLEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD	https://www.ncbi.nlm.nih.gov/protein/2ESV_D		ASSQDRDTQY			91	100		SGHDN			25	29		FVKESK			47	52

import pandas as pd

class IEDBAdapter:
    def __init__(self, file_path: str = None) -> None:
        iedb = pd.read_csv(file_path, header=[0,1], index_col=0)
        iedb.columns = iedb.columns.map(' '.join)
        self.iedb = iedb
        self.cdr3_alpha = iedb[['Chain 1 Type','Chain 1 CDR3 Calculated']].dropna().drop_duplicates()
        self.cdr3_beta = iedb[['Chain 2 Type','Chain 2 CDR3 Calculated']].dropna().drop_duplicates()
        self.epitopes = iedb[['Epitope Name']]
        self.alpha_beta_edges = iedb[['Chain 1 CDR3 Calculated','Chain 2 CDR3 Calculated']].dropna().drop_duplicates()
        self.alpha_epitope_edges = iedb[['Chain 1 Type','Chain 1 CDR3 Calculated', 'Epitope Name']].dropna().drop_duplicates()
        self.beta_epitope_edges = iedb[['Chain 2 Type','Chain 2 CDR3 Calculated', 'Epitope Name']].dropna().drop_duplicates()

        
        
    def get_nodes(self):

        for row in self.cdr3_alpha.itertuples():
            _id = "_".join(["TRA", row[1]])
            _type = "TRA"
            _props = {
                'v_call' : self.iedb['Chain 1 Calculated V Gene'],
                'j_call' : self.iedb['Chain 1 Calculated J Gene'],
                'CDR1' : self.iedb['Chain 1 CDR1 Calculated'],
                'CDR2' : self.iedb['Chain 1 CDR2 Calculated'],
                'species' : self.iedb['Chain 1 Organism IRI']
            }

            yield (_id, _type, _props)

        for row in self.cdr3_beta.itertuples():
            _id = "_".join(["TRB", row[1]])
            _type = "TRB"
            _props = {
                'v_call' : self.iedb['Chain 2 Calculated V Gene'],
                'j_call' : self.iedb['Chain 2 Calculated J Gene'],
                'CDR1' : self.iedb['Chain 2 CDR1 Calculated'],
                'CDR2' : self.iedb['Chain 2 CDR2 Calculated'],
                'species' : self.iedb['Chain 2 Organism IRI']
            }

            yield (_id, _type, _props)

        for row in self.epitopes.itertuples():
            _id = row[1]
            _type = "Epitope"
            _props = {
                'protein' : self.iedb['Epitope Source Molecule'],
                'species' : self.iedb['Epitope Source Organism'],
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
