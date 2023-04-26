from biocypher import BioCypher
from tcr_epitope.adapters.vdjdb_adapter import VDJDBAdapter

bc = BioCypher()

vdjdb_adapter = VDJDBAdapter(
    file_path='data/SearchTable-2023-04-26 14_23_16.112.tsv'
)

# Create a knowledge graph from the adapter
bc.write_nodes(vdjdb_adapter.get_nodes())
bc.write_edges(vdjdb_adapter.get_edges())

bc.summary()