from biocypher import BioCypher
from tcr_epitope.adapters.vdjdb_adapter import VDJDBAdapter
from tcr_epitope.adapters.mcpas_adapter import MCPASAdapter

bc = BioCypher()


vdjdb_adapter = VDJDBAdapter()
mcpas_adapter = MCPASAdapter()


# Create a knowledge graph from the adapter
bc.write_nodes(vdjdb_adapter.get_nodes())
bc.write_edges(vdjdb_adapter.get_edges())

bc.write_nodes(mcpas_adapter.get_nodes())
bc.write_edges(mcpas_adapter.get_edges())

bc.write_import_call()
bc.summary()