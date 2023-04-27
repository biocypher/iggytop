from biocypher import BioCypher
from tcr_epitope.adapters.iedb_adapter import IEDBAdapter
from tcr_epitope.adapters.vdjdb_adapter import VDJDBAdapter
from tcr_epitope.adapters.mcpas_adapter import MCPASAdapter

bc = BioCypher()

adapters = [
    VDJDBAdapter(),
    # MCPASAdapter(),
    IEDBAdapter(),
]

for adapter in adapters:
    bc.write_nodes(adapter.get_nodes())
    bc.write_edges(adapter.get_edges())

bc.write_import_call()
bc.summary()
