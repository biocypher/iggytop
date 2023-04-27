from biocypher import BioCypher
from tcr_epitope.adapters.iedb_adapter import IEDBAdapter
from tcr_epitope.adapters.vdjdb_adapter import VDJDBAdapter
from tcr_epitope.adapters.mcpas_adapter import MCPASAdapter

bc = BioCypher()

adapters = [
    VDJDBAdapter(test=True),
    MCPASAdapter(test=True),
    IEDBAdapter(test=True),
]

for adapter in adapters:
    bc.add(adapter.get_nodes())
    bc.add(adapter.get_edges())

print(bc.to_df())

entities = bc.to_df()
bc.show_ontology_structure()
