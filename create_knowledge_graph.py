"""
This script creates a knowledge graph from various immunological databases
with receptor-epitope matching information and saves it in JSON format.
"""

from argparse import ArgumentParser

from biocypher import BioCypher

from tcr_epitope.adapters.cedar_adapter import CEDARAdapter
from tcr_epitope.adapters.iedb_adapter import IEDBAdapter
from tcr_epitope.adapters.mcpas_adapter import MCPASAdapter
from tcr_epitope.adapters.neotcr_adapter import NeoTCRAdapter
from tcr_epitope.adapters.tcr3d_adapter import TCR3DAdapter
from tcr_epitope.adapters.trait_adapter import TRAITAdapter
from tcr_epitope.adapters.utils import save_airr_cells_json
from tcr_epitope.adapters.vdjdb_adapter import VDJDBAdapter

parser = ArgumentParser()
parser.add_argument("--test", default=False)
parser.add_argument("--cache_dir", default=None)

args = parser.parse_args()

bc = BioCypher(cache_directory=args.cache_dir)

adapters = [
    VDJDBAdapter(bc, args.test),
    MCPASAdapter(bc, args.test),
    TRAITAdapter(bc, args.test),
    IEDBAdapter(bc, args.test),
    VDJDBAdapter(bc, args.test),
    MCPASAdapter(bc, args.test),
    TCR3DAdapter(bc, args.test),
    NeoTCRAdapter(bc, args.test),
    CEDARAdapter(bc, args.test),
]

for adapter in adapters:
    bc.add(adapter.get_nodes())
    bc.add(adapter.get_edges())

airr_cells = bc.get_kg()
bc.summary()

# This step required the final kg to be in the airr format (dbms specified in the biocypher config)
save_airr_cells_json(airr_cells, args.cache_dir)

# Usage example:
# python3 create_knowledge_graph.py --test True --cache_dir ./cache
