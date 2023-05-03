from argparse import ArgumentParser

from biocypher import BioCypher
from tcr_epitope.adapters.iedb_adapter import IEDBAdapter
from tcr_epitope.adapters.vdjdb_adapter import VDJDBAdapter
from tcr_epitope.adapters.mcpas_adapter import MCPASAdapter

parser = ArgumentParser()
parser.add_argument("--test", default=False)
parser.add_argument("--cache_dir", default=None)

args = parser.parse_args()
kwargs = {
    "test": args.test,
    "cache_dir": args.cache_dir,
}

bc = BioCypher()

adapters = [
    VDJDBAdapter(**kwargs),
    MCPASAdapter(**kwargs),
    IEDBAdapter(**kwargs),
]

for adapter in adapters:
    bc.add(adapter.get_nodes())
    bc.add(adapter.get_edges())  # TODO no edges?

print(bc.to_df())
bc.summary()
