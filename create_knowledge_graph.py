from argparse import ArgumentParser

from biocypher import BioCypher

from tcr_epitope.adapters.mcpas_adapter import MCPASAdapter
from tcr_epitope.adapters.utils import kg_pd_to_anndata
from tcr_epitope.adapters.vdjdb_adapter import VDJDBAdapter

parser = ArgumentParser()
parser.add_argument("--test", default=False)
parser.add_argument("--cache_dir", default=None)

args = parser.parse_args()

# kwargs = {
#     "test": args.test,
#     "cache_dir": args.cache_dir,
# }

bc = BioCypher(cache_directory=args.cache_dir)

adapters = [
    VDJDBAdapter(bc, args.test),
    MCPASAdapter(bc, args.test),
    # IEDBAdapter(bc, **kwargs),
]

for adapter in adapters:
    bc.add(adapter.get_nodes())
    bc.add(adapter.get_edges())

print(bc.to_df())

df = bc.to_df()

bc.summary()

# Convert DataFrame to AnnData
adata = kg_pd_to_anndata(df)
