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
    VDJDBAdapter(bc, **kwargs),
    MCPASAdapter(bc, **kwargs),
    IEDBAdapter(bc, **kwargs),
]

for adapter in adapters:
    bc.add(adapter.get_nodes())
    bc.add(adapter.get_edges())


# for adapter in adapters:
#     bc.write_nodes(adapter.get_nodes())
#     bc.write_edges(adapter.get_edges())

# bc.write_import_call()


print(bc.to_df())

df = bc.to_df()

# nodes_df = bc.nodes
# edges_df = bc.edges

# print(nodes_df.head())
# print(edges_df.head())

# try:
#     df.to_csv('test_out.csv', index=False)
# except:
#     type(df)
#     # print(head(df))
#     df = pd.DataFrame.from_dict(bc.to_df(), index=[0])
#     df.to_csv('test_out.csv', index=[0])


bc.summary()