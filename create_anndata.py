import platformdirs
import importlib.resources
from biocypher import BioCypher

from iggytop.adapters.cedar_adapter import CEDARAdapter
from iggytop.adapters.iedb_adapter import IEDBAdapter
from iggytop.adapters.mcpas_adapter import MCPASAdapter
from iggytop.adapters.neotcr_adapter import NeoTCRAdapter
from iggytop.adapters.tcr3d_adapter import TCR3DAdapter
from iggytop.adapters.trait_adapter import TRAITAdapter
from iggytop.adapters.vdjdb_adapter import VDJDBAdapter
from iggytop.adapters.utils import _set_up_config

import os
import time

start_time = time.time()


cache_dir = platformdirs.user_cache_dir("iggytop_anndata_fresh")
test_mode = False
adapters_to_include = ["VDJDB", "IEDB", "MCPAS",  "TCR3D", "NeoTCR", "CEDAR","TRAIT"]
output_format = "airr"  # doesnt matter here

config_path = _set_up_config(output_format, cache_dir)

schema_config_path = os.path.join(cache_dir, 'schema_config.yaml')
with importlib.resources.open_text('iggytop.config', 'schema_config.yaml') as schema_file:
    with open(schema_config_path, 'w') as cache_schema_file:
        cache_schema_file.write(schema_file.read())

bc = BioCypher(biocypher_config_path=config_path,
                schema_config_path=schema_config_path,
                cache_directory=cache_dir)


adapter_classes = {
    "VDJDB": VDJDBAdapter,
    "MCPAS": MCPASAdapter,
    "TRAIT": TRAITAdapter,
    "IEDB": IEDBAdapter,
    "TCR3D": TCR3DAdapter,
    "NeoTCR": NeoTCRAdapter,
    "CEDAR": CEDARAdapter,
}


selected_adapters = [
    adapter_classes[name]
    for name in adapters_to_include
    if name in adapter_classes
]

print(f"Initialization took {time.time() - start_time:.2f} seconds")

for AdapterClass in selected_adapters:
    start_time = time.time()
    adapter = AdapterClass(bc, cache_dir, test_mode)
    adapter.create_anndata()
    print(f"Execution took {time.time() - start_time:.2f} seconds")

