import platformdirs
import importlib.resources
from biocypher import BioCypher

import logging
logger = logging.getLogger(__name__) #inherits the log config from biocypher


from iggytop.adapters.vdjdb_adapter import VDJDBAdapter
import yaml
import os
import time

def _set_up_config(output_format, cache_dir):
    # Load the base configuration
    with importlib.resources.open_text('iggytop.config', 'biocypher_config.yaml') as file:
        config = yaml.safe_load(file)
    if output_format:
        # Check if the output format is allowed
        allowed_formats = ['airr', 'networkx', 'neo4j', 'docker']
        if output_format not in allowed_formats:
            raise ValueError(f"Invalid output_format: {output_format}. Allowed formats are: {allowed_formats}")
        # Modify the configuration 
        if output_format == 'neo4j' or output_format == 'networkx':
            config['biocypher']['online_mode'] = False

        config['biocypher']['dbms'] = output_format
        if output_format == 'docker':
                with importlib.resources.open_text('iggytop.config', 'biocypher_docker_config.yaml') as d_file:
                    config = yaml.safe_load(d_file)


    # Ensure the cache directory exists
    os.makedirs(cache_dir, exist_ok=True)

    # Save the modified configuration to the cache directory
    modified_config_path = os.path.join(cache_dir, 'biocypher_config.yaml')
    with open(modified_config_path, 'w') as file:
        yaml.safe_dump(config, file)

    return modified_config_path

start_time = time.time()

cache_dir = platformdirs.user_cache_dir("iggytop_anndata")
test_mode = False
adapters_to_include = ["VDJDB"]
output_format = "airr"  # Options: "airr", "networkx", "neo4j"

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
    #"MCPAS": MCPASAdapter,
    #"TRAIT": TRAITAdapter,
    #"IEDB": IEDBAdapter,
    #"TCR3D": TCR3DAdapter,
    #"NeoTCR": NeoTCRAdapter,
    #"CEDAR": CEDARAdapter,
}


selected_adapters = [
    adapter_classes[name]
    for name in adapters_to_include
    if name in adapter_classes
]

for AdapterClass in selected_adapters:
    adapter = AdapterClass(bc,cache_dir, True)

adapter.create_anndata()

logger.info(f"Execution took {time.time() - start_time:.2f} seconds")