""" this module contains the main pipeline for knowledge graph creation using iggytop adapters """

import importlib.resources
import math
from typing import List, Optional

import networkx as nx
import pandas as pd
import platformdirs
from biocypher import BioCypher

import logging
logger = logging.getLogger(__name__) #inherits the log config from biocypher

from iggytop.adapters.cedar_adapter import CEDARAdapter
from iggytop.adapters.iedb_adapter import IEDBAdapter
from iggytop.adapters.mcpas_adapter import MCPASAdapter
from iggytop.adapters.neotcr_adapter import NeoTCRAdapter
from iggytop.adapters.tcr3d_adapter import TCR3DAdapter
from iggytop.adapters.trait_adapter import TRAITAdapter
from iggytop.adapters.utils import save_airr_cells_json
from iggytop.adapters.vdjdb_adapter import VDJDBAdapter
import yaml
import os


def create_knowledge_graph(
    cache_dir: str = platformdirs.user_cache_dir("iggytop"),
    test_mode: bool = False,
    adapters_to_include: Optional[List[str]] = ["VDJDB", "MCPAS", "TRAIT", "IEDB", "TCR3D", "NeoTCR", "CEDAR"],
    output_format: str | None = None,
):
    """
    Generates the knowledge graph using specified adapters and saves it in the requested format.

    Args:
        cache_dir (str, optional): Directory to store cache and output files. Includes raw datasets and generated knowledge graphs (see logs for filenames). Defaults to user cache directory.
        test_mode (bool, optional): Test mode will use only 1% of the data for faster execution. Defaults to False.
        adapters_to_include (List[str], optional): List of adapter names to run.
                             Available adapters: ["VDJDB", "MCPAS", "TRAIT", "IEDB", "TCR3D", "NeoTCR", "CEDAR"].
                             Defaults to providing all available adapters.
        output_format (str, optional): Output format, currently either 'airr','neo4j' or 'networkx'
    """
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

    for AdapterClass in selected_adapters:
        adapter = AdapterClass(bc, test_mode)
        bc.add(adapter.get_nodes())
        bc._add_edges(adapter.get_edges())
        logger.info(f"Added data from {AdapterClass.__name__}")

    bc.summary()

    if output_format == "airr":
        airr_cells = bc.get_kg()
        # This step required the final kg to be in the airr format (dbms specified in the biocypher config)
        save_airr_cells_json(airr_cells, cache_dir)
    elif output_format == "networkx":
        #this is more of a hack, updates on BioCypher needed to make it work properly
        bc._in_memory_kg = bc._writer.in_memory_networkx_kg
        # Remove the 'nodes' attribute from the BioCypher instance if it exists, this is a workaround
        if hasattr(bc, '_nodes'):
            del bc._nodes
        iggytop_di_graph = bc.to_networkx()
        
        bc.summary()
        # Save the NetworkX DiGraph to a file in GraphML format
        output_file = f"{cache_dir}/iggytop_knowledge_graph_vgene.graphml"

        for n, data in iggytop_di_graph.nodes(data=True):
            # Get all keys where the value is None, else the file can not be saved
            keys_to_del = [k for k, v in data.items() if v is None]
            for k in keys_to_del:
                del data[k]

        for s, t, data in iggytop_di_graph.edges(data=True):
            # Get all keys where the value is None, else the file can not be saved
            keys_to_del = [k for k, v in data.items() if v is None or (isinstance(v, float) and math.isnan(v))]
            for k in keys_to_del:
                del data[k]
        
        nx.write_graphml(iggytop_di_graph, output_file)
        print(f"Knowledge graph saved to {output_file}")
    elif output_format == "neo4j" or output_format == "docker":
        bc.write_import_call()
        bc.summary()

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

<<<<<<< HEAD
        config['biocypher']['dbms'] = output_format
        if output_format == 'docker':
                with importlib.resources.open_text('iggytop.config', 'biocypher_docker_config.yaml') as d_file:
                    config = yaml.safe_load(d_file)
        if output_format == 'networkx' or output_format == 'neo4j':
            config['biocypher']['offline'] = True

    # Ensure the cache directory exists
    os.makedirs(cache_dir, exist_ok=True)

    # Save the modified configuration to the cache directory
    modified_config_path = os.path.join(cache_dir, 'biocypher_config.yaml')
    with open(modified_config_path, 'w') as file:
        yaml.safe_dump(config, file)

    return modified_config_path
