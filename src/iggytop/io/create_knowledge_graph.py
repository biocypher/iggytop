import logging
import os
from typing import List, Optional

import networkx as nx
import platformdirs
from biocypher import BioCypher

from iggytop.adapters.cedar_adapter import CEDARAdapter
from iggytop.adapters.iedb_adapter import IEDBAdapter
from iggytop.adapters.itrap_adapter import ITRAPAdapter
from iggytop.adapters.mcpas_adapter import MCPASAdapter
from iggytop.adapters.neotcr_adapter import NEOTCRAdapter
from iggytop.adapters.tcr3d_adapter import TCR3DAdapter
from iggytop.adapters.trait_adapter import TRAITAdapter
from iggytop.adapters.utils import _set_up_config, _set_up_schema, save_airr_cells_json
from iggytop.adapters.vdjdb_adapter import VDJDBAdapter

logger = logging.getLogger(__name__)  # inherits the log config from biocypher


def create_knowledge_graph(
    cache_dir: str = platformdirs.user_cache_dir("iggytop"),
    test_mode: bool = False,
    receptors_to_include: Optional[List[str]] = ["TCR", "BCR"],
    adapters_to_include: Optional[List[str]] = [
        "VDJDB",
        "MCPAS",
        "TRAIT",
        "IEDB",
        "TCR3D",
        "NEOTCR",
        "CEDAR",
        "ITRAP",
    ],
    output_format: str | None = None,
    filter_10x: bool = False,
):
    """
    Generates the knowledge graph using specified adapters and saves it in the requested format.

    Args:
        cache_dir (str, optional): Directory to store cache and output files.
            Includes raw datasets and generated knowledge graphs (see logs for filenames).
            Defaults to user cache directory.
        test_mode (bool, optional): Test mode will use only 1% of the data for faster execution.
            Defaults to False.
        receptors_to_include (List[str], optional): List of receptor types to include in the knowledge graph.
            Available receptor types: ["TCR", "BCR"].
            Defaults to including both TCR and BCR.
        adapters_to_include (List[str], optional): List of adapter names to run.
            Available adapters: ["VDJDB", "MCPAS", "TRAIT", "IEDB", "TCR3D", "NEOTCR", "CEDAR", "ITRAP"].
            Defaults to providing all available adapters.
        output_format (str, optional): Output format, currently either 'airr','neo4j' or 'networkx'
        filter_10x (bool, optional): Whether to filter out 10X Genomics datasets. Defaults to False.
    """
    os.makedirs(cache_dir, exist_ok=True)
    config_path = _set_up_config(output_format, cache_dir)

    schema_config_path = _set_up_schema(cache_dir)

    bc = BioCypher(biocypher_config_path=config_path, schema_config_path=schema_config_path, cache_directory=cache_dir)

    adapter_classes = {
        "VDJDB": VDJDBAdapter,
        "MCPAS": MCPASAdapter,
        "TRAIT": TRAITAdapter,
        "IEDB": IEDBAdapter,
        "TCR3D": TCR3DAdapter,
        "ITRAP": ITRAPAdapter,
        "NEOTCR": NEOTCRAdapter,
        "CEDAR": CEDARAdapter,
    }

    selected_adapters = [adapter_classes[name] for name in adapters_to_include if name in adapter_classes]
    selected_adapters = [a for a in selected_adapters if any(receptor in receptors_to_include for receptor in a.available_receptors)]

    for AdapterClass in selected_adapters:
        adapter = AdapterClass(bc, cache_dir, receptors_to_include, test_mode, filter_10x)
        bc.add(adapter.get_nodes())
        bc._add_edges(adapter.get_edges())  # or bc.add(adapter.get_edges()) if in online mode
        logger.info(f"Added data from {AdapterClass.__name__}")

    bc.summary()

    if output_format == "airr":
        airr_cells = bc.get_kg()
        # This step required the final kg to be in the airr format (dbms specified in the biocypher config)
        save_airr_cells_json(airr_cells, cache_dir)
    elif output_format == "networkx":
        bc._in_memory_kg = bc._writer.in_memory_networkx_kg
        # Remove the 'nodes' attribute from the BioCypher instance if it exists, this is a workaround
        if hasattr(bc, "_nodes"):
            del bc._nodes
        iggytop_di_graph = bc.to_networkx()

        bc.summary()
        # Save the NetworkX DiGraph to a file in GraphML format
        output_file = f"{cache_dir}/iggytop_knowledge_graph_vgene.graphml"

        for n, data in iggytop_di_graph.nodes(data=True):
            # Get all keys where the value is None, else the file can not be saved
            keys_to_del = [k for k, v in data.items() if v is None]
            print(f"Removing node attributes with None values for node {n}: {keys_to_del}")
            for k in keys_to_del:
                del data[k]

        nx.write_graphml(iggytop_di_graph, output_file)
        print(f"Knowledge graph saved to {output_file}")
    elif output_format == "neo4j" or output_format == "docker":
        bc.write_import_call()
        bc.summary()
