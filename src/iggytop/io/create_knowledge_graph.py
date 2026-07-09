import logging
import os
from typing import List, Optional

import networkx as nx
import platformdirs
from biocypher import BioCypher

from iggytop.adapters.batcave_adapter import BATCAVEAdapter
from iggytop.adapters.cedar_adapter import CEDARAdapter
from iggytop.adapters.iedb_adapter import IEDBAdapter
from iggytop.adapters.itrap_adapter import ITRAPAdapter
from iggytop.adapters.mcpas_adapter import MCPASAdapter
from iggytop.adapters.neotcr_adapter import NEOTCRAdapter
from iggytop.adapters.tcr3d_adapter import TCR3DAdapter
from iggytop.adapters.trait_adapter import TRAITAdapter
from iggytop.adapters.utils import _TT_LOG_PATH, _flush_tt_warnings, _set_up_config, _set_up_schema, save_airr_cells_json
from iggytop.adapters.vdjdb_adapter import VDJDBAdapter

logger = logging.getLogger(__name__)  # inherits the log config from biocypher

# Single source of truth for the set of available adapters. Add new adapters here only;
# DEFAULT_ADAPTERS (and every script/function that lists "all adapters") derives from it.
ADAPTER_CLASSES = {
    "VDJDB": VDJDBAdapter,
    "MCPAS": MCPASAdapter,
    "TRAIT": TRAITAdapter,
    "IEDB": IEDBAdapter,
    "TCR3D": TCR3DAdapter,
    "ITRAP": ITRAPAdapter,
    "NEOTCR": NEOTCRAdapter,
    "CEDAR": CEDARAdapter,
    "BATCAVE": BATCAVEAdapter,
}

DEFAULT_ADAPTERS = list(ADAPTER_CLASSES.keys())


def build_adapters(
    cache_dir: str,
    test_mode: bool = False,
    receptors_to_include: Optional[List[str]] = None,
    adapters_to_include: Optional[List[str]] = None,
    bc: BioCypher | None = None,
):
    """
    Instantiates each selected adapter exactly once, triggering (lazily, on first access) the
    expensive per-source table build (`read_table`/`harmonize_sequences`, including IEDB lookups).

    The returned adapters can be driven into any number of downstream outputs (AnnData, AIRR JSON,
    knowledge graph via `write_knowledge_graph`) without repeating that work, since `table`,
    `airr_cells` and `_get_ontoweaver_kg` are all memoized per adapter instance.

    Args:
        cache_dir: Directory to store cache files.
        test_mode: Test mode will use only 1% of the data for faster execution. Defaults to False.
        receptors_to_include: List of receptor types to include. Defaults to including both TCR and BCR.
        adapters_to_include: List of adapter names to run. Defaults to providing all available adapters.
        bc: An existing BioCypher instance to reuse (e.g. built by the caller for a specific dbms).
            If not given, a default one is created purely for source-download caching purposes.

    Returns:
        A tuple of (BioCypher instance, list of instantiated adapters).
    """
    receptors_to_include = receptors_to_include or ["TCR", "BCR"]
    adapters_to_include = adapters_to_include or DEFAULT_ADAPTERS

    os.makedirs(cache_dir, exist_ok=True)

    if bc is None:
        config_path = _set_up_config(None, cache_dir)
        schema_config_path = _set_up_schema(cache_dir)
        bc = BioCypher(biocypher_config_path=config_path, schema_config_path=schema_config_path, cache_directory=cache_dir)
    print(f"Tidytcells standardization warnings: {_TT_LOG_PATH}")

    selected_adapters = [ADAPTER_CLASSES[name] for name in adapters_to_include if name in ADAPTER_CLASSES]
    selected_adapters = [a for a in selected_adapters if any(receptor in receptors_to_include for receptor in a.available_receptors)]

    adapters = [AdapterClass(bc, cache_dir, receptors_to_include, test_mode) for AdapterClass in selected_adapters]
    return bc, adapters


def write_knowledge_graph(adapters, cache_dir: str, output_format: str = "neo4j", output_directory: str | None = None) -> None:
    """
    Drives already-built adapters into a fresh, dbms-appropriate BioCypher instance and writes out
    the requested output format. Since `adapters` are pre-built, `adapter.get_nodes()`/`get_edges()`
    only replay already-memoized OntoWeaver results; no source table is re-read.

    For 'neo4j'/'networkx'/'airr', BioCypher's own CSV+import-call output is written to a predictable
    `<cache_dir>/knowledge_graph` directory (unless overridden), so it sits next to the release's
    other outputs and can be packaged as a release asset without importing into a running Neo4j
    instance. The 'docker' format keeps the fixed path baked into the Docker build instead.

    Args:
        adapters: Adapters previously created via `build_adapters`.
        cache_dir: Directory to store cache and output files.
        output_format: Output format, currently either 'airr', 'neo4j', 'networkx' or 'docker'.
        output_directory: Where to write the graph output. Defaults to `<cache_dir>/knowledge_graph`
            for every format except 'docker', which relies on its own fixed config-driven path.
    """
    config_path = _set_up_config(output_format, cache_dir)
    schema_config_path = _set_up_schema(cache_dir)
    if output_directory is None and output_format != "docker":
        output_directory = os.path.join(cache_dir, "knowledge_graph")
    bc = BioCypher(
        biocypher_config_path=config_path,
        schema_config_path=schema_config_path,
        cache_directory=cache_dir,
        output_directory=output_directory,
    )

    for adapter in adapters:
        bc.add(adapter.get_nodes())
        bc._add_edges(adapter.get_edges())  # or bc.add(adapter.get_edges()) if in online mode
        logger.info(f"Added data from {adapter.DB_NAME}")

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


def create_knowledge_graph(
    cache_dir: str = platformdirs.user_cache_dir("iggytop"),
    test_mode: bool = False,
    receptors_to_include: Optional[List[str]] = ["TCR", "BCR"],
    adapters_to_include: Optional[List[str]] = DEFAULT_ADAPTERS,
    output_format: str = "neo4j",
):
    """
    Generates the knowledge graph using specified adapters and saves it in the requested format.

    Thin convenience wrapper around `build_adapters` + `write_knowledge_graph` for standalone,
    graph-only runs. Callers that also need the AnnData/AIRR outputs from the same adapter tables
    (e.g. `create_release.py`) should call `build_adapters` once and pass the result to
    `write_knowledge_graph` directly, to avoid rebuilding the source tables.

    Args:
        cache_dir (str, optional): Directory to store cache and output files.
            Includes raw datasets and generated knowledge graphs (see logs for filenames).
            Defaults to user cache directory.
        test_mode (bool, optional): Test mode will use only 1% of the data for faster execution.
            Defaults to False.
        receptors_to_include (List[str], optional): List of receptor types to include in the knowledge graph.
            Available receptor types: ["TCR", "BCR"].
            Defaults to including both TCR and BCR.
        adapters_to_include (List[str], optional): List of adapter names to run. See `ADAPTER_CLASSES`
            for the available adapters. Defaults to `DEFAULT_ADAPTERS` (all of them).
        output_format (str, optional): Output format, currently either 'airr', 'neo4j', 'networkx' or 'docker'.
            Defaults to 'neo4j'.
    """
    _bc, adapters = build_adapters(cache_dir, test_mode, receptors_to_include, adapters_to_include)
    write_knowledge_graph(adapters, cache_dir, output_format)
    _flush_tt_warnings()
