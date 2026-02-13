"""
This script creates a knowledge graph from various immunological databases
with receptor-epitope matching information and saves it in JSON format.
"""

import platformdirs

from iggytop.io.create_knowledge_graph import create_knowledge_graph

create_knowledge_graph(cache_dir=platformdirs.user_cache_dir("iggytop_dev"), adapters_to_include = ["VDJdb", "McPAS", "IEDB", "TCR3d", "NeoTCR", "CEDAR", "TRAIT"])
