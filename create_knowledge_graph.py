"""
This script creates a knowledge graph from various immunological databases
with receptor-epitope matching information and saves it in JSON format.
"""

import platformdirs

from iggytop.io.create_knowledge_graph import create_knowledge_graph

create_knowledge_graph(cache_dir=platformdirs.user_cache_dir("iggytop_dev"), adapters_to_include = ["VDJDB", "MCPAS", "IEDB", "TCR3D", "NEOTCR", "CEDAR", "TRAIT"])
