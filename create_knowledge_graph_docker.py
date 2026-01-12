"""
This script creates a knowledge graph from various immunological databases
with receptor-epitope matching information and saves it in JSON format.
"""

import platformdirs

from iggytop.io.create_knowledge_graph import create_knowledge_graph

create_knowledge_graph(adapters_to_include= ["TCR3D"], output_format="docker")
