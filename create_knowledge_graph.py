"""
This script creates a knowledge graph from various immunological databases
with receptor-epitope matching information and saves it in JSON format.
"""

import argparse

import platformdirs
from iggytop.io.create_knowledge_graph import create_knowledge_graph


def main():
    parser = argparse.ArgumentParser(description="Create a knowledge graph from immunological databases.")
    parser.add_argument("--test-mode", action="store_true", help="Run in test mode with a small subset of data.")
    parser.add_argument(
        "--cache-dir",
        type=str,
        default=platformdirs.user_cache_dir("iggytop_dev"),
        help="Directory for caching results.",
    )
    args = parser.parse_args()

    create_knowledge_graph(
        cache_dir=args.cache_dir,
        test_mode=args.test_mode,
        adapters_to_include=["VDJDB", "MCPAS", "IEDB", "TCR3D", "NEOTCR", "CEDAR", "TRAIT", "ITRAP"],
    )


if __name__ == "__main__":
    main()
