"""
This script acts as a wrapper to create a knowledge graph from various immunological databases
with receptor-epitope matching information using Iggytop .
"""

import argparse

import platformdirs

from iggytop.io.create_knowledge_graph import DEFAULT_ADAPTERS, create_knowledge_graph


def main():
    parser = argparse.ArgumentParser(description="Create a knowledge graph from immunological databases.")
    parser.add_argument("--test-mode", action="store_true", default=False, help="Run in test mode with a small subset of data.")
    parser.add_argument(
        "--cache-dir",
        type=str,
        default=platformdirs.user_cache_dir("iggytop_airr"),
        help="Directory for caching results.",
    )
    parser.add_argument(
        "--adapters",
        nargs="+",
        default=DEFAULT_ADAPTERS,
        help="List of adapters to include (e.g., --adapters VDJDB CEDAR). Defaults to all.",
    )
    parser.add_argument(
        "--output-format",
        type=str,
        default="neo4j",
        choices=["airr", "neo4j", "networkx"],
        help="Format to write the knowledge graph in (default: neo4j).",
    )
    args = parser.parse_args()

    create_knowledge_graph(
        cache_dir=args.cache_dir,
        test_mode=args.test_mode,
        adapters_to_include=args.adapters,
        output_format=args.output_format,
    )


if __name__ == "__main__":
    main()
