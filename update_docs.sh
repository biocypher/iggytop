#!/bin/bash

# Regenerate API documentation, excluding specific files
uv run sphinx-apidoc -o docs/source/api . --exclude create_knowledge_graph.py --exclude create_knowledge_graphx.py

# Build the documentation
uv run sphinx-build -b html docs/source docs/build