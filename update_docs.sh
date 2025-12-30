#!/bin/bash

# Regenerate API documentation, pointing to the package source
uv run sphinx-apidoc -f -o docs/source/api src/iggytop

# Build the documentation
uv run sphinx-build -b html docs/source docs/build