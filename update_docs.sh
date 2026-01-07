#!/bin/bash

# Build the documentation
rm -rf docs/build
uv run sphinx-build -b html docs/source docs/build
