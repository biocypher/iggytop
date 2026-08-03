#!/bin/bash -eu

# Imports a pre-built knowledge graph (CSVs produced by create_release.py, output-format
# "neo4j") that the Dockerfile downloaded from a GitHub release and baked into the image at
# /kg. No adapters run here: the graph is taken as-is from the release asset.
#
# /kg is deliberately outside of $NEO4J_HOME (/var/lib/neo4j): the base image's entrypoint
# unconditionally chowns and chmods 700 everything under $NEO4J_HOME on startup (when running
# as root), which would otherwise silently rewrite ownership/permissions on the CSV directory.

KG_DIR="/kg"
SOURCE_SCRIPT="$KG_DIR/neo4j-admin-import-call.sh"
IMPORT_SCRIPT="import-call.sh"  # writable copy: $KG_DIR is baked into a read-only image layer
MARKER="/data/.kg_imported"

if [ -f "$MARKER" ]; then
  echo "Knowledge graph already imported (found $MARKER), skipping neo4j-admin import."
  echo "To force a re-import, run 'docker compose down -v'."
  exit 0
fi

if [ ! -f "$SOURCE_SCRIPT" ]; then
  echo "No import script found at $SOURCE_SCRIPT - did the image build's knowledge graph download fail?"
  exit 1
fi

echo "Importing knowledge graph from $KG_DIR ..."
# The script ships with the absolute cache-dir path baked in from wherever it was generated;
# repoint every "<...>/knowledge_graph/" prefix at our mounted location instead.
sed -E "s#[^\",]+/knowledge_graph/#${KG_DIR}/#g" "$SOURCE_SCRIPT" > "$IMPORT_SCRIPT"
bash "$IMPORT_SCRIPT"
touch "$MARKER"
