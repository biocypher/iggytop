# Where the baked-in knowledge graph comes from:
#   release (default) - a knowledge_graph.tar.gz asset from an IggyTop GitHub data release
#   local             - a `knowledge_graph` directory built on this machine, handed in as the
#                       `local-kg` build context (see docker-compose.yml)
ARG IGGYTOP_KG_SOURCE=release

FROM docker.io/neo4j:5.26-community AS base

COPY docker/* ./
RUN cat entrypoint_patch.sh | cat - /startup/docker-entrypoint.sh > docker-entrypoint.sh && \
    mv docker-entrypoint.sh /startup/ && \
    chmod +x /startup/docker-entrypoint.sh

# ---- KG source "release": download and unpack a data release's knowledge_graph.tar.gz.
# Done in a throwaway curl image rather than in the neo4j image below, so that image needs no
# extra packages installed and keeps neither curl nor the downloaded archive in its layers. ----
FROM docker.io/curlimages/curl:8.22.0 AS kg-release
ARG IGGYTOP_REPO=biocypher/iggytop
ARG IGGYTOP_RELEASE_TAG=latest

# The image's default `curl` user can't write at /, and this throwaway stage contributes
# nothing but /kg to the final image.
USER root

RUN if [ "$IGGYTOP_RELEASE_TAG" = "latest" ]; then \
      asset_url="https://github.com/${IGGYTOP_REPO}/releases/latest/download/knowledge_graph.tar.gz"; \
    else \
      asset_url="https://github.com/${IGGYTOP_REPO}/releases/download/${IGGYTOP_RELEASE_TAG}/knowledge_graph.tar.gz"; \
    fi && \
    curl -fsSL "$asset_url" -o /tmp/knowledge_graph.tar.gz && \
    mkdir -p /kg && \
    tar -xzf /tmp/knowledge_graph.tar.gz -C /kg --strip-components=1 && \
    rm /tmp/knowledge_graph.tar.gz

# ---- KG source "local": a knowledge_graph directory produced by create_knowledge_graph.py.
# Taken from a *named build context* rather than the main one, so it can sit anywhere on the host
# -- including under .dockerignore's `biocypher-*`, where the default cache dir puts it. ----
FROM scratch AS kg-local
COPY --from=local-kg . /kg/

# Resolves to `kg-release` or `kg-local`. BuildKit builds only the selected one, so a release
# build neither needs the `local-kg` context nor downloads anything for the unused branch.
FROM kg-${IGGYTOP_KG_SOURCE} AS kg

# ---- The image itself: neo4j plus the knowledge graph at /kg, which docker/import.sh loads on
# first start. /kg is deliberately outside $NEO4J_HOME; see that script for why. ----
FROM base AS deploy-stage
COPY --from=kg /kg /kg

# Plugin and procedure defaults, baked in so a plain `docker run` of the published image behaves
# like `docker compose up`. Any of them can still be overridden with `-e` or an env file. Plugins
# are fetched by the neo4j entrypoint on container start, which needs network access.
ENV NEO4J_PLUGINS='["graph-data-science", "apoc"]' \
    NEO4J_dbms_security_procedures_unrestricted='gds.*,apoc.*' \
    NEO4J_dbms_security_procedures_allowlist='gds.*,apoc.*' \
    NEO4J_apoc_export_file_enabled=true \
    NEO4J_apoc_import_file_enabled=true

# Catch an empty or mis-pointed /kg here rather than at the container's first start.
RUN test -f /kg/neo4j-admin-import-call.sh || { \
      echo "ERROR: /kg contains no neo4j-admin-import-call.sh."; \
      echo "  With IGGYTOP_KG_SOURCE=local, check that the local-kg build context points at a"; \
      echo "  'knowledge_graph' directory written by create_knowledge_graph.py."; \
      exit 1; \
    }
