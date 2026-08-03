FROM docker.io/neo4j:4.4-enterprise AS deploy-stage

# Bake the knowledge graph from an IggyTop data release into the image, so `docker compose up`
# needs neither a local run of the release pipeline nor a host-mounted ./knowledge_graph directory.
ARG IGGYTOP_REPO=biocypher/iggytop
ARG IGGYTOP_RELEASE_TAG=latest

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl ca-certificates && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir -p /kg && \
    if [ "$IGGYTOP_RELEASE_TAG" = "latest" ]; then \
      asset_url="https://github.com/${IGGYTOP_REPO}/releases/latest/download/knowledge_graph.tar.gz"; \
    else \
      asset_url="https://github.com/${IGGYTOP_REPO}/releases/download/${IGGYTOP_RELEASE_TAG}/knowledge_graph.tar.gz"; \
    fi && \
    curl -fsSL "$asset_url" -o /tmp/knowledge_graph.tar.gz && \
    tar -xzf /tmp/knowledge_graph.tar.gz -C /kg --strip-components=1 && \
    rm /tmp/knowledge_graph.tar.gz

COPY docker/* ./
RUN cat entrypoint_patch.sh | cat - /startup/docker-entrypoint.sh > docker-entrypoint.sh && \
    mv docker-entrypoint.sh /startup/ && \
    chmod +x /startup/docker-entrypoint.sh
