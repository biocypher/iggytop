FROM python:3.11-slim AS setup-stage

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    autoconf \
    automake \
    libtool \
    m4 \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Sets default for BIOCYPHER_CONFIG
ARG BIOCYPHER_CONFIG=src/iggytop/config/biocypher_docker_config.yaml
ENV USED_BIOCYPHER_CONFIG=$BIOCYPHER_CONFIG

WORKDIR /usr/app/
COPY pyproject.toml uv.lock ./

# Install dependencies
RUN uv sync --frozen --no-dev --no-install-project

COPY . ./

# Install the project and set the environment path
RUN uv sync --frozen --no-dev
ENV PATH="/usr/app/.venv/bin:$PATH"

RUN cp ${USED_BIOCYPHER_CONFIG} src/iggytop/config/biocypher_config.yaml
RUN uv run create_knowledge_graph_docker.py

FROM docker.io/neo4j:4.4-enterprise AS deploy-stage
COPY --from=setup-stage /usr/app/biocypher-out/ /var/lib/neo4j/import/
COPY docker/* ./
RUN cat biocypher_entrypoint_patch.sh | cat - /startup/docker-entrypoint.sh > docker-entrypoint.sh && \
    mv docker-entrypoint.sh /startup/ && \
    chmod +x /startup/docker-entrypoint.sh
