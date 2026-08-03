# IggyTop: **I**mmunolo**g**ical **G**raph **Y**ielding **Top** receptor-epitope pairings

[![Python Version](https://img.shields.io/badge/python-3.13+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![docs](https://img.shields.io/readthedocs/iggytop)](https://img.shields.io/readthedocs/iggytop)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20037575.svg)](https://zenodo.org/records/20037575)
![figure1](./overview.png)

This repository uses the [BioCypher](https://biocypher.org) framework to harmonize databases with existing immunoreceptor-epitope matching information. The aggregated data is provided in tabular and graph formats.

BioCypher provides a modular architecture where each data source is processed through dedicated transformation scripts called adapters. These adapters are the interface between raw data sources and the IggyTop dataset. We use [tidytcells](https://tidytcells.readthedocs.io/en/stable/index.html) to standardize gene names, cdr3 sequences and junction sequences for both BCR and TCR datasets. We aim to follow [AIRR](https://docs.airr-community.org/en/stable/) and [IMGT](https://www.imgt.org/) standards. This project provides adapters for the following databases:

- [IEDB](https://www.iedb.org/)
- [VDJdb](https://github.com/antigenomics/vdjdb-db)
- [McPAS-TCR](https://friedmanlab.weizmann.ac.il/McPAS-TCR/)
- [CEDAR](https://cedar.iedb.org/home_v3.php)
- [ITRAP](https://github.com/mnielLab/ITRAP_benchmark/blob/main/ITRAP.csv) (data filtered from [10X Genomics Dataset](https://www.10xgenomics.com/library/a14cde))
- [TRAIT](https://pgx.zju.edu.cn/traitdb/)
- [TCR3d](https://tcr3d.ibbr.umd.edu/)
- [NeoTCR](https://github.com/lyotvincent/NeoTCR?tab=readme-ov-file)
- [BATCAVE](https://github.com/meyer-lab-cshl/BATMAN-paper)

The aggregated data from all adapters is available in the (bimonthly) [releases](https://github.com/biocypher/iggytop/releases). These releases are used by [Scirpy](https://scirpy.scverse.org/en/v0.24.0/index.html), which provides a simple interface to access the data in [anndata](https://anndata.readthedocs.io/en/latest/) format. You can also rebuild the DB using the provided code with custom parameters. On a consumer laptop, building the full DB typically takes 20-30 minutes.

## Quick start (Scirpy)

To use the database in Python, install [Scirpy](https://scirpy.scverse.org/en/v0.24.0/index.html) (>=v0.24.0). Then use [this function](https://scirpy.scverse.org/en/v0.24.0/generated/scirpy.datasets.iggytop.html) to load the dataset in [anndata](https://anndata.readthedocs.io/en/latest/) (tabular) format:

```python
import scirpy as ir

iggytop = ir.datasets.iggytop(deduplicated=True, tag='latest')

# Or (e.g. for VDJdb only)
vdjdb = scirpy.datasets.vdjdb(tag='latest')
```


## Graphs vs. Tables
Two paths are covered: A tabular path that stacks source databases into a large table (used in Scirpy and Releases), and a knowledge graph path that converts the source data into a graph. Both are documented in the [documentation](https://iggytop.readthedocs.io/en/latest/). For details, see [Graph Data Structure](https://iggytop.readthedocs.io/en/latest/graph_data_structure.html) and [Tabular Data Structure](https://iggytop.readthedocs.io/en/latest/tabular_data_structure.html).

## Prerequisites
To run DB/graph generation locally:
- [uv](https://docs.astral.sh/uv/): for dependency management
- [docker](https://www.docker.com/get-started/): optional for Neo4j (see below)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/biocypher/iggytop.git
   cd iggytop
   ```

2. Install dependencies using uv:
   ```bash
   # Core installation (includes dev dependencies)
   uv sync

   # Include documentation and Jupyter tools
   uv sync --group docs
   ```

3. You are ready to go!
   ```bash
   uv run create_knowledge_graph.py
   ```
   or

   ```bash
   uv run create_release.py
   ```
More information can be found in the [documentation](https://iggytop.readthedocs.io/en/latest/).

## Documentation
We use [Sphinx](https://www.sphinx-doc.org/en/master/) for documentation (see `./docs`).
The full documentation is available online via [Read the Docs](https://iggytop.readthedocs.io/en/latest/).

## Pipeline

- `create_release.py`: builds every selected adapter's source table exactly once, then emits the raw merged dataset, the filtered+deduplicated dataset (in [anndata](https://anndata.readthedocs.io/en/stable/index.html)/ AIRR JSON format, for integration into [Scirpy](https://scirpy.scverse.org/en/latest/)), and the neo4j knowledge graph files from that same run, avoiding rebuilding the tables multiple times. Each output can be toggled independently:
```bash
uv run create_release.py --adapters VDJDB CEDAR --filter-10x
```
`--filter-10x` only affects the deduplicated dataset. Use `--not_merged`, `--not_deduplicate` or `--not_graph` to skip any of the three outputs, and `--graph-output-format` to choose the graph's BioCypher output format (default `neo4j`).

- `create_knowledge_graph.py`: a standalone, graph-only entrypoint into `io.create_knowledge_graph()`, useful when you don't need the AnnData/AIRR outputs. Use the `--adapters` flag to select specific source databases and `--output-format` to choose the output (`neo4j`, `networkx` or `airr`; default `neo4j`):
```bash
uv run create_knowledge_graph.py --adapters VDJDB CEDAR --output-format networkx
```

- `src/iggytop/adapters` contains modules that define each data source adapter.

- `src/iggytop/config/schema_config.yaml`: defines the schema of the knowledge graph. It is used by BioCypher to map the data source to the knowledge representation based on ontology (see [this part of the BioCypher tutorial](https://biocypher.org/tutorial-ontology.html)).

- `src/iggytop/config/biocypher_config.yaml`: defines BioCypher parameters such as the mode, separators, and other options. More on its use can be found in the [documentation](https://biocypher.org/BioCypher/reference/biocypher-config/).

## Releases, Testing and CI/CD

IggyTop uses GitHub Actions to automate **bimonthly data releases** and ensure data integrity through continuous integration.
Currently this only involves the tabular part of IggyTop (`create_release.py`).
Check out the latest release [here](https://github.com/biocypher/iggytop/releases/latest).
### Bimonthly Data Releases
Find the releases [here](https://github.com/biocypher/iggytop/releases)
- **Frequency**: Automated releases on the **1st day of every 2nd month** (first scheduled on **May 1, 2026**).
- **Release Assets**:
   Check out the release notes for more information on the released datasets.

### Automated Testing
Before any data is released, the CI pipeline (based on GitHub Actions) runs a validation suite to catch breaking changes in upstream databases.

**How to run tests locally:**
```bash
# Install all dependencies (including docs for notebook testing)
uv sync --all-groups

# Install Jupyter kernel for notebook execution (one-time setup)
uv run python -m ipykernel install --user --name python3

# Run all tests (including notebook validation)
uv run pytest tests/
```
**Why the kernel installation?** The test suite includes validation of Jupyter notebooks (tutorials and database summaries) to ensure they execute without errors. This requires a Jupyter kernel registered with the name "python3" to match the notebooks' configuration. The installation is a one-time setup per environment.

**How to test the CI pipelines**
Ensure you have Docker and `act` installed, then run:
```bash
# Run the workflow
act workflow_dispatch -W .github/workflows/ci_ingestion.yml
```

## Graph visualization using Neo4j on Docker

This repo contains a `docker compose` workflow that spins up a Neo4j
instance pre-loaded with the full IggyTop knowledge graph from a
[data release](https://github.com/biocypher/iggytop/releases), no local
Python/BioCypher build required. To run it, execute
```bash
docker compose up -d --build
```
in the root directory of the project.

At build time, the image downloads `knowledge_graph.tar.gz` (the neo4j-admin
import CSVs + import script produced by `create_release.py`) from a pinned
IggyTop release and bakes it into the image. On first container start, it
runs a `neo4j-admin database import` — the fast bulk-loader Neo4j recommends
for populating an empty database — directly against those CSVs, rather than
loading the graph row-by-row via Cypher. The imported data then lives in the
named volume `biocypher_neo4j_volume`, so subsequent restarts skip the
import entirely.

By default the build uses whatever release GitHub currently marks
"latest". To pin a specific release instead, override the build arg:
```bash
IGGYTOP_RELEASE_TAG=data-2026.08.03.063907 docker compose up -d --build
```
To force a re-import after switching release tags, first tear down the volume:
```bash
docker compose down -v
```

This will start up a single (detached) Docker container with a Neo4j
instance containing the full knowledge graph, which you can connect to and
browse at localhost:7474 (note that this can take a couple minutes as the graph is large). Authentication is set to `neo4j/neo4jpassword` by
default and can be modified in the `docker-variables.env` file.

Open http://localhost:7474 to access the Neo4j database. You can now run
queries against the database — the full graph has over a million nodes, so
scope exploratory queries with a `LIMIT`, e.g.:
```
MATCH (n) RETURN n LIMIT 100
```

To get a feel for the schema, a good starting point is to expand outward
from a single `Binding` node — the central node type linking a TCR/BCR
receptor to the epitope it binds — and follow its relationships a few hops
out:
```
MATCH path = (n:Binding)-[*1..4]->(m) WHERE n.complete = true RETURN path LIMIT 14
```
For more information on the Graph schema, check out [Graph Data Structure](https://iggytop.readthedocs.io/en/latest/graph_data_structure.html).
The result should look similar to this:


![knowledge graph example](./docs/source/graph.png)

## Related work
If you find a dataset (e.g. training data for a model) and would like to find the source of the records using the IggyTop dataset, check out [this tool](https://github.com/RaphaelDeGottardi/TCR_source_detection).

This project is built on (and part of) the [BioCypher](https://github.com/biocypher/biocypher) framework. Make sure to check out this cool project!


## Contributing

Contributions are welcome! Please feel free to submit a Pull Request or create
an Issue if you discover any problems.

## License

This project is licensed under the MIT License - see the LICENSE file for
details.

## Citation

If you use IggyTop in your research, please cite it using the following DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20037575.svg)](https://zenodo.org/records/20037575)

You can find the full citation details on the [Zenodo page](https://zenodo.org/records/20037575).

We also provide a CITATION.cff file for customized citations.
