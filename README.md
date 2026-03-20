# IggyTop: **I**mmunolo**g**ical **G**raph **Y**ielding **Top** receptor-epitope pairings

[![Python Version](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

![figure1](./overview.png)

This repository uses [BioCypher](https://biocypher.org) framework for harmonization of databases with existing immunoreceptor-epitope matching information.

BioCypher is designed to facilitate the standardized integration of heterogeneous data sources through a regulated framework. The BioCypher framework implements a modular architecture where each data source is processed through dedicated transformation scripts called adapters. These adapters serve as the primary interface between raw data sources and the BioCypher knowledge graph infrastructure. This project provides adapters for the following databases:

- [IEDB](https://www.iedb.org/)
- [VDJdb](https://github.com/antigenomics/vdjdb-db)
- [McPAS-TCR](https://friedmanlab.weizmann.ac.il/McPAS-TCR/)
- [CEDAR](https://cedar.iedb.org/home_v3.php)
- [ITRAP](https://github.com/mnielLab/ITRAP_benchmark/blob/main/ITRAP.csv) (data filtered from [10X Genomics Dataset](https://www.10xgenomics.com/library/a14cde))
- [TRAIT](https://pgx.zju.edu.cn/traitdb/)
- [TCR3d](https://tcr3d.ibbr.umd.edu/)
- [NeoTCR](https://github.com/lyotvincent/NeoTCR?tab=readme-ov-file)

These include data from both, original sources, extracting data directly from studies, such es McPAS-TCR, and from already pulled sources such as TRAIT.
A script is provided to build a knowledge graph with all these adapters. On a consumer laptop, building the full graph typically takes 20-30 mins.

The final output is the **IggyTop** database, which integrates immunoreceptor-epitope matching information from all supported data sources in the unified list of [AIRR cells](https://scirpy.scverse.org/en/stable/generated/scirpy.io.AirrCell.html).

## Node and Edge Types
### Nodes
- tra sequence
- trb sequence
- igh sequence
- igl sequence
- epitope

### Edges
- alpha sequence to beta sequence association
- heavy sequence to light sequence association
- t cell receptor sequence to epitope association
- b cell receptor sequence to epitope association

## Prerequisites

- [uv](https://docs.astral.sh/uv/) for dependency management
- [docker](https://www.docker.com/get-started/) optional for containerization

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
   uv run create_anndata.py
   ```
More information can be found in the documentation (see below).
## Pipeline

- `create_knowledge_graph.py`: the main script that orchestrates the pipeline.
It brings together the BioCypher package with the data sources. It calls the `io.create_knowledge_graph()` function which creates a knowledge graph including all available databases and saves it to airr format in a json file. use the `--adapters` flag to select single source databases
```bash
uv run create_knowledge_graph.py --adapters VDJDB CEDAR
```
- `create_anndata.py`: this script can be used to obtain the harmonized, merged (and deduplicated) data from all (or selected) available databases in [anndata](https://anndata.readthedocs.io/en/stable/index.html) format.
It will initialize the adapters but not generate the knowledge graph. The main purpose is integration of the available data into [Scirpy](https://scirpy.scverse.org/en/latest/). You can specify which adapters to include:
```bash
uv run create_anndata.py --adapters VDJDB CEDAR
```

- `src/iggytop/adapters` contains modules that define the adapter to the data source.

- `src/iggytop/config/schema_config.yaml`: a configuration file
that defines the schema of the knowledge graph. It is used by BioCypher to map
the data source to the knowledge representation on the basis of ontology (see
[this part of the BioCypher tutorial](https://biocypher.org/tutorial-ontology.html)).

- `src/iggytop/config/biocypher_config.yaml`: a configuration file that defines some BioCypher parameters, such as the mode, the
separators used, and other options. More on its use can be found in the
[Documentation](https://biocypher.org/BioCypher/reference/biocypher-config/).

## Documentation

This repository uses [Sphinx](https://www.sphinx-doc.org/) for documentation.

### Building the Documentation

To build the documentation, ensure you have the `docs` dependency group installed:

```bash
uv sync --group docs
```

Then, execute the following command:

```bash
uv run ./update_docs.sh
```

This will generate the documentation in the `docs/build` directory.

### Hosting the Documentation Locally

To host the documentation locally, run:

```bash
uv run python3 -m http.server --directory docs/build 8000
```

You can then access the documentation in your browser at `http://localhost:8000`.

## Testing and CI/CD

IggyTop uses GitHub Actions to automate **bimonthly data releases** and ensure data integrity through continuous integration.

### Bimonthly Data Releases
- **Frequency**: Automated releases on the **1st day of every 2nd month**. (first scheduled on **May 1,2026**)
- **Release Assets**:
    - `deduplicated_anndata.h5ad`: Harmonized dataset for [Scirpy](https://scirpy.scverse.org/).
    - `deduplicated_airr_cells.json.gz`: Standardized AIRR JSON export.
    - `metadata.json`: Source versions, download dates, SHA-256 checksums, and change indicators.
    - `RELEASE_NOTES.md`: A human-readable record of data source changes.

### Automated Testing
Before any data is released, the CI pipeline runs a validation suite to catch breaking changes in upstream databases.

**How to run tests locally:**
```bash
uv run pytest tests/
```

### Local Dev Workflows (using `act`)
To simulate the ci pipeline locally without pushing to GitHub, you can use **[act](https://github.com/nektos/act)**. Ensure you have Docker and `act` installed, then run:

```bash
# Run the workflow
act workflow_dispatch -W .github/workflows/ci_ingestion.yml
```

### Metadata & Verification
Each release contains a `metadata.json` file. This metadata is also embedded within the data files:
- **AnnData**: Access via `adata.uns["iggytop_metadata"]`.
- **JSON**: Found in the top-level `"metadata"` key of the `.json.gz` file.
The metadata includes a `has_changed` flag (⚠️) for each source, which compares the currently downloaded version with the version from the previous GitHub release.

## 🐳 Docker

This repo also contains a `docker compose` workflow to create the example
database using BioCypher and load it into a dockerised Neo4j instance
automatically. To run it, simply execute
```
docker compose up -d --build
```
in the root directory of the project. The example instance consists of the TCR3d database only as it is small enough to visualize, for other database compositions, just edit the `create_knowledge_graph_docker.py` script to your needs. This will start up a single (detached) docker
container with a Neo4j instance that contains the knowledge graph built by
BioCypher as the DB `docker`, which you can connect to and browse at
localhost:7474. Authentication is set to `neo4j/neo4jpassword` by default
and can be modified in the `docker_variables.env` file.

Open http://localhost:7474 to access the neo4j database. You can now run queries against the database.
To get a visual representation of the tcr3d knowledge grraph constructed by iggytop, run the following CYPHER query:
```
MATCH (n) return n
```

The `biocypher_docker_config.yaml` file is used instead of the
`biocypher_config.yaml`. Everything else is the same as in the local setup. The
first container installs and runs the BioCypher pipeline, and the second
container installs and runs Neo4j. The files created by BioCypher in the first
container are copied and automatically imported into the DB in the second
container.

## Scirpy integration

This project helps generating the anndata versions of all the Scirpy reference databases supported by Iggytop. The Anndata objects are stored in h5ad file format.
This can be replicated by running the `create_anndata` script (while selecting the databases of interest using the variable `adapters_to_include`)

```bash
uv run create_anndata.py
```

Note: this is a wip

## Related work:
If you find a dataset (eg training data for a model) and would like to find the source of the records using the IggyTop dataset, check out [this tool](https://github.com/RaphaelDeGottardi/TCR_source_detection).

This project is built on (and part of) the [BioCypher](https://github.com/biocypher/biocypher) framework. Make sure to check out this cool project!


## Contributing

Contributions are welcome! Please feel free to submit a Pull Request or create
an Issue if you discover any problems.

## License

This project is licensed under the MIT License - see the LICENSE file for
details.
