# BioCypher project template
A quick way to set up a BioCypher-driven knowledge graph pipeline.

## Using the GitHub Template functionality
You can use this template in GitHub directly. Just select 
`biocypher/project-template` as your template when creating a new repository
on GitHub.

## ⚙️ Installation (local, for docker see below)
1. Clone this repository and rename to your project name.
```{bash}
git clone https://github.com/biocypher/project-template.git
mv project-template my-project
cd my-project
```
2. Make the repository your own.
```{bash}
rm -rf .git
git init
git add .
git commit -m "Initial commit"
# (you can add your remote repository here)
```
3. Install the dependencies using [Poetry](https://python-poetry.org/). (Or feel
free to use your own dependency management system. We provide a `pyproject.toml`
to define dependencies.)
```{bash}
poetry install
```
4. You are ready to go!
```{bash}
poetry shell
python create_knowledge_graph.py
```


## 🛠 Usage

### Structure
The project template is structured as follows:
```
.
│  # Project setup
│
├── LICENSE
├── README.md
├── pyproject.toml
│
│  # Docker setup
│
├── Dockerfile
├── docker
│   ├── biocypher_entrypoint_patch.sh
│   ├── create_table.sh
│   └── import.sh
├── docker-compose.yml
├── docker-variables.env
│
│  # Project pipeline
│
├── create_knowledge_graph.py
├── config
│   ├── biocypher_config.yaml
│   ├── biocypher_docker_config.yaml
│   └── schema_config.yaml
└── template_package
    └── adapters
        └── example_adapter.py
```

The main components of the BioCypher pipeline are the
`create_knowledge_graph.py`, the configuration in the `config` directory, and
the adapter module in the `template_package` directory. The latter can be used
to publish your own adapters (see below). You can also use other adapters from
anywhere on GitHub, PyPI, or your local machine.

**The BioCypher ecosystem relies on the collection of adapters (planned, in
development, or already available) to inform the community about the available
data sources and to facilitate the creation of knowledge graphs. If you think
your adapter could be useful for others, please create an issue for it on the
[main BioCypher repository](https://github.com/biocypher/biocypher/issues).**

In addition, the docker setup is provided to run the pipeline (from the same
python script) in a docker container, and subsequently load the knowledge graph
into a Neo4j instance (also from a docker container). This is useful if you want
to run the pipeline on a server, or if you want to run it in a reproducible
environment.

### Running the pipeline

`python create_knowledge_graph.py` will create a knowledge graph from the
example data included in this repository (borrowed from the [BioCypher
tutorial](https://biocypher.org/tutorial.html)). To do that, it uses the
following components:

- `create_knowledge_graph.py`: the main script that orchestrates the pipeline.
It brings together the BioCypher package with the data sources. To build a 
knowledge graph, you need at least one adapter (see below). For common 
resources, there may already be an adapter available in the BioCypher package or
in a separate repository. You can also write your own adapter, should none be
available for your data.

- `example_adapter.py` (in `template_package.adapters`): a module that defines
the adapter to the data source. In this case, it is a random generator script.
If you want to create your own adapters, we recommend to use the example adapter
as a blueprint and create one python file per data source, approproately named.
You can then import the adapter in `create_knowledge_graph.py` and add it to
the pipeline. This way, you ensure that others can easily install and use your 
adapters.

- `schema_config.yaml`: a configuration file (found in the `config` directory)
that defines the schema of the knowledge graph. It is used by BioCypher to map
the data source to the knowledge representation on the basis of ontology (see
[this part of the BioCypher 
tutorial](https://biocypher.org/tutorial-ontology.html)).

- `biocypher_config.yaml`: a configuration file (found in the `config` 
directory) that defines some BioCypher parameters, such as the mode, the 
separators used, and other options. More on its use can be found in the
[Documentation](https://biocypher.org/installation.html#configuration).

### Publishing your own adapters
After adding your adapter(s) to the `adapters` directory, you may want to
publish them for easier reuse. To create a package to distribute your own
adapter(s), we recommend using [Poetry](https://python-poetry.org/). Poetry,
after setup, allows you to publish your package to PyPI using few simple
commands. To set up your package, rename the `template_package` directory to
your desired package name and update the `pyproject.toml` file accordingly. Most
importantly, update the `name`,`author`, and `version` fields. You can also add
a `description` and a `license`.  Then, you can publish your package to PyPI
using the following commands:

```{bash}
poetry build
poetry publish
```

If you don't want to publish your package to PyPI, you can also install it from
GitHub using the respective functions of poetry or pip.

## 🐳 Docker

This repo also contains a `docker compose` workflow to create the example
database using BioCypher and load it into a dockerised Neo4j instance
automatically. To run it, simply execute `docker compose up -d` in the root 
directory of the project. This will start up a single (detached) docker
container with a Neo4j instance that contains the knowledge graph built by
BioCypher as the DB `docker`, which you can connect to and browse at 
localhost:7474 (don't forget to switch the DB to `docker` instead of the 
standard `neo4j`). Authentication is set to `neo4j/neo4jpassword` by default
and can be modified in the `docker_variables.env` file.

By using the `BIOCYPHER_CONFIG` environment variable in the Dockerfile, the
`biocypher_docker_config.yaml` file is used instead of the 
`biocypher_config.yaml`. Everything else is the same as in the local setup. The
first container installs and runs the BioCypher pipeline, and the second
container installs and runs Neo4j. The files created by BioCypher in the first
container are copied and automatically imported into the DB in the second
container.
