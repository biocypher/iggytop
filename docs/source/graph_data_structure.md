# Knowledge Graph Data structure

Iggytop is built on top of BioCypher by providing a set of adapters as well as an ontology to generate knowledge graphs for TCR-epitope datasets.

The outputs aim to be compatible with the [AIRR standards](https://docs.airr-community.org/en/latest/datarep/rearrangements.html).

We aim to be integrated into scirpy, which offers a standardized way of analyzing T cell receptor (TCR) or B cell receptor (BCR) repertoires.

## Data Generation Process
The generation process relies on the `create_knowledge_graph.py` script. The pipeline follows these key steps:

### 1. Source Data Harmonization
Similar to the [tabular data structure](./tabular_data_structure.md), IggyTop leverages **BioCypher adapters** to read and harmonize data:
- **Mapping**: Source formats are mapped to internal registry keys.
- **Gene Normalization**: V(D)J genes are aligned with IMGT standards.
- **Sequence Processing**: Harmonization of CDR3 and epitope sequences.

### 2. Graph Construction
Instead of just stacking tables, the pipeline uses the BioCypher framework to:
- Instantiate **Nodes** for sequences (TRA, TRB, IGH, IGL) and epitopes based on the [Ontology](#ontology).
- Create **Edges** representing the associations between receptors and epitopes.
- The resulting graph can be exported to various formats (Neo4j, NetworkX, GraphML).

### 3. Processing Options
Users can customize the graph generation using several flags in `create_knowledge_graph.py`:

- **Receptor Types**: Specify which receptors to include (e.g., `--receptors TCR BCR`).
- **Adapter Selection**: Choose specific databases to include in the graph.
- **10X Data Filtering**: To address concerns regarding the confidence of some large-scale datasets, users can use the `--filter-10x` flag to exclude data originating from the [10X Genomics dataset](https://www.10xgenomics.com/library/a14cde). This will remove records stored in the source databases which stem from this dataset. This flag is also set for the released dataset (deduplicated_anndata.h5ad)

    **Note** The ITRAP dataset contains data from this dataset. The ITRAP data are the (5k out of 60k) pairs that have passed the ITRAP qc filtering and are therefore considered high quality. These records are not filtered out. If you want to completely exclude 10X data, consider excluding ITRAP from the pipeline.

## Design Choices
(ontology)=
### Ontology

BioCypher uses the Biolink ontology and allows custom modifications. This is done using configuration files.
The ontology used for iggytop is defined in `config/schema_config.yaml`. This includes defining the node and edge types and their relationships (hierarchy).

```text
entity
├── association
│   ├── alpha sequence to beta sequence association
│   ├── b cell receptor sequence to epitope association
│   ├── heavy sequence to light sequence association
│   └── t cell receptor sequence to epitope association
└── named thing
    └── biological entity
        └── polypeptide
            ├── epitope
            └── immune receptor sequence
                ├── b cell receptor sequence
                │   ├── igh sequence
                │   └── igl sequence
                └── t cell receptor sequence
                    ├── tra sequence
                    └── trb sequence
```


### Node and Edge Types
#### Nodes
- tra sequence
- trb sequence
- igh sequence
- igl sequence
- epitope

#### Edges
- alpha sequence to beta sequence association
- heavy sequence to light sequence association
- t cell receptor sequence to epitope association
- b cell receptor sequence to epitope association

(uniquenes)=
### Uniqueness

Immune receptor sequences are represented as nodes labeled according to their type (`tra`, `trb`, `igh`, `igl`): CDR3 sequence: and if available their V gene (see [base-adapter](./generated/iggytop.adapters.base_adapter.BaseAdapter.rst)).

Example node ID: `trb:CASSFTDTQYF:TRBV6-2`

Epitopes are represented as nodes labeled according to their type (`epitope`): (`iedb:` IRI if available or `seq:` amino acid sequence else) see [harmonize_sequences()](./generated/iggytop.adapters.utils.rst). The IRIs are retrieved using the [IEDB Database Query API](https://help.iedb.org/hc/en-us/articles/4402872882189-Immune-Epitope-Database-Query-API-IQ-API#h_01F8C6C8SN9CDMBWQ41MWES31A), see [get_iedb_ids_batch()](./generated/iggytop.adapters.utils.rst).

Example node IDs: `epitope:iedb:37257`, `epitope:seq:SLSNRLYYL`

Edges link between two nodes; their ID is: source node - target node ID.

Example edges: `tra:CAVTTDSWGKLQF:TRAV12-2-trb:CASRPGLAGGRPEQYF:TRBV6-5`, `tra:CAVTTDSWGKLQF:TRAV12-2-epitope:iedb:37257`

## Output Formats and Availability

The knowledge graph can be exported in several ways:
1. **Neo4j**: Optimized for graph database queries. Check out the Docker guide in the README.
2. **NetworkX / GraphML**: Useful for Python-based graph analysis and visualization in tools like Cytoscape.
3. **AIRR JSON**: While natively a graph, output can be converted back to the AIRR format (tabular).

### Bimonthly Releases
Knowledge graph exports (e.g., in GraphML) are not yet provided in bimonthly releases.

## Creating Your Own Graph

You can run the graph generation locally to create custom subsets or use specific versions of the data:

```bash
python create_knowledge_graph.py --adapters VDJDB MCPAS --filter-10x
```
Note that some parameters are defined in the `config/biocypher_config.yaml`. Check out this file and change it for more control (eg defining output type).

### Assumptions

During construction of the graph, redundant data can be neglected (e.g., pairs reported in multiple databases); however, some information is also lost.
See [this issue](https://github.com/biocypher/iggytop/issues/31).
