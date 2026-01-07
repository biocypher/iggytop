# Data structure and principles

## Compatibility

Iggytop is built on top of BioCypher by providing a set of adapters as well as an ontology to generate knowledge graphs for TCR-epitope datasets.

The outputs aim to be compatible with the AIRR standards.

We aim to be integrated into scirpy, which offers a standardized way of analyzing T cell receptor (TCR) or B cell receptor (BCR) repertoires.

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
(uniquenes)=
### Uniqueness

Immune receptor sequences are represented as nodes labeled according to their type (`tra`, `trb`, `igh`, `igl`): CDR3 sequence: and if available their V gene (see [base-adapter](./generated/iggytop.adapters.base_adapter.BaseAdapter.rst)).

Example node ID: `trb:CASSFTDTQYF:TRBV6-2`

Epitopes are represented as nodes labeled according to their type (`epitope`): (`iedb:` IRI if available or `seq:` amino acid sequence else) see [harmonize_sequences()](./generated/iggytop.adapters.utils.rst). The IRIs are retrieved using the [IEDB Database Query API](https://help.iedb.org/hc/en-us/articles/4402872882189-Immune-Epitope-Database-Query-API-IQ-API#h_01F8C6C8SN9CDMBWQ41MWES31A), see [get_iedb_ids_batch()](./generated/iggytop.adapters.utils.rst).

Example node IDs: `epitope:iedb:37257`, `epitope:seq:SLSNRLYYL`

Edges link between two nodes; their ID is: source node - target node ID.

Example edges: `tra:CAVTTDSWGKLQF:TRAV12-2-trb:CASRPGLAGGRPEQYF:TRBV6-5`, `tra:CAVTTDSWGKLQF:TRAV12-2-epitope:iedb:37257`

### Assumptions

During construction of the graph, redundant data can be neglected (e.g., pairs reported in multiple databases); however, some information is also lost.
See [this issue](https://github.com/biocypher/iggytop/issues/31).

