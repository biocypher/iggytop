from biocypher import BioCypher
from template_package.adapters.example_adapter import (
    ExampleAdapter,
    ExampleAdapterNodeType,
    ExampleAdapterEdgeType,
    ExampleAdapterProteinField,
    ExampleAdapterDiseaseField,
)

# Instantiate the BioCypher interface
# You can use `config/biocypher_config.yaml` to configure the framework or
# supply settings via parameters below
bc = BioCypher(
    biocypher_config_path="config/biocypher_config.yaml",
)

# Take a look at the ontology structure of the KG according to the schema
bc.show_ontology_structure()

# Choose node types to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_types = [
    ExampleAdapterNodeType.PROTEIN,
    ExampleAdapterNodeType.DISEASE,
]

# Choose protein adapter fields to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_fields = [
    # Proteins
    ExampleAdapterProteinField.ID,
    ExampleAdapterProteinField.SEQUENCE,
    ExampleAdapterProteinField.DESCRIPTION,
    ExampleAdapterProteinField.TAXON,
    # Diseases
    ExampleAdapterDiseaseField.ID,
    ExampleAdapterDiseaseField.NAME,
    ExampleAdapterDiseaseField.DESCRIPTION,
]

edge_types = [
    ExampleAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION,
    ExampleAdapterEdgeType.PROTEIN_DISEASE_ASSOCIATION,
]

# Create a protein adapter instance
adapter = ExampleAdapter(
    node_types=node_types,
    node_fields=node_fields,
    edge_types=edge_types,
    # we can leave edge fields empty, defaulting to all fields in the adapter
)


# Create a knowledge graph from the adapter
bc.write_nodes(adapter.get_nodes())
bc.write_edges(adapter.get_edges())

# Write admin import statement
bc.write_import_call()

# Check output
bc.log_duplicates()
bc.log_missing_bl_types()
