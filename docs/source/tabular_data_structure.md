# Tabular Data Structure

IggyTop provides a pipeline to generate harmonized, tabular datasets of immunoreceptor-epitope pairings. Instead of building a knowledge graph, this approach stacks harmonized tables from the adapters into a single large dataset, following the [AIRR standards](https://docs.airr-community.org/en/stable/datarep/rearrangements.html).

## Data Generation Process
The generation process relies on the `create_anndata.py` script. The pipeline follows these key steps:

### 1. Source Data Harmonization
IggyTop leverages **BioCypher adapters** to read data from multiple immunological databases.
- **Mapping**: Each database's unique format is mapped to a set of internal registry keys.
- **Gene Normalization**: V(D)J genes are normalized to IMGT standards.
- **Sequence Processing**: CDR3 sequences and epitopes are harmonized.
for details on the harmonization please refer to the `self.read_table` function of any given adapter.

### 2. Tabular Stacking
Once each source is processed into a harmonized format:
- The individual tables (in anndata format) are merged ("stacked") into a single, unified dataset.
- The resulting dataset is stored as an `AnnData` object, which is highly compatible with the Python-based single-cell ecosystem (e.g., [Scirpy](https://scirpy.readthedocs.io/)).
- We also provide the same data in json.gz format for interoperability (eg with R).

### 3. Processing Options

Users can customize the dataset generation using several flags in `create_anndata.py`:

- **Receptor Types**: Specify which receptors to include (e.g., `--receptors TCR BCR`).
- **Adapter Selection**: Choose specific databases to include in the graph.
- **Deduplication**: By default, the script can deduplicate entries across different source databases based on matching CDR3 sequences and epitopes. This can be toggled using `--not_deduplicate`.

- **10X Data Filtering**: To address concerns regarding the confidence of some large-scale datasets, users can use the `--filter-10x` flag to exclude data originating from the [10X Genomics dataset](https://www.10xgenomics.com/library/a14cde). This will remove records stored in the source databases which stem from this dataset. This flag is also set for the released dataset (deduplicated_anndata.h5ad)

    **Note** The ITRAP dataset contains data from this dataset. The ITRAP data are the (5k out of 60k) pairs that have passed the ITRAP qc filtering and are therefore considered high quality. These records are not filtered out. If you want to completely exclude 10X data, consider excluding ITRAP from the pipeline.

## Output Formats and Availability

IggyTop datasets are made available in two main formats:
1. **AnnData (.h5ad)**: Optimized for analysis with Scanpy and Scirpy.
2. **AIRR JSON (.json.gz)**: A standardized, compressed format for cross-tool compatibility.

### Bimonthly Releases
We provide pre-generated, bimonthly releases of these datasets. These releases include:
- A harmonized, comprehensive merged dataset (all sources). `merged`
- A additionally deduplicated, and filtered (10X) "IggyTop" dataset.

### Technical Applications: Scirpy Integration
One major application of this tabular pipeline is its potential integration into **Scirpy**. By providing a standardized data retrieval option, IggyTop enables researchers to easily download and use these integrated reference databases directly within their TCR/BCR analysis workflows.

## Creating Your Own Dataset

You can run the pipeline locally to create custom subsets or use specific versions of the data:

```bash
python create_anndata.py --adapters VDJDB MCPAS --filter-10x
```
