### IggyTop Data Release

IggyTop (**I**mmunolo**g**ical **G**raph **Y**ielding **Top** receptor-epitope pairings) harmonizes and integrates various immunoreceptor-epitope databases using the BioCypher framework. This release includes updated data from multiple sources in AnnData and AIRR formats.

For more information, please visit the [IggyTop GitHub repository](https://github.com/biocypher/iggytop).

#### Output Datasets

This release includes the following datasets:

Datasets are provided in AnnData (.h5ad) and AIRR JSON (.json.gz) formats.

- **Merged Dataset**:

    A comprehensive dataset combining all sources. Harmonized, no deduplication applied.

    `merged_anndata.h5ad`

    `merged_airr_cells.json.gz`

- **Deduplicated Dataset**: Considered the IggyTop dataset.

    This corresponds to the merged dataset after some filtering steps:
      - Data originating from the 10X dataset is filtered out (except data in the ITRAP dataset).
      - Records are deduplicated (see `create_anndata.py`)

    Find out more about the deduplication and filtering in the [IggyTop documentation](https://iggytop.readthedocs.io/en/latest/).

    `deduplicated_anndata.h5ad`

    `deduplicated_airr_cells.json.gz`

#### Data Source Information

Both datasets (merged and deduplicated) use the same source datasets, reported here. The full source
information is available in `metadata.json`.

{{SOURCE_TABLE}}
