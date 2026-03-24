import gzip
import json
import os

import pytest
from iggytop.adapters.utils import deduplicate_and_aggregate, get_file_checksum, save_airr_cells_json
from scirpy.io._datastructures import AirrCell


def test_get_file_checksum(tmp_path):
    test_file = tmp_path / "test.txt"
    test_file.write_text("hello world")

    checksum = get_file_checksum(str(test_file))
    # SHA256 for "hello world"
    expected = "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
    assert checksum == expected


def test_save_airr_cells_json_with_metadata(tmp_path):
    # Mock data
    cell = AirrCell("id1")
    cells = [cell]
    metadata = {"test": "data"}

    save_airr_cells_json(cells, directory=str(tmp_path), filename="test_output", metadata=metadata)

    output_path = tmp_path / "test_output.json.gz"
    assert os.path.exists(output_path)

    with gzip.open(output_path, "rt") as f:
        data = json.load(f)
        assert data["metadata"] == metadata
        assert len(data["cells"]) == 1
        assert data["cells"][0]["cell_id"] == "id1"


def test_deduplication_on_real_data_subset():
    """Test deduplication logic using the real data subset created for testing."""
    import anndata as ad

    test_subset_path = "tests/data/test_subset.h5ad"
    if not os.path.exists(test_subset_path):
        pytest.skip("Test subset data not found.")

    adata = ad.read_h5ad(test_subset_path)

    # We'll duplicate the first 5 rows to ensure we have something to deduplicate
    original_n_obs = adata.n_obs
    dup_adata = adata[:5].copy()
    # Modify PMIDs in duplicates to test aggregation
    dup_adata.obs["PMID"] = "NEW_PMID"

    combined_adata = ad.concat([adata, dup_adata])
    combined_adata.obs_names_make_unique()  # Ensure unique obs names for the duplicate

    subset_cols = ["VJ_1_junction_aa", "VJ_1_v_call", "VDJ_1_v_call", "VDJ_1_junction_aa", "iedb_iri"]
    agg_cols = ["PMID", "source"]

    # Run deduplication
    dedup = deduplicate_and_aggregate(combined_adata, subset_cols, agg_cols)

    # The number of observations should be exactly what we started with
    assert dedup.n_obs == original_n_obs

    # Check if aggregation worked for one of the duplicates
    # We check if the NEW_PMID was joined with the original PMID
    for i in range(5):
        pmid_val = dedup.obs["PMID"].iloc[i]
        assert "NEW_PMID" in pmid_val
