import platformdirs
from pathlib import Path
import importlib.resources
from biocypher import BioCypher
import anndata as ad
import scirpy as ir


from iggytop.adapters.cedar_adapter import CEDARAdapter
from iggytop.adapters.iedb_adapter import IEDBAdapter
from iggytop.adapters.mcpas_adapter import MCPASAdapter
from iggytop.adapters.neotcr_adapter import NeoTCRAdapter
from iggytop.adapters.tcr3d_adapter import TCR3DAdapter
from iggytop.adapters.trait_adapter import TRAITAdapter
from iggytop.adapters.vdjdb_adapter import VDJDBAdapter
from iggytop.adapters.utils import _set_up_config

import os
import time

merge = True  # whether to merge all datasets into a single AnnData
deduplicate = False  # whether to deduplicate the merged AnnData, set unique id's below

start_time = time.time()


cache_dir = Path(platformdirs.user_cache_dir("iggytop_anndata_fresh"))
test_mode = False
adapters_to_include = ["VDJDB", "IEDB", "MCPAS",  "TCR3D", "NeoTCR", "CEDAR", "TRAIT"]
output_format = "airr"  # doesnt matter here

config_path = _set_up_config(output_format, cache_dir)

schema_config_path = os.path.join(cache_dir, 'schema_config.yaml')
with importlib.resources.open_text('iggytop.config', 'schema_config.yaml') as schema_file:
    with open(schema_config_path, 'w') as cache_schema_file:
        cache_schema_file.write(schema_file.read())

bc = BioCypher(biocypher_config_path=config_path,
                schema_config_path=schema_config_path,
                cache_directory=cache_dir)


adapter_classes = {
    "VDJDB": VDJDBAdapter,
    "MCPAS": MCPASAdapter,
    "TRAIT": TRAITAdapter,
    "IEDB": IEDBAdapter,
    "TCR3D": TCR3DAdapter,
    "NeoTCR": NeoTCRAdapter,
    "CEDAR": CEDARAdapter,
}


selected_adapters = [
    adapter_classes[name]
    for name in adapters_to_include
    if name in adapter_classes
]

print(f"Initialization took {time.time() - start_time:.2f} seconds")

for AdapterClass in selected_adapters:
    start_time = time.time()
    adapter = AdapterClass(bc, cache_dir, test_mode)
    adapter.create_anndata()
    print(f"Execution took {time.time() - start_time:.2f} seconds")

if merge:
    # List of files to load
    files = [
        "VDJdb",
        "IEDB",
        "McPAS-TCR",
        "TCR3d",
        "NeoTCR",
        "CEDAR",
        "TRAIT"
    ]
    adatas = {}

    for f in files:
        file_path = cache_dir / (f+"_anndata.h5ad")
        if file_path.exists():
            print(f"Loading {f}...")
            adatas[f] = ad.read_h5ad(file_path)
            print(f"  Shape: {adatas[f].shape}")
            print(f"  Columns: {adatas[f].obs.columns.tolist()}")
        else:
            print(f"Warning: {f} not found in {cache_dir}")

    # Quick summary table
    for f in adatas:
        print({
        "File": f,
        "Observations": adatas[f].n_obs,
        "Variables": adatas[f].n_vars,
        "Columns": len(adatas[f].obs.columns)
    } )

    common_cols = set.intersection(*(set(adatas[f].obs.columns) for f in adatas))
    print(f"Common columns: {common_cols}")

    # Concatenate all AnnData objects
    merged_adata = ad.concat(adatas, label="source", index_unique="_")
    print(f"Merged shape: {merged_adata.shape}")
    if deduplicate:
        # Deduplicate: Keep one entry per clonotype. 
        # Note: We might want to keep other info, but by default we take the first.
        with ir.get.airr_context(merged_adata, ["v_call", "junction_aa"], chain=["VJ_1", "VDJ_1"]) as m:
            deduplicated_adata = merged_adata[~m.obs.duplicated(subset=['VJ_1_junction_aa', 'VJ_1_v_call', 'VDJ_1_v_call', 'VDJ_1_junction_aa', 'iedb_iri']), :].copy()

        print(f"Number of entries after deduplication: {deduplicated_adata.n_obs}")
        # Convert PMID column to string to avoid serialization errors with mixed types
        deduplicated_adata.obs['PMID'] = deduplicated_adata.obs['PMID'].astype(str)
        deduplicated_adata.write_h5ad(cache_dir / "deduplicated_anndata.h5ad")
        print(f"Merged AnnData saved to {cache_dir / 'deduplicated_anndata.h5ad'}")
    else:
        # Convert PMID column to string to avoid serialization errors with mixed types
        merged_adata.obs['PMID'] = merged_adata.obs['PMID'].astype(str)
        merged_adata.write_h5ad(cache_dir / "merged_anndata.h5ad")
        print(f"Merged AnnData saved to {cache_dir / 'merged_anndata.h5ad'}")
