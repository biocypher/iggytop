from biocypher import BioCypher
from pathlib import Path
import os
from datetime import datetime
import platformdirs
import pandas as pd
import anndata as ad
import scirpy as ir
from scirpy.pp import index_chains

from iggytop.adapters.cedar_adapter import CEDARAdapter
from iggytop.adapters.iedb_adapter import IEDBAdapter
from iggytop.adapters.itrap_adapter import ITRAPAdapter
from iggytop.adapters.mcpas_adapter import MCPASAdapter
from iggytop.adapters.neotcr_adapter import NEOTCRAdapter
from iggytop.adapters.tcr3d_adapter import TCR3DAdapter
from iggytop.adapters.trait_adapter import TRAITAdapter

from iggytop.adapters.vdjdb_adapter import VDJDBAdapter
from iggytop.adapters.utils import (
    _set_up_config,
    _set_up_schema,
    get_previous_release_metadata,
    save_airr_cells_json,
    deduplicate_and_aggregate,
)


merge = True  # whether to merge all datasets into a single AnnData
deduplicate = True  # whether to deduplicate the merged AnnData, set unique id's below
save_single_adapter_data = False  # whether to save individual AnnData files for each adapter (in addition to the merged version)
save_airr_json = True  # whether to save the merged and deduplicated AnnData as AIRR JSON
receptors_to_include = ["TCR", "BCR"]
cache_dir = platformdirs.user_cache_dir("iggytop_airr") #cachedir must be a string (biocypher requirement)
os.makedirs(cache_dir, exist_ok=True)

adapters_to_include = ["ITRAP","VDJDB", "MCPAS", "IEDB", "TCR3D", "NEOTCR", "CEDAR", "TRAIT"]

test_mode = False
output_format = "airr"  # doesnt matter here


config_path = _set_up_config(output_format, cache_dir)
schema_config_path = _set_up_schema(cache_dir)

bc = BioCypher(biocypher_config_path=config_path,
                schema_config_path=schema_config_path,
                cache_directory=cache_dir)


adapter_classes = {
    "VDJDB": VDJDBAdapter,
    "MCPAS": MCPASAdapter,
    "TRAIT": TRAITAdapter,
    "ITRAP": ITRAPAdapter,
    "IEDB": IEDBAdapter,
    "TCR3D": TCR3DAdapter,
    "NEOTCR": NEOTCRAdapter,
    "CEDAR": CEDARAdapter,
}


selected_adapters = [
    adapter_classes[name]
    for name in adapters_to_include
    if name in adapter_classes
]
selected_adapters = [a for a in selected_adapters if any(receptor in receptors_to_include for receptor in a.available_receptors)]

# Fetch previous release metadata for change detection
prev_metadata = get_previous_release_metadata() or {}
prev_sources = prev_metadata.get("sources", {})

global_metadata = {
    "iggytop_version": "latest", # Could be fetched from bumpversion or package
    "release_date": datetime.now().isoformat(),
    "sources": {}
}

adapters = []
for AdapterClass in selected_adapters:
    adapter = AdapterClass(bc, cache_dir, receptors_to_include, test_mode)
    adapters.append(adapter)
    
    # Update adapter metadata with change information
    prev_source_version = prev_sources.get(adapter.db_name, {}).get("version")
    adapter.set_metadata(previous_version=prev_source_version)
    global_metadata["sources"][adapter.db_name] = adapter.metadata

    adapter.create_anndata()
    if save_airr_json and save_single_adapter_data:
        save_airr_cells_json(adapter.airr_cells, directory=cache_dir, filename=f"{adapter.db_name}_airr_cells", metadata=adapter.metadata)

cache_dir = Path(cache_dir)

if merge:
    adatas = {}

    for adapter in adapters:
        db_name = adapter.db_name
        file_path = cache_dir / (db_name + "_anndata.h5ad")
        if os.path.exists(file_path):
            print(f"Loading {db_name}...")
            adatas[db_name] = ad.read_h5ad(file_path)
            print(f"  Shape: {adatas[db_name].shape}")
            print(f"  Columns: {adatas[db_name].obs.columns.tolist()}")
        else:
            print(f"Warning: {db_name} not found in {cache_dir}")

    # Ensure all obs columns exist across datasets
    all_cols = set().union(*(set(a.obs.columns) for a in adatas.values()))
    for a in adatas.values():
        missing = all_cols - set(a.obs.columns)
        for col in missing:
            a.obs[col] = None

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
    print(f"Number of entries: {merged_adata.n_obs}")

    # Convert object columns to string to avoid serialization issues with h5py (e.g. for PMID)
    for col in merged_adata.obs.columns:
        if merged_adata.obs[col].dtype == object:
            merged_adata.obs[col] = merged_adata.obs[col].astype(str)

    merged_adata.write_h5ad(cache_dir / "merged_anndata.h5ad")
    print(f"Merged AnnData saved to {cache_dir / 'merged_anndata.h5ad'}")

    
    if deduplicate:
        # Deduplicate and aggregate specific attributes
        subset_cols = ['VJ_1_junction_aa', 'VJ_1_v_call', 'VDJ_1_v_call', 'VDJ_1_junction_aa', 'epitope_sequence'] #epitope IRI can be ambiguous
        agg_cols = ['PMID', 'source']
        
        try:
            deduplicated_adata = deduplicate_and_aggregate(merged_adata, subset_cols, agg_cols)
        except (ValueError, KeyError) as e:
            print(f"Deduplication failed due to unexpected data state: {e}")
            raise

        print(f"Number of entries after deduplication: {deduplicated_adata.n_obs}")
        index_chains(deduplicated_adata)
        
        # Store metadata in AnnData
        deduplicated_adata.uns["iggytop_metadata"] = global_metadata
        
        deduplicated_adata.write_h5ad(cache_dir / "deduplicated_anndata.h5ad")
        print(f"Deduplicated AnnData saved to {cache_dir / 'deduplicated_anndata.h5ad'}")

    # Optional: Export to AIRR JSON format
    if save_airr_json:
        merged_airr_list = ir.io.to_airr_cells(merged_adata)
        save_airr_cells_json(merged_airr_list, directory=cache_dir, filename="merged_airr_cells", metadata=global_metadata)
        if deduplicate:
            deduplicated_airr_list = ir.io.to_airr_cells(deduplicated_adata)
            save_airr_cells_json(deduplicated_airr_list, directory=cache_dir, filename="deduplicated_airr_cells", metadata=global_metadata)


