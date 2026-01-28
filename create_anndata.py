from biocypher import BioCypher
from pathlib import Path
import platformdirs
import importlib.resources
import pandas as pd
import anndata as ad
import scirpy as ir

from iggytop.adapters.cedar_adapter import CEDARAdapter
from iggytop.adapters.iedb_adapter import IEDBAdapter
from iggytop.adapters.mcpas_adapter import MCPASAdapter
from iggytop.adapters.neotcr_adapter import NeoTCRAdapter
from iggytop.adapters.tcr3d_adapter import TCR3DAdapter
from iggytop.adapters.trait_adapter import TRAITAdapter
from iggytop.adapters.vdjdb_adapter import VDJDBAdapter
from iggytop.adapters.utils import _set_up_config, _set_up_schema

import os
import time
import warnings

merge = True  # whether to merge all datasets into a single AnnData
deduplicate = True  # whether to deduplicate the merged AnnData, set unique id's below
remove_entries_without_epitope  = True  # whether to remove entries without epitope info, this mainly affects MCPAS data (~24k entries)
rebuild = True  # whether to rebuild the individual AnnData files from source
cache_dir = platformdirs.user_cache_dir("iggytop_anndata") #cachedir must be a string (biocypher requirement)

adapters_to_include = ["VDJdb", "IEDB", "McPAS",  "TCR3d", "NeoTCR", "CEDAR", "TRAIT"]

test_mode = False
output_format = "airr"  # doesnt matter here


def aggregate_unique_joined(series, separator='|'):
    """
    Helper function to aggregate unique values into a joined string.
    Warns if string 'nan' are found.
    """ 
    values = set()
    for v in series:
        if pd.isna(v) or str(v).lower() == 'nan':
            continue
            
        s_v = str(v).strip()  
        if s_v:
            values.update([s_v])
    if len(values) == 0:
        values.update(['nan'])
    return separator.join(sorted(values))


def deduplicate_and_aggregate(adata, subset_cols, agg_cols, separator='|'):
    """
    Deduplicates AnnData based on subset_cols and aggregates values in agg_cols.
    Uses scirpy airr_context to access TCR-specific columns if needed.
    """
    with ir.get.airr_context(adata, ["v_call", "junction_aa"], chain=["VJ_1", "VDJ_1"]) as m:
        obs_df = m.obs.copy()
        
        # Verify columns exist
        for col in subset_cols + agg_cols:
            if col not in obs_df.columns:
                raise KeyError(f"Column '{col}' not found in AnnData observations.")

        # Define aggregation map
        agg_map = {col: lambda s: aggregate_unique_joined(s, separator) for col in agg_cols}
        
        # Group and aggregate
        # observed=True is used to avoid 'Product space too large' errors
        aggs = obs_df.groupby(subset_cols, observed=True, dropna=False, sort=False).agg(agg_map).reset_index()
        
        # Identify first occurrences to keep as base records
        # Use the same placeholder logic for duplication identification
        is_first = ~obs_df.duplicated(subset=subset_cols)

        # Map aggregated info back to the deduplicated AnnData
        first_entries_keys = obs_df.loc[is_first, subset_cols].reset_index()
        final_info = first_entries_keys.merge(aggs, on=subset_cols, how='left')

    deduplicated_adata = adata[is_first, :].copy()
    for col in agg_cols:
        deduplicated_adata.obs[col] = final_info[col].values
        
    return deduplicated_adata

config_path = _set_up_config(output_format, cache_dir)
schema_config_path = _set_up_schema(cache_dir)

bc = BioCypher(biocypher_config_path=config_path,
                schema_config_path=schema_config_path,
                cache_directory=cache_dir)


adapter_classes = {
    "VDJdb": VDJDBAdapter,
    "McPAS": MCPASAdapter,
    "TRAIT": TRAITAdapter,
    "IEDB": IEDBAdapter,
    "TCR3d": TCR3DAdapter,
    "NeoTCR": NeoTCRAdapter,
    "CEDAR": CEDARAdapter,
}


selected_adapters = [
    adapter_classes[name]
    for name in adapters_to_include
    if name in adapter_classes
]
if rebuild:
    for AdapterClass in selected_adapters:
        start_time = time.time()
        adapter = AdapterClass(bc, cache_dir, test_mode)
        adapter.create_anndata()
        print(f"Execution took {time.time() - start_time:.2f} seconds")

cache_dir = Path(cache_dir)
if merge:
    # List of files to load
    adatas = {}

    for f in adapters_to_include:
        file_path = cache_dir / (f+"_anndata.h5ad")
        if os.path.exists(file_path):
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
    print(f"Number of entries: {merged_adata.n_obs}")
    # Convert PMID column to string to avoid serialization errors with mixed types
    merged_adata.obs['PMID'] = merged_adata.obs['PMID'].astype(str)

    if remove_entries_without_epitope:
        initial_count = merged_adata.n_obs
        has_epitope = merged_adata.obs['epitope_sequence'] != 'nan'
        merged_adata = merged_adata[has_epitope, :].copy()
        removed_count = initial_count - merged_adata.n_obs
        print(f"Removed {removed_count} entries without epitope information. Remaining entries: {merged_adata.n_obs}")

    merged_adata.write_h5ad(cache_dir / "merged_anndata.h5ad")
    print(f"Merged shape: {merged_adata.shape}")
    if deduplicate:
        # Deduplicate and aggregate specific attributes
        subset_cols = ['VJ_1_junction_aa', 'VJ_1_v_call', 'VDJ_1_v_call', 'VDJ_1_junction_aa', 'iedb_iri']
        agg_cols = ['PMID', 'source']
        
        try:
            deduplicated_adata = deduplicate_and_aggregate(merged_adata, subset_cols, agg_cols)
        except (ValueError, KeyError) as e:
            print(f"Deduplication failed due to unexpected data state: {e}")
            raise

        print(f"Number of entries after deduplication: {deduplicated_adata.n_obs}")
        deduplicated_adata.write_h5ad(cache_dir / "deduplicated_anndata.h5ad")
        print(f"Merged AnnData saved to {cache_dir / 'deduplicated_anndata.h5ad'}")