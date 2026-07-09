import argparse
import os
from datetime import datetime
from pathlib import Path

import anndata as ad
import platformdirs
import scirpy as ir
from scirpy.pp import index_chains

from iggytop.adapters.utils import (
    _flush_tt_warnings,
    deduplicate_and_aggregate,
    get_previous_release_metadata,
    save_airr_cells_json,
)
from iggytop.io.create_knowledge_graph import DEFAULT_ADAPTERS, build_adapters, write_knowledge_graph


def _save_adata(adata: ad.AnnData, path: Path, *, name: str, metadata: dict):
    index_chains(adata)
    adata.uns["DB"] = {"name": name, "date_created": datetime.now().isoformat(), "version": metadata["iggytop_version"]}
    adata.uns["iggytop_metadata"] = metadata
    # Opt in to writing pd.arrays.StringArray in obs/var.
    ad.settings.allow_write_nullable_strings = True
    adata.write_h5ad(path, compression="gzip")
    print(f"{path.name} saved to {path}")


def _filter_10x(adata: ad.AnnData) -> ad.AnnData:
    """Filter out the 10X Genomics dataset, which has been criticized for poor confidence."""
    pmid = adata.obs["PMID"].astype("string")
    is_10x = (pmid == "no_pmid_1036521").fillna(False) | pmid.str.contains("https://www.10xgenomics.com", na=False)
    return adata[(~is_10x).to_numpy(dtype=bool)].copy()


def main():
    parser = argparse.ArgumentParser(
        description="Build the full iggytop release: merged AnnData/AIRR, deduplicated AnnData/AIRR, and the neo4j knowledge graph, "
        "from a single build of each adapter's source table."
    )
    parser.add_argument("--test-mode", action="store_true", default=False, help="Run in test mode with a small subset of data.")
    parser.add_argument(
        "--cache-dir",
        type=str,
        default=platformdirs.user_cache_dir("iggytop_airr"),
        help="Directory for caching results.",
    )
    parser.add_argument(
        "--adapters",
        nargs="+",
        default=DEFAULT_ADAPTERS,
        help="List of adapters to include (e.g., --adapters VDJDB CEDAR). Defaults to all.",
    )
    parser.add_argument(
        "--not_merged",
        action="store_true",
        default=False,
        help="Skip saving the merged (non-deduplicated) AnnData/AIRR output.",
    )
    parser.add_argument(
        "--not_deduplicate",
        action="store_true",
        default=False,
        help="Skip generating the deduplicated AnnData/AIRR output.",
    )
    parser.add_argument(
        "--not_graph",
        action="store_true",
        default=False,
        help="Skip generating the neo4j knowledge graph output.",
    )
    parser.add_argument(
        "--graph-output-format",
        type=str,
        default="neo4j",
        choices=["airr", "neo4j", "networkx", "docker"],
        help="BioCypher output format for the graph output (default: neo4j).",
    )
    parser.add_argument(
        "--adata_only",
        action="store_true",
        default=False,
        help="Whether not to save the AnnData as AIRR JSON.",
    )
    parser.add_argument(
        "--filter-10x",
        action="store_true",
        dest="filter_10x",
        default=False,
        help="Filter out 10X Genomics records from the deduplicated dataset only (default: False).",
    )
    parser.add_argument(
        "--tag",
        type=str,
        default="latest",
        help="Release tag to embed in metadata (e.g. data-YYYY.MM.DD.HHMMSS).",
    )

    args = parser.parse_args()

    adapters_to_include = args.adapters
    merge = len(adapters_to_include) > 1
    save_merged = not args.not_merged
    deduplicate = not args.not_deduplicate
    include_graph = not args.not_graph
    save_single_adapter_data = False
    save_airr_json = not args.adata_only
    receptors_to_include = ["TCR", "BCR"]
    cache_dir = args.cache_dir

    _bc, adapters = build_adapters(
        cache_dir=cache_dir,
        test_mode=args.test_mode,
        receptors_to_include=receptors_to_include,
        adapters_to_include=adapters_to_include,
    )

    # Fetch previous release metadata for change detection
    prev_metadata = get_previous_release_metadata() or {}
    prev_sources = prev_metadata.get("sources", {})

    global_metadata = {
        "iggytop_version": args.tag,
        "release_date": datetime.now().isoformat(),
        "sources": {},
    }

    for adapter in adapters:
        # Update adapter metadata with change information
        prev_source = prev_sources.get(adapter.db_name, {})
        prev_source_version = prev_source.get("version")
        prev_source_checksum = prev_source.get("checksum")
        adapter.set_metadata(previous_version=prev_source_version, previous_checksum=prev_source_checksum)
        global_metadata["sources"][adapter.db_name] = adapter.metadata

        adapter.create_anndata()
        if save_airr_json and save_single_adapter_data:
            save_airr_cells_json(
                adapter.airr_cells,
                directory=cache_dir,
                filename=f"{adapter.db_name}_airr_cells",
                metadata=adapter.metadata,
            )

    if include_graph:
        # Reuses the already-built adapters (and their memoized tables), so the source data isn't re-read.
        write_knowledge_graph(adapters, cache_dir, output_format=args.graph_output_format)

    cache_dir = Path(cache_dir)  # Biocypher likes paths as strings

    if merge and (save_merged or deduplicate):
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
        for f, adata in adatas.items():
            print(
                {
                    "File": f,
                    "Observations": adata.n_obs,
                    "Variables": adata.n_vars,
                    "Columns": len(adata.obs.columns),
                }
            )

        common_cols = set.intersection(*(set(adatas[f].obs.columns) for f in adatas))
        print(f"Common columns: {common_cols}")

        # Concatenate all AnnData objects
        merged_adata = ad.concat(adatas, label="source", index_unique="_")
        print(f"Number of entries: {merged_adata.n_obs}")

        # Convert object columns to string to avoid serialization issues with h5py (e.g. for PMID)
        for col in merged_adata.obs.columns:
            if merged_adata.obs[col].dtype == object:
                merged_adata.obs[col] = merged_adata.obs[col].astype("string")

        if deduplicate:
            # 10X filtering only ever applies to the deduplicated output, not the merged/raw one.
            dedup_source_adata = _filter_10x(merged_adata) if args.filter_10x else merged_adata

            # Drop records where junction is missing for both chains — they can't be meaningfully deduplicated
            with ir.get.airr_context(dedup_source_adata, ["junction_aa"], chain=["VJ_1", "VDJ_1"]) as m:
                has_junction = m.obs["VJ_1_junction_aa"].notna() | m.obs["VDJ_1_junction_aa"].notna()
            merged_adata_for_dedup = dedup_source_adata[has_junction].copy()

            # Deduplicate and aggregate specific attributes
            subset_cols = [
                "VJ_1_junction_aa",
                "VJ_1_v_call",
                "VDJ_1_v_call",
                "VJ_1_j_call",
                "VDJ_1_j_call",
                "VDJ_1_junction_aa",
                "epitope_sequence",
            ]  # epitope IRI can be ambiguous
            agg_cols = ["PMID", "source"]

            try:
                deduplicated_adata = deduplicate_and_aggregate(merged_adata_for_dedup, subset_cols, agg_cols)
            except (ValueError, KeyError) as e:
                print(f"Deduplication failed due to unexpected data state: {e}")
                raise

            print(f"Number of entries after deduplication: {deduplicated_adata.n_obs}")

        # Save result to AnnData
        if save_merged:
            _save_adata(merged_adata, cache_dir / "merged_anndata.h5ad", name="iggytop_merged", metadata=global_metadata)
        if deduplicate:
            _save_adata(deduplicated_adata, cache_dir / "deduplicated_anndata.h5ad", name="iggytop_deduplicated", metadata=global_metadata)

        # Optional: Export to AIRR JSON format
        if save_airr_json:
            if save_merged:
                merged_airr_list = ir.io.to_airr_cells(merged_adata)
                save_airr_cells_json(merged_airr_list, directory=cache_dir, filename="merged_airr_cells", metadata=global_metadata)
            if deduplicate:
                deduplicated_airr_list = ir.io.to_airr_cells(deduplicated_adata)
                save_airr_cells_json(
                    deduplicated_airr_list,
                    directory=cache_dir,
                    filename="deduplicated_airr_cells",
                    metadata=global_metadata,
                )

    _flush_tt_warnings()


if __name__ == "__main__":
    main()
