import argparse
import gzip
import json
import os


def generate_release_assets(cache_dir, release_dir):
    """
    Processes metadata from cache and generates RELEASE_NOTES.md and metadata.json.
    """
    # Find metadata from cache
    metadata = {}
    for f in os.listdir(cache_dir):
        if f.startswith("deduplicated_airr_cells") and f.endswith(".json.gz"):
            file_path = os.path.join(cache_dir, f)
            with gzip.open(file_path, "rt") as gz:
                data = json.load(gz)
                metadata = data.get("metadata", {})
                break

    if not metadata:
        print("No metadata found.")
        return False

    # Ensure release directory exists
    os.makedirs(release_dir, exist_ok=True)

    # Write release notes and metadata directly to release folder
    table = "| Data Source | Version | Changed | Checksum (SHA256) |\n"
    table += "| --- | --- | --- | --- |\n"
    for name, info in metadata.get("sources", {}).items():
        changed_str = "⚠️ YES" if info.get("has_changed") else "No"
        table += f"| {name} | {info.get('version', 'N/A')} | {changed_str} | `{info.get('checksum', 'N/A')[:10]}...` |\n"

    release_notes_path = os.path.join(release_dir, "RELEASE_NOTES.md")
    with open(release_notes_path, "w") as f:
        f.write("### IggyTop Data Release\n\n")
        f.write(
            "IggyTop (**I**mmunolo**g**ical **G**raph **Y**ielding **Top** receptor-epitope pairings) "
            "harmonizes and integrates various immunoreceptor-epitope databases using the BioCypher framework. "
            "This release includes updated data from multiple sources in AnnData and AIRR formats.\n\n"
        )
        f.write("For more information, please visit the [IggyTop GitHub repository](https://github.com/biocypher/iggytop).\n\n")
        f.write("#### Data Source Information\n\n")
        f.write(table)

    metadata_path = os.path.join(release_dir, "metadata.json")
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=4)

    print(f"Release assets generated in {release_dir}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Generate IggyTop release assets.")
    parser.add_argument("--cache-dir", default="/tmp/iggytop_cache", help="Directory containing cached data")
    parser.add_argument("--release-dir", default="/tmp/iggytop_release", help="Directory to write release assets")

    args = parser.parse_args()

    success = generate_release_assets(args.cache_dir, args.release_dir)
    if not success:
        exit(1)


if __name__ == "__main__":
    main()
