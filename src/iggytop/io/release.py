import argparse
import gzip
import json
from datetime import datetime
from importlib import resources
from pathlib import Path

from iggytop.adapters.utils import get_file_checksum


def generate_release_assets(release_dir):
    """
    Processes metadata from cache and generates RELEASE_NOTES.md and metadata.json.
    This information is the same for the merged and deduplicated datasets.
    """

    release_dir = Path(release_dir)
    metadata = {}
    assets = {}
    for f in release_dir.iterdir():
        # get checksums of asset files that will be attached to the release
        if f.is_file():
            assets[f.name] = get_file_checksum(f)
        # get source file metadata from merged_airr_cells json
        if f.name.startswith("merged_airr_cells") and f.name.endswith(".json.gz"):
            with gzip.open(f, "rt") as gz:
                data = json.load(gz)
                metadata = data.get("metadata", {})
                break

    if not metadata:
        print("No metadata found.")
        return False

    metadata["assets"] = assets

    # Ensure release directory exists
    release_dir.mkdir(parents=True, exist_ok=True)

    # Write release notes and metadata directly to release folder
    table = "| Data Source | Version | Date downloaded | Changed | Checksum (SHA256) |\n"
    table += "| --- | --- | --- | --- | --- |\n"
    for name, info in metadata.get("sources", {}).items():
        changed_str = "⚠️ YES" if info.get("has_changed") else "No"
        try:
            download_date = datetime.fromisoformat(info["download_date"]).strftime("%Y-%m-%d")
        except (ValueError, TypeError, KeyError):
            download_date = "N/A"

        # Check for both 'checksums' and 'checksum' keys
        checksums = info.get("checksums") or info.get("checksum") or "N/A"

        if isinstance(checksums, str):
            checksums_list = [checksums]
        elif isinstance(checksums, list):
            checksums_list = checksums
        else:
            checksums_list = ["N/A"]

        # for the table, let's show the first 8 chars only. Full information is in metadata
        checksums_str = "|".join([str(c)[:8] for c in checksums_list if c])
        if not checksums_str:
            checksums_str = "N/A"

        table += f"| {name} | {info.get('version', 'N/A')} | {download_date} | {changed_str} | {checksums_str} |\n"

    template_text = resources.files("iggytop.io").joinpath("RELEASE_NOTES_TEMPLATE.md").read_text()
    release_notes_content = template_text.replace("{{SOURCE_TABLE}}", table)

    release_notes_path = release_dir / "RELEASE_NOTES.md"
    with open(release_notes_path, "w") as f:
        f.write(release_notes_content)

    # asset information in `metadata.json` will be used downstream by scirpy to fetch the latest release.
    # Be careful when modifying this
    metadata_path = release_dir / "metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=4)

    print(f"Release assets generated in {release_dir}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Generate IggyTop release assets.")
    parser.add_argument("--release-dir", default="/tmp/iggytop_release", help="Directory to write release assets")

    args = parser.parse_args()

    success = generate_release_assets(args.release_dir)
    if not success:
        exit(1)


if __name__ == "__main__":
    main()
