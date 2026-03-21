import argparse
import gzip
import json
import os
from importlib import resources


def generate_release_assets(release_dir):
    """
    Processes metadata from cache and generates RELEASE_NOTES.md and metadata.json.
    This information is the same for the merged and deduplicated datasets.
    """

    metadata = {}
    for f in os.listdir(release_dir):
        if f.startswith("merged_airr_cells") and f.endswith(".json.gz"):
            file_path = os.path.join(release_dir, f)
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

    template_text = resources.files("iggytop.io").joinpath("RELEASE_NOTES_TEMPLATE.md").read_text()
    release_notes_content = template_text.replace("{{SOURCE_TABLE}}", table)

    release_notes_path = os.path.join(release_dir, "RELEASE_NOTES.md")
    with open(release_notes_path, "w") as f:
        f.write(release_notes_content)

    metadata_path = os.path.join(release_dir, "metadata.json")
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
