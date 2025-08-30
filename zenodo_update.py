import os
import shutil
import tempfile
from argparse import ArgumentParser

import requests
from biocypher import BioCypher
from zenodo_client import update_zenodo

from tcr_epitope.adapters.mcpas_adapter import MCPASAdapter
from tcr_epitope.adapters.utils import save_airr_cells_json


def delete_existing_files(deposit_id, access_token, sandbox=True):
    """Delete all existing files from a Zenodo deposit"""

    base_url = "https://sandbox.zenodo.org/api" if sandbox else "https://zenodo.org/api"
    params = {"access_token": access_token}

    # Get the draft deposit to access files
    res = requests.get(f"{base_url}/deposit/depositions/{deposit_id}", params=params)
    res.raise_for_status()
    deposit = res.json()

    # Delete each file
    for file_info in deposit.get("files", []):
        file_id = file_info["id"]
        delete_url = f"{base_url}/deposit/depositions/{deposit_id}/files/{file_id}"

        print(f"Deleting file: {file_info['filename']}")
        delete_res = requests.delete(delete_url, params=params)
        delete_res.raise_for_status()
        print(f"Successfully deleted: {file_info['filename']}")


def create_new_version(deposit_id, access_token, sandbox=True):
    """Create a new version of a Zenodo deposit"""

    base_url = "https://sandbox.zenodo.org/api" if sandbox else "https://zenodo.org/api"
    params = {"access_token": access_token}

    # Create new version
    new_version_url = f"{base_url}/deposit/depositions/{deposit_id}/actions/newversion"
    print(new_version_url)
    res = requests.post(new_version_url, params=params)
    print(res["errors"])
    res.raise_for_status()

    # Extract new version ID from response
    new_version_id = res.json()["links"]["latest_draft"].split("/")[-1]
    print(f"New version ID: {new_version_id}")
    return new_version_id


def create_and_upload_kg(test_mode=False, cache_dir=None, zenodo_token=None):
    """Create knowledge graph and upload directly to Zenodo."""

    # Set Zenodo token
    if zenodo_token:
        os.environ["ZENODO_TOKEN_SANDBOX"] = zenodo_token
    elif "ZENODO_TOKEN_SANDBOX" not in os.environ:
        raise ValueError("No Zenodo token provided. Set ZENODO_TOKEN_SANDBOX env var or pass zenodo_token parameter")

    # Create knowledge graph
    print("Creating knowledge graph...")
    bc = BioCypher(cache_directory=cache_dir)

    adapters = [
        # VDJDBAdapter(bc, test_mode),
        MCPASAdapter(bc, test_mode),
        # TRAITAdapter(bc, test_mode),
        # IEDBAdapter(bc, test_mode),
        # TCR3DAdapter(bc, test_mode),
        # NeoTCRAdapter(bc, test_mode),
        # CEDARAdapter(bc, test_mode),
    ]

    for adapter in adapters:
        bc.add(adapter.get_nodes())
        bc.add(adapter.get_edges())

    airr_cells = bc.get_kg()
    print("Knowledge graph created.")

    # Upload to Zenodo
    print("Uploading to Zenodo...")
    temp_dir = tempfile.mkdtemp(prefix="tcr_epitope_kg_")

    try:
        # Write AIRR data to temporary file
        save_airr_cells_json(airr_cells, temp_dir)
        json_files = [f for f in os.listdir(temp_dir) if f.endswith(".json.gz")]

        if not json_files:
            raise FileNotFoundError("No JSON file was created by save_airr_cells_json")

        json_file_path = os.path.join(temp_dir, json_files[0])
        print(f"Created file: {json_file_path}")

        SANDBOX_DEP_ID = "308567"

        # last_id = get_last_id(SANDBOX_DEP_ID, zenodo_token, sandbox=True)
        # print("Deleting existing files from Zenodo deposit...")
        # delete_existing_files(last_id, zenodo_token, sandbox=True)
        # print("Existing files deleted.")

        # paths = [json_file_path,]
        # update_zenodo(last_id, paths, sandbox=True)

        last_id = get_last_id(SANDBOX_DEP_ID, zenodo_token, sandbox=True)
        print(last_id)

        # Step 2: Create a new version from the latest version
        print("Creating new version...")
        new_version_id = create_new_version(last_id, zenodo_token, sandbox=True)

        # Step 3: Delete files from the NEW version (not the old one)
        print("Deleting existing files from new version...")
        delete_existing_files(new_version_id, zenodo_token, sandbox=True)
        print("Existing files deleted.")

        # Step 4: Upload to the new version
        paths = [
            json_file_path,
        ]
        update_zenodo(new_version_id, paths, sandbox=True)

        print("Successfully uploaded to Zenodo!")

    finally:
        # Clean up temporary directory
        try:
            shutil.rmtree(temp_dir)
        except:
            pass


def get_last_id(deposit_id, access_token, sandbox=True):
    """Get the last version ID of a Zenodo deposit"""

    base_url = "https://sandbox.zenodo.org/api" if sandbox else "https://zenodo.org/api"
    params = {"access_token": access_token}

    # Get current version
    res = requests.get(f"{base_url}/records/{deposit_id}", params=params)
    res.raise_for_status()
    last_id = res.json()["id"]
    print(f"Last version ID: {last_id}")
    return last_id


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--test", action="store_true", help="Run in test mode")
    parser.add_argument("--cache-dir", help="Cache directory for BioCypher")
    parser.add_argument("--zenodo-token", help="Zenodo access token")

    args = parser.parse_args()

    create_and_upload_kg(test_mode=args.test, cache_dir=args.cache_dir, zenodo_token=args.zenodo_token)

# export ZENODO_TOKEN=your_token
# python3 kg_to_zenodo.py --test --cache-dir ./cache
