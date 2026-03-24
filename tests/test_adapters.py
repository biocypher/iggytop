import os

import pandas as pd
import pytest
from biocypher import BioCypher
from iggytop.adapters.cedar_adapter import CEDARAdapter
from iggytop.adapters.iedb_adapter import IEDBAdapter
from iggytop.adapters.itrap_adapter import ITRAPAdapter
from iggytop.adapters.mcpas_adapter import MCPASAdapter
from iggytop.adapters.neotcr_adapter import NEOTCRAdapter
from iggytop.adapters.tcr3d_adapter import TCR3DAdapter
from iggytop.adapters.trait_adapter import TRAITAdapter
from iggytop.adapters.vdjdb_adapter import VDJDBAdapter


@pytest.fixture
def bc():
    return BioCypher(
        biocypher_config_path="src/iggytop/config/biocypher_config.yaml",
        schema_config_path="src/iggytop/config/schema_config.yaml",
    )


def test_vdjdb_availability(bc):
    adapter = VDJDBAdapter(bc, test=True)
    assert adapter.metadata["version"] is not None
    assert adapter.metadata["source_url"].startswith("https://")
    assert os.path.exists(adapter._table_path)


def test_mcpas_availability(bc):
    adapter = MCPASAdapter(bc, test=True)
    assert adapter.metadata["version"] == "latest"
    assert adapter.metadata["source_url"] == MCPASAdapter.DB_URL
    assert os.path.exists(adapter._table_path)


def test_iedb_availability(bc):
    adapter = IEDBAdapter(bc, test=True)
    assert adapter.metadata["version"] == "v3"
    tcr_path, bcr_path = adapter._table_path
    assert os.path.exists(tcr_path)
    assert os.path.exists(bcr_path)


def test_cedar_availability(bc):
    adapter = CEDARAdapter(bc, test=True)
    assert adapter.metadata["version"] == "latest"
    assert adapter.metadata["source_url"] == CEDARAdapter.DB_URL
    tcr_path, bcr_path = adapter._table_path
    assert os.path.exists(tcr_path)
    assert os.path.exists(bcr_path)


def test_neotcr_availability(bc):
    adapter = NEOTCRAdapter(bc, test=True)
    assert adapter.metadata["version"] == "latest"
    assert adapter.metadata["source_url"] == NEOTCRAdapter.RAW_URL
    assert os.path.exists(adapter._table_path)


def test_tcr3d_availability(bc):
    adapter = TCR3DAdapter(bc, test=True)
    assert adapter.metadata["version"] == "latest"
    assert adapter.metadata["source_url"] == TCR3DAdapter.DB_URL
    assert os.path.exists(adapter._table_path)


def test_trait_availability(bc):
    adapter = TRAITAdapter(bc, test=True)
    assert adapter.metadata["version"] == "latest"
    assert adapter.metadata["source_url"] == TRAITAdapter.DB_URL
    assert os.path.exists(adapter._table_path)


def test_itrap_availability(bc):
    adapter = ITRAPAdapter(bc, test=True)
    assert adapter.metadata["version"] == "latest"
    assert adapter.metadata["source_url"] == ITRAPAdapter.DB_URL
    assert os.path.exists(adapter._table_path)


def test_adapter_table_format(bc):
    # Test a single adapter for table structure
    adapter = VDJDBAdapter(bc, test=True)
    table = adapter.table
    assert isinstance(table, pd.DataFrame)
    assert not table.empty
    # Check for mandatory columns (from REGISTRY_KEYS)
    from iggytop.adapters.constants import REGISTRY_KEYS

    assert REGISTRY_KEYS.CHAIN_2_CDR3_KEY in table.columns
    assert REGISTRY_KEYS.CHAIN_1_CDR3_KEY in table.columns
    assert REGISTRY_KEYS.EPITOPE_KEY in table.columns
