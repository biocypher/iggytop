import os
import pytest
import pandas as pd
from biocypher import BioCypher
from iggytop.adapters.vdjdb_adapter import VDJDBAdapter
from iggytop.adapters.mcpas_adapter import MCPASAdapter
from iggytop.adapters.iedb_adapter import IEDBAdapter

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

def test_adapter_table_format(bc):
    # Test a single adapter for table structure
    adapter = VDJDBAdapter(bc, test=True)
    table = adapter.table
    assert isinstance(table, pd.DataFrame)
    assert not table.empty
    # Check for mandatory columns (from REGISTRY_KEYS)
    from iggytop.adapters.constants import REGISTRY_KEYS
    assert REGISTRY_KEYS.CHAIN_1_CDR3_KEY in table.columns
    assert REGISTRY_KEYS.EPITOPE_KEY in table.columns
