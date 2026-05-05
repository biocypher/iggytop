#!/usr/bin/env python3
"""
Test script to execute Jupyter notebooks and check for execution errors.
This test ensures notebooks run without failures without saving changes.
"""

import subprocess
import pytest
from pathlib import Path

# List of notebooks to test
NOTEBOOKS = [
    "docs/source/notebooks/tutorial_3k_tcr.ipynb",
    "docs/source/notebooks/database_summary.ipynb",
]

def test_notebooks_execution():
    """Test that specified notebooks execute without errors."""
    for notebook in NOTEBOOKS:
        assert Path(notebook).exists(), f"Notebook {notebook} does not exist"
        
        print(f"Executing {notebook}...")
        # Run nbconvert to execute the notebook
        # --execute: execute the notebook
        # --to notebook: output as notebook (but we discard)
        # --stdout: output to stdout (we capture and discard)
        result = subprocess.run(
            ["jupyter", "nbconvert", "--execute", notebook, "--to", "notebook", "--stdout"],
            capture_output=True,
            text=True,
            timeout=600,  # 10 minutes timeout
        )
        
        assert result.returncode == 0, f"Failed to execute {notebook}. STDERR: {result.stderr}"
        print(f"✓ Successfully executed {notebook}")