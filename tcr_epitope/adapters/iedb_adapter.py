import logging
import os
from pathlib import Path

import pandas as pd
from biocypher import BioCypher

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import get_pmids_batch, harmonize_sequences

logger = logging.getLogger(__name__)

import zipfile

import requests


class IEDBAdapter(BaseAdapter):
    """BioCypher adapter for the Immune Epitope Database (IEDB)[https://www.iedb.org/].

    Parameters
    ----------
    bc
        BioCypher instance for DB download.
    test
        If `True`, only a subset of the data will be loaded for testing purposes.
    prefer_calculated
        If `True`, calculated values are preferred over curated values. If `False`, curated values are preferred.
    """

    DB_URL = "https://www.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip"
    DB_DIR = "iedb_latest"
    TCR_FNAME = "tcr_full_v3.csv"
    BCR_FNAME = "bcr_full_v3.csv"

    def get_latest_release(self, bc: BioCypher) -> tuple[str, str]:
        # Create cache directory manually
        cache_dir = Path(bc._cache_directory) / "iedb_latest"
        cache_dir.mkdir(parents=True, exist_ok=True)

        zip_file_path = cache_dir / "receptor_full_v3.zip"
        extracted_dir = cache_dir / "receptor_full_v3_extracted"

        # Check if already downloaded and extracted
        if extracted_dir.exists():
            tcr_path = None
            bcr_path = None

            # Search for the specific files
            for root, _dirs, files in os.walk(extracted_dir):
                for file in files:
                    if file == self.TCR_FNAME:
                        tcr_path = os.path.join(root, file)
                    elif file == self.BCR_FNAME:
                        bcr_path = os.path.join(root, file)

            # If both files found, return them
            if tcr_path and os.path.exists(tcr_path) and bcr_path and os.path.exists(bcr_path):
                return tcr_path, bcr_path

        # Download with proper headers using requests
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.9",
            "Referer": "https://www.iedb.org/",
            "Accept-Encoding": "gzip, deflate, br",
            "Connection": "keep-alive",
        }

        try:
            print(f"Downloading IEDB data from {self.DB_URL}")
            response = requests.get(self.DB_URL, headers=headers, stream=True, timeout=60)
            response.raise_for_status()

            # Save the zip file
            with open(zip_file_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            print(f"Downloaded to {zip_file_path}")

            # Extract the zip file
            extracted_dir.mkdir(exist_ok=True)
            with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
                zip_ref.extractall(extracted_dir)

            print(f"Extracted to {extracted_dir}")

            # Search for the specific files (same logic as original)
            tcr_path = None
            bcr_path = None

            for root, _dirs, files in os.walk(extracted_dir):
                for file in files:
                    if file == self.TCR_FNAME:
                        tcr_path = os.path.join(root, file)
                    elif file == self.BCR_FNAME:
                        bcr_path = os.path.join(root, file)

            if not tcr_path or not os.path.exists(tcr_path):
                raise FileNotFoundError(f"TCR file '{self.TCR_FNAME}' not found in IEDB database")

            if not bcr_path or not os.path.exists(bcr_path):
                raise FileNotFoundError(f"BCR file '{self.BCR_FNAME}' not found in IEDB database")

            print(f"Found TCR file: {tcr_path}")
            print(f"Found BCR file: {bcr_path}")

            return tcr_path, bcr_path

        except requests.RequestException as e:
            raise FileNotFoundError(f"Failed to download IEDB database: {e}")
        except zipfile.BadZipFile as e:
            raise FileNotFoundError(f"Downloaded file is not a valid zip: {e}")
        except Exception as e:
            raise FileNotFoundError(f"Error processing IEDB download: {e}")

    # def get_latest_release(self, bc: BioCypher) -> str:
    #     # Download IEDB
    #     iedb_resource = FileDownload(
    #         name=self.DB_DIR,
    #         url_s=self.DB_URL,
    #         lifetime=30,
    #         is_dir=False,
    #         headers = headers,
    #     )

    #     iedb_paths = bc.download(iedb_resource)
    #     db_dir = Path(iedb_paths[0]).parent
    #     for root, _dirs, files in os.walk(db_dir):
    #         for file in files:
    #             if file == self.TCR_FNAME:
    #                 tcr_path = os.path.join(root, file)
    #             elif file == self.BCR_FNAME:
    #                 bcr_path = os.path.join(root, file)

    #     if not tcr_path or not os.path.exists(tcr_path):
    #         raise FileNotFoundError(f"Failed to download IEDB database from {self.DB_URL}")

    #     return tcr_path, bcr_path

    def read_table(
        self, bc: BioCypher, table_path: str, test: bool = False, prefer_calculated: bool = True
    ) -> pd.DataFrame:
        tcr_table_path, bcr_table_path = table_path

        tcr_table = pd.read_csv(tcr_table_path, header=[0, 1])
        tcr_table.columns = tcr_table.columns.map(" ".join)
        tcr_table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        tcr_table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY

        bcr_table = pd.read_csv(bcr_table_path, header=[0, 1])
        bcr_table.columns = bcr_table.columns.map(" ".join)
        bcr_table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.IGH_KEY
        bcr_table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.IGL_KEY

        table = pd.concat([tcr_table, bcr_table], ignore_index=True)
        if test:
            table = table.sample(frac=0.01, random_state=42)
        # Replace NaN and empty strings with None
        table = table.replace(["", "nan"], None).where(pd.notnull, None)

        # Fill columns based on preference: calculated vs curated
        if prefer_calculated:
            # Fill calculated columns with curated values if calculated is empty
            for chain_num in [1, 2]:
                table[f"Chain {chain_num} CDR3 Calculated"] = table[f"Chain {chain_num} CDR3 Calculated"].fillna(
                    table[f"Chain {chain_num} CDR3 Curated"]
                )
                table[f"Chain {chain_num} Calculated V Gene"] = table[f"Chain {chain_num} Calculated V Gene"].fillna(
                    table[f"Chain {chain_num} Curated V Gene"]
                )
                table[f"Chain {chain_num} Calculated J Gene"] = table[f"Chain {chain_num} Calculated J Gene"].fillna(
                    table[f"Chain {chain_num} Curated J Gene"]
                )
        else:
            # Fill curated columns with calculated values if curated is empty
            for chain_num in [1, 2]:
                table[f"Chain {chain_num} CDR3 Curated"] = table[f"Chain {chain_num} CDR3 Curated"].fillna(
                    table[f"Chain {chain_num} CDR3 Calculated"]
                )
                table[f"Chain {chain_num} Curated V Gene"] = table[f"Chain {chain_num} Curated V Gene"].fillna(
                    table[f"Chain {chain_num} Calculated V Gene"]
                )
                table[f"Chain {chain_num} Curated J Gene"] = table[f"Chain {chain_num} Curated J Gene"].fillna(
                    table[f"Chain {chain_num} Calculated J Gene"]
                )

        # Choose column names based on preference
        if prefer_calculated:
            cdr3_col_1 = "Chain 1 CDR3 Calculated"
            cdr3_col_2 = "Chain 2 CDR3 Calculated"
            v_gene_col_1 = "Chain 1 Calculated V Gene"
            v_gene_col_2 = "Chain 2 Calculated V Gene"
            j_gene_col_1 = "Chain 1 Calculated J Gene"
            j_gene_col_2 = "Chain 2 Calculated J Gene"
        else:
            cdr3_col_1 = "Chain 1 CDR3 Curated"
            cdr3_col_2 = "Chain 2 CDR3 Curated"
            v_gene_col_1 = "Chain 1 Curated V Gene"
            v_gene_col_2 = "Chain 2 Curated V Gene"
            j_gene_col_1 = "Chain 1 Curated J Gene"
            j_gene_col_2 = "Chain 2 Curated J Gene"

        rename_cols = {
            "Epitope Name": REGISTRY_KEYS.EPITOPE_KEY,
            "Epitope IEDB IRI": REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
            "Epitope Source Molecule": REGISTRY_KEYS.ANTIGEN_KEY,
            "Epitope Source Organism": REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            "Assay MHC Allele Names": REGISTRY_KEYS.MHC_GENE_1_KEY,
            cdr3_col_1: REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            cdr3_col_2: REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            v_gene_col_1: REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            j_gene_col_1: REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            v_gene_col_2: REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            j_gene_col_2: REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            "Chain 1 Organism IRI": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            "Chain 2 Organism IRI": REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY,
            REGISTRY_KEYS.CHAIN_1_TYPE_KEY: REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
            REGISTRY_KEYS.CHAIN_2_TYPE_KEY: REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
            "Reference IEDB IRI": REGISTRY_KEYS.PUBLICATION_KEY,
        }

        table = table.rename(columns=rename_cols)
        table = table[list(rename_cols.values())]

        # Extract iedb ID from the url
        table[REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY] = (
            "iedb:" + table[REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY].str.extract(r"/epitope/(\d+)$")[0]
        )

        # Preprocesses CDR3 sequences, epitope sequences, and gene names
        table_preprocessed = harmonize_sequences(bc, table)

        ref_urls = table_preprocessed[REGISTRY_KEYS.PUBLICATION_KEY].dropna().unique().tolist()
        ref_map = get_pmids_batch(bc, ref_urls)
        table_preprocessed[REGISTRY_KEYS.PUBLICATION_KEY] = table_preprocessed[REGISTRY_KEYS.PUBLICATION_KEY].map(
            ref_map
        )

        return table_preprocessed

    def get_nodes(self):
        # chain 1
        yield from self._generate_nodes_from_table(
            subset_cols=[
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            ],
            unique_cols=[
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            ],
            property_cols=[
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            ],
        )

        # chain 2
        yield from self._generate_nodes_from_table(
            subset_cols=[
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY,
            ],
            unique_cols=[
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            ],
            property_cols=[
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY,
            ],
        )

        # epitope
        yield from self._generate_nodes_from_table(
            subset_cols=[
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.PUBLICATION_KEY,
            ],
            unique_cols=[
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
            ],
            property_cols=[
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.PUBLICATION_KEY,
            ],
        )

    def get_edges(self):
        # chain 1 - chain 2
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            ],
            [
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            ],
            source_unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
        )

        # chain 1 - epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            ],
            REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
            source_unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
        )

        # chain 2 - epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            ],
            REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
            source_unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
        )
