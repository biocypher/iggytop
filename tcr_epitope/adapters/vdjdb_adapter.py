import os
from pathlib import Path

import numpy as np
import pandas as pd
from biocypher import BioCypher, FileDownload
from github import Github

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import get_iedb_ids_batch, harmonize_sequences


class VDJDBAdapter(BaseAdapter):
    """BioCypher adapter for the VDJdb database (https://vdjdb.cdr3.net/).

    Parameters
    ----------
    bc
        BioCypher instance for DB download.
    test
        If `True`, only a subset of the data will be loaded for testing purposes.
    """

    REPO_NAME = "antigenomics/vdjdb-db"
    DB_DIR = "vdjdb_latest"
    DB_FNAME = "vdjdb.txt"

    def get_latest_release(self, bc: BioCypher) -> str:
        repo = Github().get_repo(self.REPO_NAME)
        db_url = repo.get_latest_release().get_assets()[0].browser_download_url
        # db_url = "https://github.com/antigenomics/vdjdb-db/releases/download/pyvdjdb-2025-02-21/vdjdb-2025-02-21.zip"

        vdjdb_resource = FileDownload(
            name=self.DB_DIR,
            url_s=db_url,
            lifetime=30,
            is_dir=False,
        )

        vdjdb_paths = bc.download(vdjdb_resource)

        db_dir = Path(vdjdb_paths[0]).parent
        for root, _dirs, files in os.walk(db_dir):
            for file in files:
                if file == self.DB_FNAME:
                    db_path = os.path.join(root, file)

        if not db_path or not os.path.exists(db_path):
            raise FileNotFoundError(f"Failed to download VDJdb database from {db_url}")

        return db_path

    def read_table(self, bc: BioCypher, table_path: str, test: bool = False) -> pd.DataFrame:
        table = pd.read_csv(table_path, sep="\t")
        if test:
            table = table.sample(frac=0.01, random_state=42)
        # Replace NaN and empty strings with None
        table = table.replace(["", "nan"], None).where(pd.notnull, None)

        # WITH THIS OPTIMIZED METHOD:
        table = self._transform_paired_data_efficient(table)

        # Rest of your code stays the same:
        rename_cols = {
            "cdr3_chain_1": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,  # Note: changed from cdr3_chain_1
            "v.segm_chain_1": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,  # Note: changed from v.segm_chain_1
            "j.segm_chain_1": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,  # Note: changed from j.segm_chain_1

            "cdr3_chain_2": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,  # Note: changed from cdr3_chain_2
            "v.segm_chain_2": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,  # Note: changed from v.segm_chain_2
            "j.segm_chain_2": REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,  # Note: changed from j.segm_chain_2
            "species": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,

            "antigen.epitope": REGISTRY_KEYS.EPITOPE_KEY,
            "antigen.gene": REGISTRY_KEYS.ANTIGEN_KEY,
            "antigen.species": REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            "reference.id": REGISTRY_KEYS.PUBLICATION_KEY,

            "mhc.class": REGISTRY_KEYS.MHC_CLASS_KEY,
            "mhc.a": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "mhc.b": REGISTRY_KEYS.MHC_GENE_2_KEY,
        }
        
        # Only rename columns that exist in the table
        table = table.rename(columns=rename_cols)
        table = table[list(rename_cols.values())]

        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY]
        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY

        # Map epitope sequences to IEDB IDs
        valid_epitopes = table[REGISTRY_KEYS.EPITOPE_KEY].dropna().drop_duplicates().tolist()
        if len(valid_epitopes) > 0:
            epitope_map = get_iedb_ids_batch(bc, valid_epitopes)

        # Apply the mapping to create the IEDB ID column
        table[REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY] = table[REGISTRY_KEYS.EPITOPE_KEY].map(epitope_map)

        # Preprocesses CDR3 sequences, epitope sequences, and gene names
        table_preprocessed = harmonize_sequences(table)

        return table_preprocessed

    def _transform_paired_data_efficient(self, df):
        """Efficient transformation that handles ALL cases correctly."""
        
        # 1. Separate unpaired (complex.id == 0)
        unpaired = df[df['complex.id'] == 0].copy()
        
        # 2. Find ACTUALLY paired data (duplicated non-zero complex.ids)
        paired_mask = (df['complex.id'] != 0) & (df['complex.id'].duplicated(keep=False))
        paired = df[paired_mask].copy()
        
        # 3. Find incomplete pairs (non-zero, non-duplicated complex.ids)
        incomplete = df[(df['complex.id'] != 0) & (~paired_mask)].copy()
        
        result_parts = []
        
        # Process complete pairs
        if len(paired) > 0:
            # Check which complex.ids have both TRA and TRB
            chain_counts = paired.groupby('complex.id')['gene'].apply(set)
            complete_mask = chain_counts.apply(lambda x: {'TRA', 'TRB'}.issubset(x))
            complete_complexes = chain_counts[complete_mask].index
            
            if len(complete_complexes) > 0:
                complete_data = paired[paired['complex.id'].isin(complete_complexes)]
                tra_complete = complete_data[complete_data['gene'] == 'TRA']
                trb_complete = complete_data[complete_data['gene'] == 'TRB']
                
                merge_cols = ['complex.id', 'antigen.epitope', 'antigen.gene', 'antigen.species', 
                            'mhc.class', 'mhc.a', 'mhc.b', 'reference.id']
                
                paired_result = tra_complete.merge(
                    trb_complete[merge_cols + ['cdr3', 'v.segm', 'j.segm']], 
                    on=merge_cols, suffixes=('_chain_1', '_chain_2')
                )
                paired_result = paired_result.rename(columns={
                    'cdr3': 'cdr3_tra', 'v.segm': 'v.segm_tra', 'j.segm': 'j.segm_tra'
                })
                result_parts.append(paired_result)

        
        # Process all single chains (unpaired + incomplete pairs)
        all_singles = pd.concat([unpaired, incomplete], ignore_index=True) if len(incomplete) > 0 else unpaired
        
        if len(all_singles) > 0:
            single_tra = self._process_single_chain(all_singles[all_singles['gene'] == 'TRA'], 'tra')
            single_trb = self._process_single_chain(all_singles[all_singles['gene'] == 'TRB'], 'trb')
            result_parts.extend([single_tra, single_trb])
        
        # Combine all
        result_parts = [df for df in result_parts if len(df) > 0]
        return pd.concat(result_parts, ignore_index=True) if result_parts else df

    def _process_single_chain(self, df, chain_type):
        """Process single chain data (TRA or TRB only)."""
        if len(df) == 0:
            return df
            
        result = df.copy()
        if chain_type == 'tra':
            result['cdr3_chain_1'] = result['cdr3']
            result['v.segm_chain_1'] = result['v.segm'] 
            result['j.segm_chain_1'] = result['j.segm']
            result['cdr3_chain_2'] = None
            result['v.segm_chain_2'] = None
            result['j.segm_chain_2'] = None
        else:  # trb
            result['cdr3_chain_2'] = result['cdr3']
            result['v.segm_chain_2'] = result['v.segm']
            result['j.segm_chain_2'] = result['j.segm'] 
            result['cdr3_chain_1'] = None
            result['v.segm_chain_1'] = None
            result['j.segm_chain_1'] = None
            
        return result.drop(columns=['cdr3', 'v.segm', 'j.segm'])

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
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.MHC_GENE_2_KEY,
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
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.MHC_GENE_2_KEY,
                REGISTRY_KEYS.PUBLICATION_KEY,
            ],
        )

    def get_edges(self):
        # chain 1 to chain 2
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

        # chain 1 to epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            ],
            [REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY],
            source_unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
        )

        # chain 2 to epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            ],
            [REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY],
            source_unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
        )
        