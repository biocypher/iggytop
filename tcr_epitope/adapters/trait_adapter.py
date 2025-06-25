import os

import pandas as pd
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import harmonize_sequences


class TRAITAdapter(BaseAdapter):
    """A Comprehensive Database for T-cell Receptor-Antigen Interactions (TRAIT)[https://pgx.zju.edu.cn/traitdb/].

    Parameters
    ----------
    bc
        BioCypher instance for DB download.
    test
        If `True`, only a subset of the data will be loaded for testing purposes.
    """

    DB_URL = "https://pgx.zju.edu.cn/download.trait/Interactive_TCR-pMHC_Pairs.zip_20250312.zip"

    DB_DIR = "trait_latest"

    def get_latest_release(self, bc: BioCypher) -> str:
        trait_resource = FileDownload(
            name=self.DB_DIR,
            url_s=self.DB_URL,
            lifetime=30,
            is_dir=False,
        )

        trait_parent_path = bc.download(trait_resource)
        trait_path = os.path.join(trait_parent_path[0], os.listdir(trait_parent_path[0])[0])

        if not trait_path:
            raise FileNotFoundError(f"Failed to download TRAIT database from {self.DB_URL}")

        return trait_path

    def read_table(self, bc: BioCypher, table_path: str, test: bool = False) -> pd.DataFrame:
        table = pd.read_excel(table_path)
        if test:
            table = table.sample(frac=0.01, random_state=42)
        # Replace NaN and empty strings with None
        table = table.replace(["", "nan"], None).where(pd.notnull, None)

        rename_cols = {
            "CDR3α": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "CDR3β": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "Epitope": REGISTRY_KEYS.EPITOPE_KEY,
            "Epitope_gene": REGISTRY_KEYS.ANTIGEN_KEY,
            "Epitope_species": REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            "MHC_class": REGISTRY_KEYS.MHC_CLASS_KEY,
            "MHC_A": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "MHC_B": REGISTRY_KEYS.MHC_GENE_2_KEY,
            "TRAV": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "TRAJ": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            "TRBV": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            "TRBJ": REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            "Species": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            # "PubMed.ID": REGISTRY_KEYS.PUBLICATION_KEY,
        }

        table = table.rename(columns=rename_cols)
        table = table[list(rename_cols.values())]
        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY
        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY]

        # Preprocesses CDR3 sequences, epitope sequences, and gene names
        table_preprocessed = harmonize_sequences(bc, table)

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
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.MHC_GENE_2_KEY,
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
