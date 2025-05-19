import pandas as pd
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import get_iedb_ids_batch, process_sequence


class MCPASAdapter(BaseAdapter):
    """BioCypher adapter for the manually-curated catalogue of pathology-associated T cell
    receptor sequences (McPAS-TCR)[https://friedmanlab.weizmann.ac.il/McPAS-TCR.csv].

    Parameters
    ----------
    bc
        BioCypher instance for DB download.
    test
        If `True`, only a subset of the data will be loaded for testing purposes.
    """

    DB_URL = "https://friedmanlab.weizmann.ac.il/McPAS-TCR.csv"
    DB_DIR = "mcpas_latest"

    def get_latest_release(self, bc: BioCypher) -> str:
        mcpas_resource = FileDownload(
            name=self.DB_DIR,
            url_s=self.DB_URL,
            lifetime=30,
            is_dir=False,
        )

        mcpas_path = bc.download(mcpas_resource)

        if not mcpas_path:
            raise FileNotFoundError(f"Failed to download McPAS-TCR database from {self.DB_URL}")

        return mcpas_path[0]

    def read_table(self, bc: BioCypher, table_path: str, test: bool = False) -> pd.DataFrame:
        table = pd.read_csv(table_path, encoding="utf-8-sig")
        if test:
            table = table.sample(frac=0.01, random_state=42)
        # Replace NaN and empty strings with None
        table = table.replace(["", "nan"], None).where(pd.notnull, None)

        rename_cols = {
            "CDR3.alpha.aa": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "CDR3.beta.aa": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "Epitope.peptide": REGISTRY_KEYS.EPITOPE_KEY,
            "Antigen.protein": REGISTRY_KEYS.ANTIGEN_KEY,
            "MHC": REGISTRY_KEYS.MHC_GENE_1_KEY,
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

        sequence_cols = [
            REGISTRY_KEYS.EPITOPE_KEY,
            REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
        ]

        for col in sequence_cols:
            table[col] = table[col].apply(process_sequence)

        # Map epitope sequences to IEDB IDs
        valid_epitopes = table[REGISTRY_KEYS.EPITOPE_KEY].dropna().drop_duplicates().tolist()
        if len(valid_epitopes) > 0:
            epitope_map = get_iedb_ids_batch(bc, valid_epitopes)
        # Apply the mapping to create the IEDB ID column
        table[REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY] = table[REGISTRY_KEYS.EPITOPE_KEY].map(epitope_map)

        return table

    def get_nodes(self):
        # chain 1
        yield from self._generate_nodes_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            ],
            unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
        )

        # chain 2
        yield from self._generate_nodes_from_table(
            [
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
                REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
                REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY,
            ],
            unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
        )

        # epitope
        yield from self._generate_nodes_from_table(
            [
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            ],
            unique_cols=REGISTRY_KEYS.EPITOPE_KEY,
        )

    def get_edges(self):
        # chain 1 to chain 2
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            ],
            [
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            ],
            source_unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
        )

        # chain 1 to epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_1_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            ],
            [REGISTRY_KEYS.EPITOPE_KEY],
            source_unique_cols=REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.EPITOPE_KEY,
        )

        # chain 2 to epitope
        yield from self._generate_edges_from_table(
            [
                REGISTRY_KEYS.CHAIN_2_TYPE_KEY,
                REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            ],
            [REGISTRY_KEYS.EPITOPE_KEY],
            source_unique_cols=REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            target_unique_cols=REGISTRY_KEYS.EPITOPE_KEY,
        )
