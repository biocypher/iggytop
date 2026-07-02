import pandas as pd
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import get_github_file_last_modified, harmonize_sequences, normalize_table_strings


class BATCAVEAdapter(BaseAdapter):
    """BioCypher adapter for the `BATCAVE <https://github.com/meyer-lab-cshl/BATCAVE-paper>`_ mutational scan database."""

    REPO_NAME = "meyer-lab-cshl/BATMAN-paper"
    """GitHub repository hosting the BATCAVE database files."""
    FILE_PATH_MHCI = "results_batman/tcr_epitope_datasets/mutational_scan_datasets/database/TCR_pMHCI_mutational_scan_database.xlsx"
    """Path (within the repo) to the MHC-I BATCAVE database file."""
    FILE_PATH_MHCII = "results_batman/tcr_epitope_datasets/mutational_scan_datasets/database/TCR_pMHCII_mutational_scan_database.xlsx"
    """Path (within the repo) to the MHC-II BATCAVE database file."""
    RAW_URL_MHCI = f"https://github.com/{REPO_NAME}/raw/main/{FILE_PATH_MHCI}"
    """URL to download the MHC-I BATCAVE database."""
    RAW_URL_MHCII = f"https://github.com/{REPO_NAME}/raw/main/{FILE_PATH_MHCII}"
    """URL to download the MHC-II BATCAVE database."""
    DB_NAME = "BATCAVE"
    """Name of the database."""
    DB_DIR_MHCI = "batcave_mhci_latest"
    """Cache directory name for the MHC-I file."""
    DB_DIR_MHCII = "batcave_mhcii_latest"
    """Cache directory name for the MHC-II file."""
    available_receptors = ["TCR"]
    """Receptor types available in BATCAVE."""

    def get_latest_release(self, bc: BioCypher) -> tuple[str, str]:
        """
        Retrieves the latest release of both BATCAVE database files (MHC-I and MHC-II).

        Args:
            bc: An instance of the BioCypher class.

        Returns:
            Tuple of (mhci_path, mhcii_path).
        """
        mhci_date = get_github_file_last_modified(self.REPO_NAME, self.FILE_PATH_MHCI)
        mhcii_date = get_github_file_last_modified(self.REPO_NAME, self.FILE_PATH_MHCII)
        dates = [d for d in (mhci_date, mhcii_date) if d]
        version = max(dates) if dates else "latest"
        self.set_metadata(version=version, source_url=self.RAW_URL_MHCI)

        mhci_resource = FileDownload(
            name=self.DB_DIR_MHCI,
            url_s=self.RAW_URL_MHCI,
            lifetime=30,
            is_dir=False,
        )
        mhcii_resource = FileDownload(
            name=self.DB_DIR_MHCII,
            url_s=self.RAW_URL_MHCII,
            lifetime=30,
            is_dir=False,
        )

        mhci_path = bc.download(mhci_resource)
        mhcii_path = bc.download(mhcii_resource)

        if not mhci_path:
            raise FileNotFoundError(f"Failed to download BATCAVE MHC-I database from {self.RAW_URL_MHCI}")
        if not mhcii_path:
            raise FileNotFoundError(f"Failed to download BATCAVE MHC-II database from {self.RAW_URL_MHCII}")

        return mhci_path[0], mhcii_path[0]

    def read_table(self, bc: BioCypher, table_path: tuple[str, str], receptors: list[str], test: bool = False) -> pd.DataFrame:
        """
        Reads and processes the BATCAVE tables from both MHC-I and MHC-II database files.

        The BATCAVE database is a mutational scan resource: each TCR is tested against many
        peptide variants of a reference epitope (`index_peptide`). Only the unique
        (TCR, index_peptide) pairings are retained, as these represent the canonical
        TCR-epitope binding events.

        Args:
            bc: An instance of the BioCypher class.
            table_path: Tuple of (mhci_path, mhcii_path).
            receptors: List of receptor types to include (only TCR is available).
            test: If `True`, loads only a subset of the data for testing (default is False).

        Returns:
            A DataFrame containing the processed table data.

        Raises:
            FileNotFoundError: If either table file cannot be found.
        """
        mhci_path, mhcii_path = table_path

        mhci = pd.read_excel(mhci_path)
        mhci[REGISTRY_KEYS.MHC_CLASS_KEY] = "I"

        mhcii = pd.read_excel(mhcii_path)
        mhcii[REGISTRY_KEYS.MHC_CLASS_KEY] = "II"

        table = pd.concat([mhci, mhcii], ignore_index=True)

        # Normalize peptide activity per TCR to its own max activity
        # as done in https://github.com/meyer-lab-cshl/BATMAN-paper/blob/main/results_batman/paper_figures/figure_1/1c/plot_tcr_mutational_scan_heatmaps_for_1c.py
        table["norm_peptide_activity"] = table.groupby("tcr")["peptide_activity"].transform(lambda x: x / x.max())

        # Apply threshold of 0.2 to be sure
        # Batman paper uses 0.1 (see https://www.biorxiv.org/content/10.1101/2024.01.22.576714v3.supplementary-material)
        table = table[table["norm_peptide_activity"] > 0.2]

        if test:
            table = table.sample(frac=0.2, random_state=42)

        table = normalize_table_strings(table)

        # Gene names in BATCAVE lack the locus prefix (e.g. "12-4" → "TRBV12-4")
        # Use masked assignment instead of apply so None values are not silently converted to np.nan
        for col in ["trav", "traj", "trbv", "trbj"]:
            mask = table[col].notna()
            table.loc[mask, col] = col.upper() + table.loc[mask, col]

        rename_cols = {
            "cdr3a": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "trav": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "traj": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            "cdr3b": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "trbv": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            "trbj": REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            "index_peptide": REGISTRY_KEYS.EPITOPE_KEY,  # see below
            "peptide": REGISTRY_KEYS.EPITOPE_KEY + "_mutant",  # see below
            "peptide_type": REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            "mhc": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "pmid": REGISTRY_KEYS.PUBLICATION_KEY,
            "tcr_source_organism": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            "_mhc_class": REGISTRY_KEYS.MHC_CLASS_KEY,
        }

        table = table.rename(columns=rename_cols)
        table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY].str.lower()
        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY]

        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY

        # we use the index peptide for the ideb api wueries as the mutants are related
        table_preprocessed = harmonize_sequences(bc, table)

        # avoid running api calls which do not hit (most epitopes unknown to IEDB)
        table[REGISTRY_KEYS.EPITOPE_KEY] = table[REGISTRY_KEYS.EPITOPE_KEY + "_mutant"]
        table.drop(columns=[REGISTRY_KEYS.EPITOPE_KEY + "_mutant"], inplace=True)

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
            unique_cols=[REGISTRY_KEYS.CHAIN_1_CDR3_KEY],
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
            unique_cols=[REGISTRY_KEYS.CHAIN_2_CDR3_KEY],
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
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
                REGISTRY_KEYS.PUBLICATION_KEY,
            ],
            unique_cols=[REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY],
            property_cols=[
                REGISTRY_KEYS.MHC_CLASS_KEY,
                REGISTRY_KEYS.EPITOPE_KEY,
                REGISTRY_KEYS.EPITOPE_IEDB_ID_KEY,
                REGISTRY_KEYS.MHC_GENE_1_KEY,
                REGISTRY_KEYS.ANTIGEN_KEY,
                REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
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
