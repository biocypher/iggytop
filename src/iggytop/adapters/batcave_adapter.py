import re

import pandas as pd
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .mapping_utils import ASSAY_FUNCTIONAL, ASSAY_SCREEN, ASSAY_UNKNOWN, combine_assays, note_unmapped_assay
from .utils import get_github_file_last_modified, harmonize_sequences, normalize_table_strings

# The functional readout recorded per mutational-scan measurement in BATCAVE's
# `assay` column. Reporter / cytokine / proliferation readouts are functional
# activation; pooled genetic screens (T-Scan family, barcode/minigene) are screens.
_BATCAVE_ASSAY_CATEGORY: dict[str, str] = {
    "nfat luminescence": ASSAY_FUNCTIONAL,
    "nfat-luc2 luminescence": ASSAY_FUNCTIONAL,
    "nfat-gfp": ASSAY_FUNCTIONAL,
    "cd69-gfp": ASSAY_FUNCTIONAL,
    "cd137 expression": ASSAY_FUNCTIONAL,
    "elisa": ASSAY_FUNCTIONAL,
    "elispot": ASSAY_FUNCTIONAL,
    "tnf secretion": ASSAY_FUNCTIONAL,
    "ifng": ASSAY_FUNCTIONAL,
    "t cell proliferation": ASSAY_FUNCTIONAL,
    "t-scan": ASSAY_SCREEN,
    "tscan ii": ASSAY_SCREEN,
    "tcr-map": ASSAY_SCREEN,
    "dna barcode enrichment": ASSAY_SCREEN,
    "minigene depletion": ASSAY_SCREEN,
}

# Some BATCAVE gene calls, e.g. "14-1*00(1200.3)" or "27 (34)", carry IMGT/V-QUEST's own fallback
# format for a call it couldn't resolve to a specific allele: "*00". Stripping it lets tidytcells
# resolve the gene, which V-QUEST already identified even without a specific allele.
_GENE_CONFIDENCE_SUFFIX = re.compile(r"\s*(\*00)?\(\d+(\.\d+)?\)\Z")


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

        # Apply threshold of 0.2 to be sure (find the reasoning behind it here: https://github.com/biocypher/iggytop/pull/55#discussion_r3519732965)
        table = table[table["norm_peptide_activity"] > 0.2]

        if test:
            table = table.sample(frac=0.2, random_state=42)

        table = normalize_table_strings(table)

        for col in ["trav", "traj", "trbv", "trbj"]:
            table[col] = table[col].str.replace(_GENE_CONFIDENCE_SUFFIX, "", regex=True)

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

        # Assay method: look up the per-measurement readout in the `assay` column.
        def _batcave_assay(value):
            if value is None or pd.isna(value):
                return (None, ASSAY_UNKNOWN)
            raw = str(value).strip().lower()
            category = _BATCAVE_ASSAY_CATEGORY.get(raw)
            if category is None:
                note_unmapped_assay(self.DB_NAME, raw)
                category = ASSAY_UNKNOWN
            return combine_assays([(raw, category)])

        _assays = table["assay"].apply(_batcave_assay)
        table[REGISTRY_KEYS.ASSAY_METHOD_RAW_KEY] = _assays.apply(lambda t: t[0])
        table[REGISTRY_KEYS.ASSAY_CATEGORY_KEY] = _assays.apply(lambda t: t[1])

        table = table.rename(columns=rename_cols)
        table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY].str.lower()
        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY]

        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY

        # we use the index peptide for the ideb api wueries as the mutants are related
        table_preprocessed = harmonize_sequences(bc, table)

        # avoid running api calls which do not hit (most epitopes unknown to IEDB)
        table_preprocessed[REGISTRY_KEYS.EPITOPE_KEY] = table_preprocessed[REGISTRY_KEYS.EPITOPE_KEY + "_mutant"]
        table_preprocessed.drop(columns=[REGISTRY_KEYS.EPITOPE_KEY + "_mutant"], inplace=True)

        return table_preprocessed

    def get_nodes(self):
        """Yield BioCypher nodes generated via OntoWeaver."""
        nodes, _ = self._get_ontoweaver_kg()
        yield from nodes

    def get_edges(self):
        """Yield BioCypher edges generated via OntoWeaver."""
        _, edges = self._get_ontoweaver_kg()
        yield from edges
