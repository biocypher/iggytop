import pandas as pd
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .mapping_utils import (
    ASSAY_FUNCTIONAL,
    ASSAY_MULTIMER,
    ASSAY_TISSUE,
    ASSAY_UNKNOWN,
    combine_assays,
    note_unmapped_assay,
)
from .utils import harmonize_sequences, normalize_table_strings

# McPAS encodes the assay in the numeric "Antigen.identification.method" column:
# 1 = pMHC multimer, 2.x = in-vitro stimulation (by antigen form), 3 = isolated
# from disease tissue with no antigen-specific selection, 4 = unspecified.
# See https://friedmanlab.weizmann.ac.il/McPAS-TCR/.
_MCPAS_METHOD: dict[str, tuple[str, str]] = {
    "1": ("pmhc-multimer", ASSAY_MULTIMER),
    "2": ("in-vitro-stimulation", ASSAY_FUNCTIONAL),
    "2.1": ("in-vitro-stimulation-peptide", ASSAY_FUNCTIONAL),
    "2.2": ("in-vitro-stimulation-protein", ASSAY_FUNCTIONAL),
    "2.3": ("in-vitro-stimulation-pathogen", ASSAY_FUNCTIONAL),
    "2.4": ("in-vitro-stimulation-tumor-cells", ASSAY_FUNCTIONAL),
    "2.5": ("in-vitro-stimulation-other", ASSAY_FUNCTIONAL),
    "3": ("tissue-isolation", ASSAY_TISSUE),
    "4": ("unspecified", ASSAY_UNKNOWN),
}


class MCPASAdapter(BaseAdapter):
    """BioCypher adapter for the manually-curated catalogue of pathology-associated T cell
    receptor sequences `McPAS-TCR <https://friedmanlab.weizmann.ac.il/McPAS-TCR.csv>`_.

    """

    DB_URL = "https://friedmanlab.weizmann.ac.il/McPAS-TCR.csv"
    """URL to download the McPAS-TCR database."""
    DB_DIR = "mcpas_latest"
    """Directory name for the downloaded database."""
    DB_NAME = "MCPAS"
    """Name of the database."""
    available_receptors = ["TCR"]
    """Receptor types available in McPAS-TCR."""

    def get_latest_release(self, bc: BioCypher) -> str:
        """
        Retrieves the latest release of the McPAS-TCR database.

        Args:
            bc: An instance of the BioCypher class.

        Returns:
            Path to the latest release file.
        """
        self.set_metadata(version="latest", source_url=self.DB_URL)
        mcpas_resource = FileDownload(
            name=self.DB_DIR,
            url_s=self.DB_URL,
            lifetime=30,
            is_dir=False,
        )

        mcpas_path = bc.download(mcpas_resource)

        if not mcpas_path:
            raise FileNotFoundError(f"Failed to download MCPAS-TCR database from {self.DB_URL}")

        # mcpas_path = "../data/MCPAS-TCR1.csv"

        return mcpas_path[0]

    def read_table(self, bc: BioCypher, table_path: str, receptors: list[str], test: bool = False) -> pd.DataFrame:
        """
        Reads and processes the MCPAS table from the downloaded database file.

        Args:
            bc: An instance of the BioCypher class.
            table_path: Path to the table file.
            receptors: List of receptor types to include in the table. (Ignored as only TCR is available).
            test: If `True`, loads only a subset of the data for testing (default is False).

        Returns:
            A DataFrame containing the processed table data.

        Raises:
            FileNotFoundError: If the table file cannot be found.
        """
        table = pd.read_csv(table_path, encoding="utf-8-sig")
        if test:
            table = table.sample(frac=0.001, random_state=42)
        table = normalize_table_strings(table)

        table["Pathology"] = table.apply(
            lambda row: "HomoSapiens" if row["Category"] == "Autoimmune" else row["Pathology"],
            axis=1,
        )

        rename_cols = {
            "CDR3.alpha.aa": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "CDR3.beta.aa": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "Epitope.peptide": REGISTRY_KEYS.EPITOPE_KEY,
            "Antigen.protein": REGISTRY_KEYS.ANTIGEN_KEY,
            "Pathology": REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY,
            "MHC": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "T.Cell.Type": REGISTRY_KEYS.MHC_CLASS_KEY,
            "TRAV": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "TRAJ": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            "TRBV": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            "TRBJ": REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
            "Species": REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY,
            "Tissue": REGISTRY_KEYS.TISSUE_KEY,
            "PubMed.ID": REGISTRY_KEYS.PUBLICATION_KEY,
        }

        # Assay method: map the numeric identification code (arrives as "1.0" / "2.2").
        def _mcpas_method(code):
            if code is None or pd.isna(code):
                return (None, ASSAY_UNKNOWN)
            key = str(code).strip()
            key = key[:-2] if key.endswith(".0") else key
            if key not in _MCPAS_METHOD:
                note_unmapped_assay(self.DB_NAME, f"identification-method-code {key}")
                return (None, ASSAY_UNKNOWN)
            return combine_assays([_MCPAS_METHOD[key]])

        _assays = table["Antigen.identification.method"].apply(_mcpas_method)
        table[REGISTRY_KEYS.ASSAY_METHOD_RAW_KEY] = _assays.apply(lambda t: t[0])
        table[REGISTRY_KEYS.ASSAY_CATEGORY_KEY] = _assays.apply(lambda t: t[1])
        rename_cols[REGISTRY_KEYS.ASSAY_METHOD_RAW_KEY] = REGISTRY_KEYS.ASSAY_METHOD_RAW_KEY
        rename_cols[REGISTRY_KEYS.ASSAY_CATEGORY_KEY] = REGISTRY_KEYS.ASSAY_CATEGORY_KEY

        table = table.rename(columns=rename_cols)
        table = table[list(rename_cols.values())]
        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY
        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY]

        # Preprocesses CDR3 sequences, epitope sequences, and gene names
        table_preprocessed = harmonize_sequences(bc, table)

        return table_preprocessed

    def get_nodes(self):
        """Yield BioCypher nodes generated via OntoWeaver."""
        nodes, _ = self._get_ontoweaver_kg()
        yield from nodes

    def get_edges(self):
        """Yield BioCypher edges generated via OntoWeaver."""
        _, edges = self._get_ontoweaver_kg()
        yield from edges
