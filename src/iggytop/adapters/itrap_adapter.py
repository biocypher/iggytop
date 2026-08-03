import pandas as pd
from biocypher import BioCypher, FileDownload

from .base_adapter import BaseAdapter
from .constants import REGISTRY_KEYS
from .utils import get_github_file_last_modified, harmonize_sequences, normalize_table_strings


class ITRAPAdapter(BaseAdapter):
    """BioCypher adapter for the ITRAP benchmark dataset.

    Data from: https://github.com/mnielLab/ITRAP_benchmark
    """

    REPO_NAME = "mnielLab/ITRAP_benchmark"
    """GitHub repository hosting the ITRAP database file."""
    FILE_PATH = "ITRAP.csv"
    """Path (within the repo) to the ITRAP database file."""
    DB_URL = f"https://raw.githubusercontent.com/{REPO_NAME}/main/{FILE_PATH}"
    """URL to download the ITRAP database."""
    DB_DIR = "itrap_latest"
    """Directory name for the downloaded database."""
    DB_NAME = "ITRAP"
    """Name of the database."""
    available_receptors = ["TCR"]
    """Receptor types available in ITRAP."""

    def get_latest_release(self, bc: BioCypher) -> str:
        """
        Retrieves the latest release of the ITRAP database.

        Args:
            bc: An instance of the BioCypher class.

        Returns:
            Path to the latest release file.
        """
        version = get_github_file_last_modified(self.REPO_NAME, self.FILE_PATH) or "latest"
        self.set_metadata(version=version, source_url=self.DB_URL)
        itrap_resource = FileDownload(
            name=self.DB_DIR,
            url_s=self.DB_URL,
            lifetime=30,
            is_dir=False,
        )

        itrap_path = bc.download(itrap_resource)

        if not itrap_path:
            raise FileNotFoundError(f"Failed to download ITRAP database from {self.DB_URL}")

        return itrap_path[0]

    def read_table(self, bc: BioCypher, table_path: str, receptors: list[str], test: bool = False) -> pd.DataFrame:
        """
        Reads and processes the ITRAP table from the downloaded database file.

        Args:
            bc: An instance of the BioCypher class.
            table_path: Path to the table file.
            receptors: List of receptor types to include in the table (Ignored as only TCR is available).
            test: If `True`, loads only a subset of the data for testing (default is False).

        Returns:
            A DataFrame containing the processed table data.
        """
        table = pd.read_csv(table_path)

        if test:
            table = table.sample(frac=0.1, random_state=42)

        table = normalize_table_strings(table)

        # Mapping between peptide and antigen/source information based on 10X genomics documentation
        # https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_AN047_IP_A_New_Way_of_Exploring_Immunity_Digital.pdf
        antigen_mapping = {
            "VTEHDTLLY": ("IE-1", "CMV"),
            "KTWGQYWQV": ("gp100", "Cancer"),
            "ELAGIGILTV": ("MART-1", "Cancer"),
            "CLLWSFQTSA": ("Tyrosinase", "Cancer"),
            "IMDQVPFSV": ("gp100", "Cancer"),
            "SLLMWITQV": ("NY-ESO-1", "Cancer"),
            "KVAELVHFL": ("MAGE A3", "Cancer"),
            "KVLEYVIKV": ("MAGE-A1", "Cancer"),
            "CLLGTYTQDV": ("Kanamycin B dioxygenase", "Bacteria"),
            "LLDFVRFMGV": ("EBNA 3B", "EBV"),
            "LLMGTLGIVC": ("HPV 16E7", "HPV"),
            "CLGGLLTMV": ("LMP-2A", "EBV"),
            "YLLEMLWRL": ("LMP1", "EBV"),
            "FLYALALLL": ("LMP2A", "EBV"),
            "GILGFVFTL": ("Flu MP", "Influenza"),
            "GLCTLVAML": ("BMLF1", "EBV"),
            "NLVPMVATV": ("pp65", "CMV"),
            "ILKEPVHGV": ("RT", "HIV"),
            "FLASKIGRLV": ("Ca2+-indepen. Plip A2", "Homo sapiens"),
            "CYTWNQMNL": ("WT1", "Homo sapiens"),
            "RTLNAWVKV": ("Gag protein", "HIV"),
            "KLQCVDLHV": ("PSA", "Homo sapiens"),
            "LLFGYPVYV": ("HTLV-1", "HTLV-1"),
            "SLFNTVATL": ("Gag protein", "HIV"),
            "SLYNTVATLY": ("Gag protein", "HIV"),
            "SLFNTVATLY": ("Gag protein", "HIV"),
            "RMFPNAPYL": ("WT-1", "Homo sapiens"),
            "YLNDHLEPWI": ("BCL-X", "Homo sapiens"),
            "MLDLQPETT": ("16E7", "HPV"),
            "KLGGALQAK": ("IE-1", "CMV"),
            "RLRAEAQVK": ("EMNA 3A", "EBV"),
            "RIAAWMATY": ("BCL-2L1", "Homo sapiens"),
            "IVTDFSVIK": ("EBNA 3B", "EBV"),
            "AVFDRKSDAK": ("EBNA 3B", "EBV"),
            "IPSINVHHY": ("pp65", "CMV"),
            "AYAQKIFKI": ("IE-1", "CMV"),
            "QYDPVAALF": ("pp65", "CMV"),
            "QPRAPIRPI": ("EBNA 6", "EBV"),
            "TPRVTGGGAM": ("pp65", "CMV"),
            "RPPIFIRRL": ("EBNA 3A", "EBV"),
            "RPHERNGFTVL": ("pp65", "CMV"),
            "RAKFKQLL": ("BZLF1", "EBV"),
            "ELRRKMMYM": ("IE-1", "CMV"),
            "FLRGRAYGL": ("EBNA 3A", "EBV"),
            # Negative controls
            "SLEGGGLGY": ("Negative Control", "NC"),
            "STEGGGLAY": ("Negative Control", "NC"),
            "ALIAPVHAV": ("Negative Control", "NC"),
            "AYSSAGASI": ("Negative Control", "NC"),
            "GPAESAAGL": ("Negative Control", "NC"),
            "AAKGRGAAL": ("Negative Control", "NC"),
        }

        # Split peptide_HLA (e.g., "GILGFVFTL HLA-A*02:01")
        # Assuming format is PEPTIDE_HLA
        def split_peptide_hla(val):
            if pd.isna(val) or " " not in val:
                return val, None
            parts = val.split(" ")
            peptide = parts[0]
            hla = " ".join(parts[1:])
            return peptide, hla

        table[["peptide", "hla"]] = table["peptide_HLA"].apply(lambda x: pd.Series(split_peptide_hla(x)))

        rename_cols = {
            "peptide": REGISTRY_KEYS.EPITOPE_KEY,
            "hla": REGISTRY_KEYS.MHC_GENE_1_KEY,
            "cdr3_TRA": REGISTRY_KEYS.CHAIN_1_CDR3_KEY,
            "v_gene_TRA": REGISTRY_KEYS.CHAIN_1_V_GENE_KEY,
            "j_gene_TRA": REGISTRY_KEYS.CHAIN_1_J_GENE_KEY,
            "cdr3_TRB": REGISTRY_KEYS.CHAIN_2_CDR3_KEY,
            "v_gene_TRB": REGISTRY_KEYS.CHAIN_2_V_GENE_KEY,
            "j_gene_TRB": REGISTRY_KEYS.CHAIN_2_J_GENE_KEY,
        }

        table = table.rename(columns=rename_cols)

        # Add metadata
        table[REGISTRY_KEYS.CHAIN_1_TYPE_KEY] = REGISTRY_KEYS.TRA_KEY
        table[REGISTRY_KEYS.CHAIN_2_TYPE_KEY] = REGISTRY_KEYS.TRB_KEY
        table[REGISTRY_KEYS.CHAIN_1_ORGANISM_KEY] = "Homo sapiens"
        table[REGISTRY_KEYS.CHAIN_2_ORGANISM_KEY] = "Homo sapiens"
        table[REGISTRY_KEYS.TISSUE_KEY] = "PBMC"
        table[REGISTRY_KEYS.PUBLICATION_KEY] = "PMID: 37133356; 10XGenomics"

        # Apply antigen mapping
        def apply_mapping(peptide):
            if peptide in antigen_mapping:
                return antigen_mapping[peptide]
            return (None, None)

        mapped_info = table[REGISTRY_KEYS.EPITOPE_KEY].apply(apply_mapping)
        table[REGISTRY_KEYS.ANTIGEN_KEY] = mapped_info.apply(lambda x: x[0])
        table[REGISTRY_KEYS.ANTIGEN_ORGANISM_KEY] = mapped_info.apply(lambda x: x[1])

        # Preprocess and harmonize
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
