from typing import NamedTuple


class _REGISTRY_KEYS_NT(NamedTuple):
    EPITOPE_KEY: str = "epitope_sequence"
    EPITOPE_IEDB_ID_KEY: str = "iedb_id"
    ANTIGEN_KEY: str = "antigen_name"
    ANTIGEN_ORGANISM_KEY: str = "antigen_organism"
    CHAIN_1_CDR1_KEY: str = "chain_1_cdr1"
    CHAIN_1_CDR2_KEY: str = "chain_1_cdr2"
    CHAIN_1_CDR3_KEY: str = "chain_1_cdr3"
    CHAIN_2_CDR1_KEY: str = "chain_2_cdr1"
    CHAIN_2_CDR2_KEY: str = "chain_2_cdr2"
    CHAIN_2_CDR3_KEY: str = "chain_2_cdr3"
    CHAIN_1_V_GENE_KEY: str = "chain_1_v_call"
    CHAIN_1_J_GENE_KEY: str = "chain_1_j_call"
    CHAIN_2_V_GENE_KEY: str = "chain_2_v_call"
    CHAIN_2_J_GENE_KEY: str = "chain_2_j_call"
    CHAIN_1_ORGANISM_KEY: str = "chain_1_organism"
    CHAIN_2_ORGANISM_KEY: str = "chain_2_organism"
    CHAIN_1_TYPE_KEY: str = "chain_1_type"
    CHAIN_2_TYPE_KEY: str = "chain_2_type"
    TRA_KEY: str = "tra"
    TRB_KEY: str = "trb"
    IGH_KEY: str = "igh"
    IGL_KEY: str = "igl"
    MHC_CLASS_KEY: str = "MHC_class"
    MHC_GENE_1_KEY: str = "mhc_gene_1"
    MHC_GENE_2_KEY: str = "mhc_gene_2"


REGISTRY_KEYS = _REGISTRY_KEYS_NT()
