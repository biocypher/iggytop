from typing import NamedTuple


class _REGISTRY_KEYS_NT(NamedTuple):
    EPITOPE_KEY: str = "epitope_sequence"
    ANTIGEN_KEY: str = "antigen_name"
    ANTIGEN_ORGANISM_KEY: str = "antigen_organism"
    CHAIN_1_CDR1_KEY: str = "chain_1_cdr1"
    CHAIN_1_CDR2_KEY: str = "chain_1_cdr2"
    CHAIN_1_CDR3_KEY: str = "chain_1_cdr3"
    CHAIN_2_CDR1_KEY: str = "chain_2_cdr1"
    CHAIN_2_CDR2_KEY: str = "chain_2_cdr2"
    CHAIN_2_CDR3_KEY: str = "chain_2_cdr3"
    CHAIN_1_V_GENE_KEY: str = "chain_1_v_gene"
    CHAIN_1_J_GENE_KEY: str = "chain_1_j_gene"
    CHAIN_2_V_GENE_KEY: str = "chain_2_v_gene"
    CHAIN_2_J_GENE_KEY: str = "chain_2_j_gene"
    CHAIN_1_ORGANISM_KEY: str = "chain_1_organism"
    CHAIN_2_ORGANISM_KEY: str = "chain_2_organism"
    TRA_KEY: str = "tra"
    TRB_KEY: str = "trb"
    IGH_KEY: str = "igh"
    IGL_KEY: str = "igl"
    TYPE_KEY: str = "type"


REGISTRY_KEYS = _REGISTRY_KEYS_NT()
