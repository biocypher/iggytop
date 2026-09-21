"""Unit tests for assay-method harmonization (no network required)."""

from iggytop.adapters.mapping_utils import (
    ASSAY_CATEGORIES,
    canonicalize_assay_token,
    combine_assays,
    harmonize_assay_iedb,
    harmonize_assay_verbs,
)


def test_canonicalize_assay_token_merges_typos_and_case():
    assert canonicalize_assay_token("Tetramer.sort") == "tetramer-sort"
    assert canonicalize_assay_token("tetramer sort") == "tetramer-sort"
    assert canonicalize_assay_token("Cultured-T-cells") == "cultured-t-cells"
    assert canonicalize_assay_token("limiting-diffusion-cloning") == "limiting-dilution-cloning"
    assert canonicalize_assay_token("antigen-loaded-targed") == "antigen-loaded-targets"
    # uninformative tokens are dropped
    assert canonicalize_assay_token("cla") is None
    assert canonicalize_assay_token("") is None


def test_harmonize_assay_verbs_single_and_multi():
    assert harmonize_assay_verbs("tetramer-sort") == ("tetramer-sort", "multimer_binding")
    # comma-separated verbs -> sorted, de-duplicated raw + combined categories
    raw, cat = harmonize_assay_verbs("tetramer.sort,dextramer-sort")
    assert raw == "dextramer-sort|tetramer-sort"
    assert cat == "multimer_binding"
    raw, cat = harmonize_assay_verbs("peptide-stimulation,tetramer-sort")
    assert cat == "functional_activation|multimer_binding"


def test_harmonize_assay_verbs_unknown_and_empty():
    assert harmonize_assay_verbs(None) == (None, "unknown")
    assert harmonize_assay_verbs("") == (None, "unknown")
    # an unrecognized verb is kept as raw but categorized unknown
    assert harmonize_assay_verbs("some-new-assay") == ("some-new-assay", "unknown")


def test_harmonize_assay_verbs_list_input_spans_columns():
    # TRAIT feeds identification + affinity + structure columns together
    raw, cat = harmonize_assay_verbs(["Dextramer-Sort", "SPR", "XRD"])
    assert raw == "dextramer-sort|spr|x-ray crystallography"
    assert cat == "biophysical_affinity|multimer_binding|structural"


def test_harmonize_assay_iedb_reads_technique_field():
    assert harmonize_assay_iedb(["qualitative binding|multimer/tetramer"]) == (
        "multimer/tetramer",
        "multimer_binding",
    )
    assert harmonize_assay_iedb(["dissociation constant KD|surface plasmon resonance (SPR)|nM"])[1] == "biophysical_affinity"
    assert harmonize_assay_iedb(["3D structure|x-ray crystallography|angstroms"])[1] == "structural"
    assert harmonize_assay_iedb(["T cell binding|High throughput multiplexed assay"])[1] == "high_throughput_screen"
    # several assays on one record combine
    raw, cat = harmonize_assay_iedb(["IFNg release|ELISA", "qualitative binding|multimer/tetramer"])
    assert cat == "functional_activation|multimer_binding"
    assert harmonize_assay_iedb([]) == (None, "unknown")


def test_combine_assays_primitive():
    # used by the MCPAS / BATCAVE adapters, which own their own vocabulary
    assert combine_assays([("pmhc-multimer", "multimer_binding")]) == ("pmhc-multimer", "multimer_binding")
    raw, cat = combine_assays([("t-scan", "high_throughput_screen"), ("elisa", "functional_activation")])
    assert raw == "elisa|t-scan"
    assert cat == "functional_activation|high_throughput_screen"
    assert combine_assays([(None, "unknown")]) == (None, "unknown")


def test_every_category_emitted_is_declared():
    for value in ("tetramer-sort", "some-new-assay", None):
        _, cat = harmonize_assay_verbs(value)
        assert set(cat.split("|")) <= set(ASSAY_CATEGORIES)
