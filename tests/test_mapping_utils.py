"""Unit tests for assay-method harmonization (no network required)."""

from iggytop.adapters.mapping_utils import (
    ASSAY_CATEGORIES,
    assay_method_categories,
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
    # a technique that decides the category on its own is not qualified with its readout
    assert harmonize_assay_iedb(["dissociation constant KD|surface plasmon resonance (SPR)|nM"])[0] == "surface plasmon resonance (spr)"
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


def test_harmonize_assay_iedb_qualifies_a_generic_technique():
    # "binding assay" means different things depending on its readout, so the two senses
    # must not collapse onto one method string: a method carries exactly one category
    # through the knowledge graph, and the shared node would inherit the concrete one.
    affinity_raw, affinity_cat = harmonize_assay_iedb(["dissociation constant KD|binding assay|nM"])
    assert affinity_raw == "binding assay (dissociation constant kd)"
    assert affinity_cat == "biophysical_affinity"

    plain_raw, plain_cat = harmonize_assay_iedb(["qualitative binding|binding assay"])
    assert plain_raw == "binding assay"
    assert plain_cat == "unknown"

    assert affinity_raw != plain_raw
    assert assay_method_categories(affinity_raw) == {"binding assay (dissociation constant kd)": "biophysical_affinity"}
    assert assay_method_categories(plain_raw) == {"binding assay": "unknown"}

    # a readout-only record (no technique component) is not double-qualified
    assert harmonize_assay_iedb(["dissociation constant KD"]) == ("dissociation constant kd", "biophysical_affinity")


def test_assay_method_categories_keeps_the_pairing():
    # The two harmonized columns are independent joined sets, so the knowledge graph
    # recovers each method's own category from the method string alone.
    raw, cat = harmonize_assay_verbs("peptide-stimulation,tetramer-sort")
    assert raw == "peptide-stimulation|tetramer-sort"
    assert cat == "functional_activation|multimer_binding"
    assert assay_method_categories(raw) == {
        "peptide-stimulation": "functional_activation",
        "tetramer-sort": "multimer_binding",
    }
    # a method no source harmonized, and one that was but stayed uncategorized
    assert assay_method_categories("never-seen-this") == {"never-seen-this": "unknown"}
    raw, _ = harmonize_assay_verbs("some-new-assay")
    assert assay_method_categories(raw) == {"some-new-assay": "unknown"}
    # a record whose source reports no method at all contributes nothing
    for empty in (None, "", "  ", "nan"):
        assert assay_method_categories(empty) == {}


def test_assay_method_categories_prefers_a_concrete_category():
    # Registration order must not matter: a concrete category always wins over `unknown`.
    combine_assays([("ambiguous-assay", "unknown")])
    combine_assays([("ambiguous-assay", "structural")])
    combine_assays([("ambiguous-assay", "unknown")])
    assert assay_method_categories("ambiguous-assay") == {"ambiguous-assay": "structural"}


def test_every_category_emitted_is_declared():
    for value in ("tetramer-sort", "some-new-assay", None):
        _, cat = harmonize_assay_verbs(value)
        assert set(cat.split("|")) <= set(ASSAY_CATEGORIES)
