"""This module contains utility functions for cleaning up species terms and cleaning antigen names"""

import sys

sys.path.append("..")

import os
import re
from datetime import datetime
from urllib.parse import quote

import requests


def map_species_terms(terms: list[str], zooma: bool = False) -> dict:
    """Harmonize and normalize species terms using manual mappings and Zooma API.
    Args:
        terms: List of species terms to normalize.
        zooma: If True, use Zooma API to get labels for normalized terms.
    Returns:
        A dictionary mapping original terms to normalized terms.
    """
    terms = [x for x in terms if x is not None]
    manual_disambiguation = {
        "AdV": "Human adenovirus",
        "CMV": "Cytomegalovirus",
        "DENV": "Dengue virus",
        "EBV": "Epstein-Barr virus",
        "HCV": "Hepatitis C virus",
        "HHV": "Human herpesvirus",
        "HIV": "Human immunodeficiency virus",
        "HPV": "Human papillomavirus",
        "HTLV": "Human T-cell leukemia virus",
        "HSV": "Herpes simplex virus",
        "InfluenzaA": "Influenza A virus",
        "LCMV": "Lymphocytic choriomeningitis virus",
        "MCPyV": "Merkel cell polyomavirus",
        "McpyV": "Merkel cell polyomavirus",
        "Mtb": "Mycobacterium tuberculosis",
        "SARS-CoV1": "Severe acute respiratory syndrome coronavirus",
        "SARS-CoV2": "Severe acute respiratory syndrome coronavirus 2",
        "SARS-CoV": "Severe acute respiratory syndrome coronavirus",
        "SIV": "Simian immunodeficiency virus",
        "YFV": "Yellow fever virus",
        "Mouse": "Mus musculus",
        "Human": "Homo sapiens",
        "autoimmune": "Homo sapiens",
    }

    def normalize_species(term: str) -> str:
        """Normalize species terms by applying manual mappings for abbreviations,
        cleaning up formatting"""
        if term.startswith("http"):
            label, iri = get_label_from_semantic_tag(term)
            return label
        term = term.strip()
        term = re.sub(r"^([a-zA-Z]+)(\d+)(?![a-zA-Z])", r"\1 \2", term)
        term_lower = term.lower()
        for prefix in manual_disambiguation:
            if term_lower.startswith(prefix.lower()):
                suffix = term[len(prefix) :]
                query_term = manual_disambiguation[prefix] + suffix
                break
        else:
            query_term = term

        # Replace common separators and clean up
        query_term = query_term.replace("_", " ")
        query_term = re.sub(r"([-_/])(?=\d)", " ", query_term)
        query_term = query_term[0].upper() + query_term[1:]

        # Remove any content in parentheses or brackets and trailing strain
        query_term = re.sub(r"\s*[\(\[].*[\)\]]", "", query_term)
        # query_term = re.sub(r"\bstrain\s.*", "", query_term).strip()
        query_term = re.sub(r"\b(strain|str\.|subsp\.|variant|genotype)\s+[^\s]+", "", query_term, flags=re.IGNORECASE)

        if "-" not in query_term:
            query_term = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", query_term)

        if (
            "severe acute respiratory syndrome coronavirus 2" in query_term.lower()
            or "severe acute respiratory coronavirus 2" in query_term.lower()
        ):
            query_term = "Severe acute respiratory syndrome coronavirus 2"

        words = query_term.split()
        if words:
            normalized_words = [words[0]]
            for word in words[1:]:
                if not word.isupper() or not word.isalpha():
                    normalized_words.append(word.lower())
                else:
                    normalized_words.append(word)
        return " ".join(normalized_words)

    def get_label_from_semantic_tag(uri: str):
        """
        Get label from semantic tag URI, supporting both OBO and IEDB ontologies.

        Args:
            uri: The URI to process (e.g., 'http://purl.obolibrary.org/obo/NCBITaxon_9838'
                or 'https://ontology.iedb.org/ontology/ONTIE_0000884')

        Returns:
            tuple: (label, full_uri)
        """
        try:
            # Handle OBO format (purl.obolibrary.org)
            if "obo/" in uri:
                term = uri.split("obo/")[-1]
                ontology = term.split("_")[0].lower()
                full_uri = f"http://purl.obolibrary.org/obo/{term}"
                encoded_uri = quote(quote(full_uri, safe=""), safe="")
                ols_url = f"https://www.ebi.ac.uk/ols4/api/ontologies/{ontology}/terms/{encoded_uri}"

                res = requests.get(ols_url, timeout=10)
                res.raise_for_status()
                label = res.json().get("label")
                return label, full_uri

            # Handle IEDB format (ontology.iedb.org)
            elif "ontology.iedb.org/ontology/" in uri:
                full_uri = uri  # Use the original URI as-is
                # IEDB uses direct JSON-LD API, just append .json to the term IRI
                iedb_url = f"{uri}.json"

                res = requests.get(iedb_url, timeout=10)
                res.raise_for_status()
                data = res.json()
                label = data.get("rdfs:label")
                return label, full_uri

            else:
                return None, None
        except Exception:
            return None, None

    def get_zooma_label(term: str):
        """Get label for a species term using the Zooma API and the following parameters:
        - propertyType: "organism"
        - sources: "uniprot"
        - ontologies: "ncbitaxon"
        - accepted confidence: "HIGH" or "GOOD"
        Zooma API first checks the sources for match and then, checks the ontologies
        """
        zooma_url = "https://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate"
        sources = ["uniprot"]
        ontologies = ["ncbitaxon"]
        params = {
            "propertyValue": term,
            "propertyType": "organism",
            "ontologies": f"[{','.join(ontologies)}]",
            "filter": f"required:[{','.join(sources)}],ontologies:[{','.join(ontologies)}]",
        }
        try:
            r = requests.get(zooma_url, params=params, timeout=10)
            r.raise_for_status()
            results = r.json()
        except Exception:
            return term

        for r in results:
            if r.get("confidence", "").upper() in {"HIGH", "GOOD"}:
                tags = r.get("semanticTags", [])
                if tags:
                    label, iri = get_label_from_semantic_tag(tags[0])
                    if label:
                        return label
                    else:
                        return term
        return None

    # Step 1: Normalize all terms (deduplicate first to avoid redundant API calls for http URIs)
    unique_terms = set(t for t in terms if t)
    normalized_unique = {term: normalize_species(term) for term in unique_terms}
    normalized_terms = {term: normalized_unique[term] for term in terms if term}

    # print("Normalized terms:", normalized_terms)
    results = {}

    if zooma:
        # Step 2: Get Zooma mappings for normalized terms
        for original_term, normalized_term in normalized_terms.items():
            zooma_result = get_zooma_label(normalized_term)
            # Create final results - use Zooma output if available, otherwise use normalized term
            if zooma_result is not None:
                results[original_term] = zooma_result
            else:
                results[original_term] = normalized_terms[original_term]
    else:
        results = normalized_terms

    return results


# ---------------------------------------------------------------------------
# Assay-method harmonization
#
# Every source records "how was this receptor:epitope pairing established" in its
# own vocabulary. We collapse those onto two harmonized columns:
#   - assay_method_raw: the source verb(s), lightly canonicalized (typos/case
#     merged) and "|"-joined when a record has several.
#   - assay_category:   one or more evidence classes (see ASSAY_CATEGORIES).
# ---------------------------------------------------------------------------

ASSAY_MULTIMER = "multimer_binding"
ASSAY_AFFINITY = "biophysical_affinity"
ASSAY_STRUCTURAL = "structural"
ASSAY_SCREEN = "high_throughput_screen"
ASSAY_FUNCTIONAL = "functional_activation"
ASSAY_CLONAL = "clonal_isolation"
ASSAY_TISSUE = "tissue_isolation"
ASSAY_UNKNOWN = "unknown"

ASSAY_CATEGORIES = (
    ASSAY_MULTIMER,
    ASSAY_AFFINITY,
    ASSAY_STRUCTURAL,
    ASSAY_SCREEN,
    ASSAY_FUNCTIONAL,
    ASSAY_CLONAL,
    ASSAY_TISSUE,
    ASSAY_UNKNOWN,
)

# Spelling / case / punctuation variants collapsed to one canonical verb.
# Applied to the free-text identification verbs of VDJDB and TRAIT (they share a
# vocabulary). Keys are matched case-insensitively against stripped tokens.
_ASSAY_TOKEN_CANONICAL = {
    "tetramer.sort": "tetramer-sort",
    "tetramer sort": "tetramer-sort",
    "cultured-t-cells": "cultured-t-cells",
    "cultivated-t-cells": "cultured-t-cells",
    "antigen-loaded-target": "antigen-loaded-targets",
    "antigen-loaded-targed": "antigen-loaded-targets",
    "antigen loaded targets": "antigen-loaded-targets",
    "antigen-expressing cells": "antigen-expressing-targets",
    "limiting-diffusion-cloning": "limiting-dilution-cloning",
    "limited-dilution-cloning": "limiting-dilution-cloning",
    "cloning-by-limiting-dilution": "limiting-dilution-cloning",
    "tetramers staining": "tetramer-stain",
    "tetramer-staining": "tetramer-stain",
    "tetramer staining": "tetramer-stain",
    # affinity / structure method labels (TRAIT's Affinity_method / Structure_method)
    "xrd": "x-ray crystallography",
    "x-ray": "x-ray crystallography",
    "xray": "x-ray crystallography",
    "surface plasmon resonance": "spr",
    "surface plasmon resonance (spr)": "spr",
}

# Tokens carrying no method information -> dropped entirely.
_ASSAY_TOKEN_DROP = {
    "cla",
    "direct",
    "n.a.",
    "na",
    "n/a",
    "-",
    "other",
    "unknown",
    "flow cytometric analysis",
}

# Canonical verb -> evidence class. Verbs absent here fall through to ASSAY_UNKNOWN.
_ASSAY_TOKEN_CATEGORY = {
    # multimer / multimer-sort binding
    "tetramer-sort": ASSAY_MULTIMER,
    "dextramer-sort": ASSAY_MULTIMER,
    "pentamer-sort": ASSAY_MULTIMER,
    "streptamer-sort": ASSAY_MULTIMER,
    "multimer-sort": ASSAY_MULTIMER,
    "monomer-sort": ASSAY_MULTIMER,
    "pelimer-sort": ASSAY_MULTIMER,
    "cd8null-tetramer-sort": ASSAY_MULTIMER,
    "tetramer-magnetic-selection": ASSAY_MULTIMER,
    "tetramer-stain": ASSAY_MULTIMER,
    "multimer-stain": ASSAY_MULTIMER,
    "pentamer-stain": ASSAY_MULTIMER,
    # pooled genetic / display screens
    "phage display": ASSAY_SCREEN,
    "magnetic beads": ASSAY_SCREEN,
    "beads": ASSAY_SCREEN,
    "mhc-peptide-beads": ASSAY_SCREEN,
    "t-scan": ASSAY_SCREEN,
    "yamtad system": ASSAY_SCREEN,
    # antigen-stimulation functional readouts
    "antigen-loaded-targets": ASSAY_FUNCTIONAL,
    "antigen-expressing-targets": ASSAY_FUNCTIONAL,
    "antigen-stimulation": ASSAY_FUNCTIONAL,
    "peptide-stimulation": ASSAY_FUNCTIONAL,
    "peptide-restimulation": ASSAY_FUNCTIONAL,
    "cultured-t-cells": ASSAY_FUNCTIONAL,
    "ctl culture": ASSAY_FUNCTIONAL,
    "cd137 expression": ASSAY_FUNCTIONAL,
    "ifng capture assay": ASSAY_FUNCTIONAL,
    "nfat-luc2 luminescence": ASSAY_FUNCTIONAL,
    "nfat luminescence": ASSAY_FUNCTIONAL,
    "tnf secretion": ASSAY_FUNCTIONAL,
    "elisa": ASSAY_FUNCTIONAL,
    # clonal isolation / expansion
    "limiting-dilution-cloning": ASSAY_CLONAL,
    "enrichment": ASSAY_CLONAL,
    "ctl clone": ASSAY_CLONAL,
    # biophysical affinity
    "spr": ASSAY_AFFINITY,
    "focal molography": ASSAY_AFFINITY,
    "bli": ASSAY_AFFINITY,
    # structural
    "structural": ASSAY_STRUCTURAL,
    "x-ray crystallography": ASSAY_STRUCTURAL,
    "cryo-em": ASSAY_STRUCTURAL,
}

# IEDB / CEDAR: the query-api `assay_names` field is "readout|technique[|unit]".
# We key on the technique (2nd component); a couple of readouts disambiguate the
# generic "binding assay" technique.
_IEDB_TECHNIQUE_CATEGORY = {
    "multimer/tetramer": ASSAY_MULTIMER,
    "surface plasmon resonance (spr)": ASSAY_AFFINITY,
    "isothermal titration calorimetry": ASSAY_AFFINITY,
    "x-ray crystallography": ASSAY_STRUCTURAL,
    "cryo-electron microscopy": ASSAY_STRUCTURAL,
    "nmr": ASSAY_STRUCTURAL,
    "high throughput multiplexed assay": ASSAY_SCREEN,
    "elisa": ASSAY_FUNCTIONAL,
    "ics": ASSAY_FUNCTIONAL,
    "elispot": ASSAY_FUNCTIONAL,
    "cytometric bead array": ASSAY_FUNCTIONAL,
    "3h-thymidine": ASSAY_FUNCTIONAL,
    "cfse": ASSAY_FUNCTIONAL,
    "brdu": ASSAY_FUNCTIONAL,
    "51 chromium": ASSAY_FUNCTIONAL,
    "bioassay": ASSAY_FUNCTIONAL,
    "biological activity": ASSAY_FUNCTIONAL,
    "reporter gene assay": ASSAY_FUNCTIONAL,
    "in vitro assay": ASSAY_FUNCTIONAL,
    "in vivo assay": ASSAY_FUNCTIONAL,
    "in vivo skin test": ASSAY_FUNCTIONAL,
    "intracellular staining": ASSAY_FUNCTIONAL,
    "rna/dna detection": ASSAY_FUNCTIONAL,
}
_IEDB_READOUT_CATEGORY = {
    "dissociation constant kd": ASSAY_AFFINITY,
    "association constant ka": ASSAY_AFFINITY,
    "on rate": ASSAY_AFFINITY,
    "off rate": ASSAY_AFFINITY,
    "half life": ASSAY_AFFINITY,
}

_ASSAY_SEP = "|"

# Method string -> evidence class, filled by `combine_assays` as each source's table is
# harmonized. The two harmonized columns are independently "|"-joined *sets*, so a record
# backed by several assays no longer says which of its methods yielded which of its
# categories -- and the knowledge graph needs exactly that pairing, to link each assay
# method node to its own category.
_METHOD_CATEGORY: dict[str, str] = {}


def assay_method_categories(value) -> dict[str, str]:
    """Split one ``assay_method_raw`` cell into ``{method: its own evidence class}``.

    A method that was never categorized (or that no rule matched) maps to ``unknown``.
    """
    if value is None or str(value).strip().lower() in ("", "nan"):
        return {}
    methods = [m.strip() for m in str(value).split(_ASSAY_SEP) if m.strip()]
    return {m: _METHOD_CATEGORY.get(m.lower(), ASSAY_UNKNOWN) for m in methods}


# A concrete method string that no category rule matched falls through to
# ``unknown``; every such string is recorded here and dumped to a log file at the
# end of a run so new experiment types get noticed and mapped (mirrors the
# tidytcells `_tt_warnings` mechanism in utils.py).
_assay_warnings: set[str] = set()
_ASSAY_LOG_PATH = os.path.join("biocypher-log", f"assay_warnings_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")


def note_unmapped_assay(source: str, raw: str) -> None:
    """Record a method string that was recognized as concrete but had no category."""
    if raw and str(raw).strip():
        _assay_warnings.add(f"{source:<8} | {str(raw).strip().lower()}")


def flush_assay_warnings() -> None:
    """Write every method string that fell through to ``unknown`` to a log file."""
    if not _assay_warnings:
        return
    os.makedirs(os.path.dirname(_ASSAY_LOG_PATH), exist_ok=True)
    with open(_ASSAY_LOG_PATH, "w", encoding="utf-8") as fh:
        fh.write("# <source> | <method string> pairs that were categorized as 'unknown'.\n")
        fh.write("# Map each in the relevant table (mapping_utils.py, or the source's adapter).\n")
        for msg in sorted(_assay_warnings):
            fh.write(msg + "\n")


def _split_assay_values(values) -> list[str]:
    """Flatten mixed input (str, list, "|"/","/";"-joined str, None) into tokens."""
    if values is None:
        return []
    if isinstance(values, str):
        values = [values]
    tokens: list[str] = []
    for value in values:
        if value is None:
            continue
        for part in re.split(r"[|,;]", str(value)):
            part = part.strip()
            if part:
                tokens.append(part)
    return tokens


def canonicalize_assay_token(token: str) -> str | None:
    """Canonicalize one raw method verb (lower-cased); None if uninformative."""
    key = token.strip().lower()
    if not key or key in _ASSAY_TOKEN_DROP:
        return None
    return _ASSAY_TOKEN_CANONICAL.get(key, key)


def _join_unique(tokens) -> str | None:
    seen: list[str] = []
    for token in tokens:
        if token and token not in seen:
            seen.append(token)
    return _ASSAY_SEP.join(sorted(seen)) if seen else None


def _finalize_categories(categories) -> str:
    cats = {c for c in categories if c}
    cats.discard(ASSAY_UNKNOWN)
    if not cats:
        return ASSAY_UNKNOWN
    return _ASSAY_SEP.join(sorted(cats))


def combine_assays(pairs) -> tuple[str | None, str]:
    """Collapse ``(raw_token, category)`` pairs into the two harmonized columns.

    Shared primitive for adapters whose method vocabulary lives in the adapter
    itself (MCPAS codes, BATCAVE assays). Raw tokens are lower-cased, sorted,
    de-duplicated and ``"|"``-joined; categories likewise, with ``"unknown"``
    dropped whenever a concrete category is present. Each pair is also registered in
    ``_METHOD_CATEGORY``, so :func:`assay_method_categories` can recover which method
    yielded which category after both have been collapsed into joined sets.

    Returns:
        ``(assay_method_raw, assay_category)`` — the category is never None.
    """
    pairs = list(pairs)
    for raw_token, category in pairs:
        # A concrete category always wins over `unknown`, so the order in which sources
        # are read cannot change what `assay_method_categories` reports.
        key = str(raw_token).strip().lower() if raw_token else ""
        if key and _METHOD_CATEGORY.get(key, ASSAY_UNKNOWN) == ASSAY_UNKNOWN:
            _METHOD_CATEGORY[key] = category
    raw = _join_unique(str(r).strip().lower() for r, _ in pairs if r and str(r).strip())
    return raw, _finalize_categories(c for _, c in pairs)


def _verb_category(token: str, source: str) -> str:
    """Category for one canonical verb, recording it when nothing matches."""
    category = _ASSAY_TOKEN_CATEGORY.get(token)
    if category is None:
        note_unmapped_assay(source, token)
        return ASSAY_UNKNOWN
    return category


def harmonize_assay_verbs(values, source: str = "?") -> tuple[str | None, str]:
    """Harmonize free-text identification verbs (VDJDB / TRAIT vocabulary).

    Returns:
        (assay_method_raw, assay_category), the latter never None (``"unknown"``).
    """
    tokens = [t for t in (canonicalize_assay_token(t) for t in _split_assay_values(values)) if t]
    return combine_assays((t, _verb_category(t, source)) for t in tokens)


def harmonize_assay_iedb(assay_names, source: str = "IEDB") -> tuple[str | None, str]:
    """Harmonize IEDB/CEDAR ``assay_names`` strings ("readout|technique|unit").

    The method string is normally the technique alone. A *generic* technique (one that
    doesn't decide the category by itself, notably ``"binding assay"``) is qualified with
    its readout -- ``"binding assay (dissociation constant kd)"`` -- so the affinity-bearing
    and uninformative senses stay distinct. They otherwise collapse onto one method string,
    and since a method carries exactly one category through the knowledge graph, the whole
    node would inherit whichever sense happened to be concrete.
    """
    pairs: list[tuple[str, str]] = []
    for name in _iter_nonempty(assay_names):
        parts = [p.strip() for p in str(name).split("|") if p.strip()]
        if not parts:
            continue
        readout = parts[0].lower()
        technique = parts[1].lower() if len(parts) > 1 else ""
        raw = technique or readout
        category = _IEDB_TECHNIQUE_CATEGORY.get(technique)
        if category is None or technique == "binding assay":
            readout_category = _IEDB_READOUT_CATEGORY.get(readout)
            if readout_category is not None:
                category = readout_category
                # The technique alone didn't decide it, so keep the deciding readout in the
                # method string. Without a technique, `raw` is already the readout.
                if technique:
                    raw = f"{technique} ({readout})"
        if not category:
            note_unmapped_assay(source, raw)
            category = ASSAY_UNKNOWN
        pairs.append((raw, category))
    return combine_assays(pairs)


def _iter_nonempty(values):
    if values is None:
        return
    if isinstance(values, str):
        values = [values]
    for value in values:
        if value is not None and str(value).strip():
            yield value


def map_antigen_names(antigen_list: list[str]) -> list[str]:
    """Clean antigen names by removing bracketed species/organism info
    Args:
        antigen_list: List of antigen names to clean.

    Returns:
        Dictionary mapping original names to cleaned names.
    """
    # TODO: improve antigen names harmonization
    cleaned_map = {}
    for name in antigen_list:
        if not name:
            continue
        original = str(name).strip()

        # Remove bracketed species/organism/etc. info
        cleaned = re.sub(r"\[.*?\]", "", original)

        # Normalize whitespace
        cleaned = " ".join(cleaned.strip().split())

        cleaned_map[original] = cleaned

    return cleaned_map
