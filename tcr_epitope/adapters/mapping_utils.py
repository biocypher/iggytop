""" This module contains utility functions for cleaning up species terms and cleaning antigen names """

import sys

sys.path.append("..")

import re
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
        "EBV": "Human gammaherpesvirus 4",
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
    }

    def normalize_species(term: str) -> str:
        """Normalize species terms by applying manual mappings for abbreviations,
        cleaning up formatting"""
        if term.startswith("http"):
            label, iri = get_label_from_semantic_tag(term)
            return label
        term = term.strip()
        term = re.sub(r"^([a-zA-Z]+)(\d+)(?![a-zA-Z])", r"\1 \2", term)
        for prefix in manual_disambiguation:
            if term.startswith(prefix):
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
        except:
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
        except:
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

    # Step 1: Normalize all terms
    normalized_terms = {term: normalize_species(term) for term in terms if term}

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


def map_antigen_names(antigen_list: list[str]) -> list[str]:
    """Clean antigen names by removing bracketed species/organism info"""
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
