import sys

sys.path.append("..")

import re
from urllib.parse import quote

import requests


def map_species_terms(terms: list[str]) -> dict:
    manual_disambiguation = {
        "EBV": "human gammaherpesvirus 4",
        "CMV": "Cytomegalovirus",
        "HIV": "Human immunodeficiency virus",
        "HSV": "Herpes simplex virus",
        "HCV": "Hepatitis C virus",
        "HHV": "Human herpesvirus",
        "DENV": "Dengue virus",
        "YFV": "Yellow fever virus",
        "HPV": "Human papillomavirus",
        "Mtb": "Mycobacterium tuberculosis",
        "SIV": "Simian immunodeficiency virus",
        "AdV": "Human adenovirus",
        "MCPyV": "Merkel cell polyomavirus",
        "McpyV": "Merkel cell polyomavirus",
        "M. tuberculosis": "Mycobacterium tuberculosis",
        "SARS-CoV2": "Severe acute respiratory syndrome coronavirus 2",
        "InfluenzaA": "Influenza A virus",
        "H5N1 subtype": "Influenza A virus",
    }

    def normalize_species(term: str) -> str:
        term = term.strip()
        term = re.sub(r"^([a-zA-Z]+)(\d+)(?![a-zA-Z])", r"\1 \2", term)
        for prefix in manual_disambiguation:
            if term.startswith(prefix):
                suffix = term[len(prefix) :]
                query_term = manual_disambiguation[prefix] + suffix
                break
        else:
            query_term = term
        query_term = query_term.replace("_", " ")
        query_term = re.sub(r"\s*[\(\[].*[\)\]]", "", query_term)
        query_term = re.sub(r"\bstrain\s.*", "", query_term).strip()

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
        try:
            if "obo/" not in uri:
                return None, None
            term = uri.split("obo/")[-1]
            ontology = term.split("_")[0].lower()
            full_uri = f"http://purl.obolibrary.org/obo/{term}"
            encoded_uri = quote(quote(full_uri, safe=""), safe="")
            ols_url = f"https://www.ebi.ac.uk/ols4/api/ontologies/{ontology}/terms/{encoded_uri}"
            res = requests.get(ols_url, timeout=10)
            res.raise_for_status()
            label = res.json().get("label")
            return label, full_uri
        except:
            return None, None

    def get_zooma_label(term: str):
        if term.startswith("http://purl.obolibrary.org/"):
            label, iri = get_label_from_semantic_tag(term)
            return label
        zooma_url = "https://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate"
        sources = ["gwas", "uniprot", "ebisc"]
        params = {
            "propertyValue": term,
            "propertyType": "organism",
            "ontologies": "ncbitaxon",
            "filter": f"required:[{','.join(sources)}]",
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
    normalized_terms = {term: normalize_species(term) for term in terms}
    # print("Normalized terms:", normalized_terms)

    # Step 2: Get Zooma mappings for normalized terms
    # zooma_output = {}
    # for original_term, normalized_term in tqdm(normalized_terms.items(), desc="Getting Zooma mappings"):
    #     zooma_result = get_zooma_label(normalized_term)
    #     zooma_output[original_term] = zooma_result

    # # print("Zooma output:", zooma_output)

    # # Step 3: Create final results - use Zooma output if available, otherwise use normalized term
    # results = {}
    # for original_term in terms:
    #     zooma_result = zooma_output[original_term]
    #     if zooma_result is not None:
    #         results[original_term] = zooma_result
    #     else:
    #         results[original_term] = normalized_terms[original_term]

    return normalized_terms
