"""Guard: the Zenodo concept DOI must agree across all metadata files."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).parents[1]
CONCEPT_DOI = "10.5281/zenodo.24043645"


def test_citation_cff_has_concept_doi():
    text = (ROOT / "CITATION.cff").read_text()
    assert re.search(rf"^doi:\s*{re.escape(CONCEPT_DOI)}\s*$", text, re.M)


def test_zenodo_json_cites_concept_doi_and_no_placeholder():
    raw = (ROOT / ".zenodo.json").read_text()
    json.loads(raw)
    assert CONCEPT_DOI in raw
    assert "CONCEPT_DOI_PENDING" not in raw


def test_readme_badge_and_bibtex_use_concept_doi():
    text = (ROOT / "README.md").read_text()
    assert "DOI_PENDING" not in text
    assert text.count(CONCEPT_DOI) >= 3
