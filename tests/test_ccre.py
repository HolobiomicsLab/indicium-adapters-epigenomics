"""ENCODE SCREEN cCRE classification index + accession recogniser tests.

Covers:
- Direct ``CCRE_CLASS_INDEX.lookup()`` for the 7 SCREEN class labels
  across notational variants (``CA-CTCF`` / ``CA_CTCF`` / ``CA+CTCF`` /
  ``CA CTCF``, ``Promoter-Like Signature`` → PLS, etc.).
- ``EH38E\\d{7}`` (hg38) and ``EM10E\\d{7}`` (mm10) cCRE-accession
  regex routing through ``accession_to_curie`` to ``ENCODE:`` with the
  ``ENCODE SCREEN`` repository label.
- End-to-end resolution: ``classified_as_ccre`` qualifier routes a
  class-label object to ``EpiCCRE:``; an ``EH38E0000001`` subject
  resolves to ``ENCODE:EH38E0000001`` regardless of qualifier (accession
  outranks qualifier preference).
- Disjointness — the new qualifier belongs to bulk only; sibling
  qualifier sets are untouched.
"""
from __future__ import annotations

import re

import pytest

from indicium_adapters_epigenomics import (
    BulkEpigenomicsAdapter,
    MethylationArrayEpigenomicsAdapter,
    SingleCellEpigenomicsAdapter,
)
from indicium_adapters_epigenomics._ccre import (
    CCRE_CLASS_INDEX,
    CCREClassIndex,
    _load_ccre_class_index_rows,
)
from indicium_adapters_epigenomics._registries import accession_to_curie

# --- Index integrity --------------------------------------------------------


def test_ccre_class_index_loads_seven_canonical_classes():
    assert isinstance(CCRE_CLASS_INDEX, CCREClassIndex)
    # ENCODE SCREEN ships exactly 7 functional classes at v0.2.0.
    assert len(CCRE_CLASS_INDEX) == 7


def test_every_row_has_epiccre_curie_with_clean_localpart():
    snake_compat = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
    rows = _load_ccre_class_index_rows()
    bad: list[tuple[str, str]] = []
    for row in rows:
        curie = row["curie"]
        if not curie.startswith("EpiCCRE:"):
            bad.append((row["symbol"], curie))
            continue
        localpart = curie.split(":", 1)[1]
        if not snake_compat.match(localpart):
            bad.append((row["symbol"], curie))
    assert not bad, f"malformed EpiCCRE CURIE local parts: {bad}"


def test_every_row_type_is_snake_case():
    snake = re.compile(r"^[a-z][a-z0-9_]*$")
    rows = _load_ccre_class_index_rows()
    offenders = [row["symbol"] for row in rows if not snake.match(row["type"])]
    assert not offenders, f"non-snake_case cCRE class types: {offenders}"


# --- Class-label lookup variants -------------------------------------------


@pytest.mark.parametrize(
    "text, expected_curie",
    [
        # Canonical symbols.
        ("PLS", "EpiCCRE:PLS"),
        ("pELS", "EpiCCRE:pELS"),
        ("dELS", "EpiCCRE:dELS"),
        ("CA-CTCF", "EpiCCRE:CA_CTCF"),
        ("CA-H3K4me3", "EpiCCRE:CA_H3K4me3"),
        ("CA-TF", "EpiCCRE:CA_TF"),
        ("CA-only", "EpiCCRE:CA_only"),
        # Case-insensitive.
        ("pls", "EpiCCRE:PLS"),
        ("ca-ctcf", "EpiCCRE:CA_CTCF"),
        # Dash / plus / underscore / whitespace normalisation.
        ("CA_CTCF", "EpiCCRE:CA_CTCF"),
        ("CA+CTCF", "EpiCCRE:CA_CTCF"),
        ("CA CTCF", "EpiCCRE:CA_CTCF"),
        ("CA only", "EpiCCRE:CA_only"),
        # Long-form aliases.
        ("Promoter-Like Signature", "EpiCCRE:PLS"),
        ("proximal Enhancer-Like Signature", "EpiCCRE:pELS"),
        ("distal Enhancer-Like Signature", "EpiCCRE:dELS"),
    ],
)
def test_ccre_class_lookup_resolves_known_variants(text, expected_curie):
    match = CCRE_CLASS_INDEX.lookup(text)
    assert match is not None, f"{text!r} should resolve"
    assert match["curie"] == expected_curie


def test_ccre_class_lookup_returns_none_for_unknown():
    assert CCRE_CLASS_INDEX.lookup("not-a-ccre-class") is None


def test_ccre_class_lookup_returns_none_for_non_string():
    assert CCRE_CLASS_INDEX.lookup(None) is None  # type: ignore[arg-type]


# --- EH38E / EM10E accession recogniser ------------------------------------


@pytest.mark.parametrize(
    "accession",
    [
        "EH38E0000001",
        "EH38E1234567",
        "EH38E9999999",
        "EM10E0000001",
        "EM10E5555555",
    ],
)
def test_ccre_accession_routes_to_encode_screen(accession):
    result = accession_to_curie(accession)
    assert result is not None
    assert result["curie"] == f"ENCODE:{accession}"
    assert result["repository"] == "ENCODE SCREEN"
    assert result["bare_id"] == accession


@pytest.mark.parametrize(
    "non_accession",
    [
        "EH38E000001",      # 6 digits — too short
        "EH38E00000001",    # 8 digits — too long
        "EH38X0000001",     # wrong type tag (X vs E)
        "EM10F0000001",     # wrong type tag (F vs E)
        "EH38E123456A",     # alphanumeric digit slot
    ],
)
def test_invalid_ccre_accessions_do_not_match(non_accession):
    """The recogniser is anchored and width-strict — anything that's not
    exactly ``[EH38E|EM10E]\\d{7}`` must fall through to the rest of the
    accession cascade (and in these cases match nothing).
    """
    assert accession_to_curie(non_accession) is None


# --- End-to-end resolution -------------------------------------------------


def test_classified_as_ccre_resolves_class_label_to_epiccre():
    bulk = BulkEpigenomicsAdapter()
    enriched = bulk.enrich_claim(
        {
            "id": "c-ccre-1",
            "subject": "EH38E0000001",
            "qualifier": "classified_as_ccre",
            "object": "dELS",
        }
    )
    # Subject (accession) → ENCODE:.
    assert enriched["ontology_terms"]["subject"] == "ENCODE:EH38E0000001"
    assert enriched["xrefs"]["subject"]["source"] == "accession"
    assert enriched["xrefs"]["subject"]["repository"] == "ENCODE SCREEN"
    # Object (class label) → EpiCCRE: via the qualifier preference.
    assert enriched["ontology_terms"]["object"] == "EpiCCRE:dELS"
    assert enriched["xrefs"]["object"]["source"] == "ccre"


def test_ccre_accession_outranks_any_qualifier_preference():
    """Accessions always win first. Even under ``classified_as_ccre``
    (which prefers CCRE), an EH38E accession in the subject still
    resolves through the accession cascade, never through CCRE_CLASS.
    """
    bulk = BulkEpigenomicsAdapter()
    enriched = bulk.enrich_claim(
        {
            "id": "c-ccre-2",
            "subject": "EH38E1234567",
            "qualifier": "classified_as_ccre",
            "object": "PLS",
        }
    )
    assert enriched["xrefs"]["subject"]["source"] == "accession"


# --- Disjointness invariant ------------------------------------------------


def test_classified_as_ccre_is_disjoint_across_sibling_adapters():
    assert "classified_as_ccre" in BulkEpigenomicsAdapter.qualifiers
    assert "classified_as_ccre" not in SingleCellEpigenomicsAdapter.qualifiers
    assert (
        "classified_as_ccre" not in MethylationArrayEpigenomicsAdapter.qualifiers
    )
