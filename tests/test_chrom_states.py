"""ChromHMM / Segway chromatin-state index behaviour tests.

Covers:
- Direct ``CHROM_STATE_INDEX.lookup()`` for canonical labels, aliases,
  case folding, and the underscore / slash / apostrophe / whitespace
  normalisation.
- End-to-end ``enrich_claim`` resolution under the two new qualifiers
  (``chromatin_state_of`` on bulk, ``marker_state_of_cell_type`` on
  single-cell) — verifying the resolver wiring + the qualifier preference
  table actually route to ``CHROM_STATE_INDEX``.
- Disjointness — confirming the new qualifiers don't bleed into other
  adapter sets.
- Vocabulary structure — every row has a populated ``EpiState:`` CURIE,
  a known model membership token, and a snake_case ``type`` category.
"""
from __future__ import annotations

import re

import pytest

from indicium_adapters_epigenomics import (
    BulkEpigenomicsAdapter,
    MethylationArrayEpigenomicsAdapter,
    SingleCellEpigenomicsAdapter,
)
from indicium_adapters_epigenomics._chrom_states import (
    CHROM_STATE_INDEX,
    ChromStateIndex,
    _load_chrom_state_index_rows,
)

# --- Index integrity --------------------------------------------------------


def test_chrom_state_index_is_loaded_and_nonempty():
    assert isinstance(CHROM_STATE_INDEX, ChromStateIndex)
    assert len(CHROM_STATE_INDEX) >= 30  # we ship ~42 curated rows


def test_every_row_has_epistate_curie_with_snake_compatible_localpart():
    """``EpiState:<localpart>`` — local part must avoid characters that
    don't round-trip through CURIE expansion (slashes, apostrophes).
    """
    rows = _load_chrom_state_index_rows()
    bad: list[tuple[str, str]] = []
    for row in rows:
        curie = row["curie"]
        if not curie.startswith("EpiState:"):
            bad.append((row["symbol"], curie))
            continue
        localpart = curie.split(":", 1)[1]
        if re.search(r"[/'\s]", localpart):
            bad.append((row["symbol"], curie))
    assert not bad, f"malformed EpiState CURIE local parts: {bad}"


def test_every_row_models_field_uses_known_tokens():
    """``models`` column is pipe-separated. Every token must come from
    the known segmentation-model vocabulary; a typo would silently
    misclassify a state.
    """
    known_models = {"chromhmm_15", "chromhmm_18", "chromhmm_25", "segway"}
    offenders: list[tuple[str, str]] = []
    for row in _load_chrom_state_index_rows():
        for token in row["models"].split("|"):
            token = token.strip()
            if token and token not in known_models:
                offenders.append((row["symbol"], token))
    assert not offenders, f"unknown model tokens: {offenders}"


def test_every_row_type_is_snake_case():
    snake = re.compile(r"^[a-z][a-z0-9_]*$")
    offenders = [
        row["symbol"]
        for row in _load_chrom_state_index_rows()
        if not snake.match(row["type"])
    ]
    assert not offenders, f"non-snake_case state categories: {offenders}"


# --- Lookup variants --------------------------------------------------------


@pytest.mark.parametrize(
    "text, expected_curie",
    [
        # Canonical symbols.
        ("TssA", "EpiState:TssA"),
        ("EnhA1", "EpiState:EnhA1"),
        ("ReprPC", "EpiState:ReprPC"),
        ("Quies", "EpiState:Quies"),
        # Case-insensitive.
        ("tssa", "EpiState:TssA"),
        ("REPRPC", "EpiState:ReprPC"),
        # Aliases — Segway / underscore notation.
        ("Active_TSS", "EpiState:TssA"),
        ("ActiveTSS", "EpiState:TssA"),
        ("Active TSS", "EpiState:TssA"),
        # Slash / underscore / dash normalisation.
        ("ZNF/Rpts", "EpiState:ZNF_Rpts"),
        ("ZNF_Rpts", "EpiState:ZNF_Rpts"),
        ("ZnfRpts", "EpiState:ZNF_Rpts"),
        # Apostrophe / whitespace stripping.
        ("Tx5'", "EpiState:Tx5"),
        ("Tx 5", "EpiState:Tx5"),
        # Long-form match via aliases.
        ("Heterochromatin", "EpiState:Het"),
        ("Polycomb", "EpiState:ReprPC"),
        ("Quiescent/Low", "EpiState:Quies"),
    ],
)
def test_chrom_state_lookup_resolves_known_variants(text, expected_curie):
    match = CHROM_STATE_INDEX.lookup(text)
    assert match is not None, f"{text!r} should resolve"
    assert match["curie"] == expected_curie


def test_chrom_state_lookup_returns_none_for_unknown():
    assert CHROM_STATE_INDEX.lookup("definitely-not-a-state") is None


def test_chrom_state_lookup_returns_none_for_non_string():
    assert CHROM_STATE_INDEX.lookup(None) is None  # type: ignore[arg-type]


def test_chrom_state_lookup_payload_carries_models_field():
    match = CHROM_STATE_INDEX.lookup("TssA")
    assert match is not None
    assert "models" in match
    # TssA appears in all three ChromHMM models per the bundled TSV.
    for token in ("chromhmm_15", "chromhmm_18", "chromhmm_25"):
        assert token in match["models"]


# --- End-to-end: qualifier-aware resolution wires to CHROM_STATE_INDEX -----


def test_bulk_chromatin_state_of_resolves_to_epistate():
    """``chromatin_state_of`` is registered as a CHROM_STATE-preferring
    qualifier in ``_resolver._QUALIFIER_PREFERENCES``. End-to-end through
    ``enrich_claim``, a state-label object must come back as an
    ``EpiState:`` CURIE.
    """
    bulk = BulkEpigenomicsAdapter()
    enriched = bulk.enrich_claim(
        {
            "id": "c-state-1",
            "subject": "K562_peak_42",
            "qualifier": "chromatin_state_of",
            "object": "TssA",
        }
    )
    assert enriched["ontology_terms"]["object"] == "EpiState:TssA"
    assert enriched["xrefs"]["object"]["source"] == "chrom_state"
    # The full ChromHMM-model trace is preserved into xrefs via the
    # **match splat.
    assert "models" in enriched["xrefs"]["object"]


def test_singlecell_marker_state_of_cell_type_resolves_to_epistate():
    sc = SingleCellEpigenomicsAdapter()
    enriched = sc.enrich_claim(
        {
            "id": "c-state-2",
            "subject": "CD8+ T cell",
            "qualifier": "marker_state_of_cell_type",
            "object": "EnhA1",
        }
    )
    assert enriched["ontology_terms"]["object"] == "EpiState:EnhA1"
    assert enriched["xrefs"]["object"]["source"] == "chrom_state"


def test_bulk_chromatin_state_of_falls_through_for_unknown_state_label():
    """The qualifier preference is a hint, not a route: an unknown label
    under ``chromatin_state_of`` must still fall through to the default
    cascade rather than wedging at CHROM_STATE.
    """
    bulk = BulkEpigenomicsAdapter()
    enriched = bulk.enrich_claim(
        {
            "id": "c-state-3",
            "subject": "K562_peak_99",
            "qualifier": "chromatin_state_of",
            "object": "not-a-known-state-label-XYZ",
        }
    )
    assert "object" not in enriched["ontology_terms"]


# --- Disjointness invariant -------------------------------------------------


def test_state_qualifiers_are_disjoint_across_sibling_adapters():
    """Each new state qualifier belongs to exactly one adapter."""
    bulk_qs = BulkEpigenomicsAdapter.qualifiers
    sc_qs = SingleCellEpigenomicsAdapter.qualifiers
    arr_qs = MethylationArrayEpigenomicsAdapter.qualifiers

    assert "chromatin_state_of" in bulk_qs
    assert "chromatin_state_of" not in sc_qs
    assert "chromatin_state_of" not in arr_qs

    assert "marker_state_of_cell_type" in sc_qs
    assert "marker_state_of_cell_type" not in bulk_qs
    assert "marker_state_of_cell_type" not in arr_qs
