"""DNA-modification chemistry index behaviour tests.

Covers:
- Direct ``DNA_MOD_INDEX.lookup()`` for canonical symbols, aliases,
  case folding, and notational variants like
  ``5-(hydroxymethyl)cytosine``.
- CURIE provenance — entries with a verified upstream emit a ``CHEBI:``
  CURIE and carry ``verified == "chebi"``; entries without one fall back
  to ``EpiMod:`` and carry ``verified == "local"``.
- End-to-end ``enrich_claim`` resolution under the methylation
  qualifiers (``methylated_in`` on bulk, ``methylated_in_cell_type`` on
  single-cell, ``differentially_methylated_at_probe`` on array) —
  verifying the resolver wiring + methylation qualifier preferences.
- Disjointness — histone-mark lookups no longer return DNA-mod rows
  (a regression guard for the split out of ``_histone_marks.py``).
- Targets retag — the DNA-methylation machinery enzymes (TETs, MBDs)
  carry their corrected functional class.
"""
from __future__ import annotations

import re

import pytest

from indicium_adapters_epigenomics import (
    BulkEpigenomicsAdapter,
    MethylationArrayEpigenomicsAdapter,
    SingleCellEpigenomicsAdapter,
)
from indicium_adapters_epigenomics._dna_mods import (
    DNA_MOD_INDEX,
    DNAModIndex,
    _load_dna_mod_index_rows,
)
from indicium_adapters_epigenomics._enrichment import (
    HISTONE_INDEX,
    TARGET_INDEX,
)

# --- Index integrity --------------------------------------------------------


def test_dna_mod_index_is_loaded_and_nonempty():
    assert isinstance(DNA_MOD_INDEX, DNAModIndex)
    assert len(DNA_MOD_INDEX) == 5  # v0.2.0: 5mC / 5hmC / 5fC / 5caC / 6mA


def test_every_row_has_known_verified_value():
    rows = _load_dna_mod_index_rows()
    offenders = [
        row["symbol"]
        for row in rows
        if row["verified"] not in {"chebi", "local"}
    ]
    assert not offenders, f"unknown verified values: {offenders}"


def test_verified_chebi_rows_carry_chebi_prefix():
    """The verified field is a load-bearing provenance discriminator —
    a ``verified == "chebi"`` row that doesn't carry a ``CHEBI:`` CURIE
    would lie to downstream consumers about the upstream authority.
    """
    rows = _load_dna_mod_index_rows()
    mismatched = [
        (row["symbol"], row["curie"])
        for row in rows
        if row["verified"] == "chebi" and not row["curie"].startswith("CHEBI:")
    ]
    assert not mismatched, (
        f"verified=chebi rows missing CHEBI: prefix: {mismatched}"
    )


def test_local_rows_carry_epimod_prefix():
    rows = _load_dna_mod_index_rows()
    mismatched = [
        (row["symbol"], row["curie"])
        for row in rows
        if row["verified"] == "local"
        and not row["curie"].startswith("EpiMod:")
    ]
    assert not mismatched, (
        f"verified=local rows missing EpiMod: prefix: {mismatched}"
    )


def test_every_type_is_snake_case():
    snake = re.compile(r"^[a-z][a-z0-9_]*$")
    rows = _load_dna_mod_index_rows()
    offenders = [row["symbol"] for row in rows if not snake.match(row["type"])]
    assert not offenders, f"non-snake_case DNA-mod types: {offenders}"


# --- Lookup variants --------------------------------------------------------


@pytest.mark.parametrize(
    "text, expected_curie",
    [
        # Canonical symbols.
        ("5mC", "CHEBI:27551"),
        ("5hmC", "CHEBI:76792"),
        ("5fC", "CHEBI:76794"),
        ("5caC", "EpiMod:5caC"),
        ("6mA", "EpiMod:6mA"),
        # Case-insensitive.
        ("5MC", "CHEBI:27551"),
        ("5hmc", "CHEBI:76792"),
        # Dash / paren variants.
        ("5-mC", "CHEBI:27551"),
        ("5-hmC", "CHEBI:76792"),
        ("5-fC", "CHEBI:76794"),
        # Long-form aliases.
        ("5-methylcytosine", "CHEBI:27551"),
        ("5-hydroxymethylcytosine", "CHEBI:76792"),
        ("5-(hydroxymethyl)cytosine", "CHEBI:76792"),
        ("5-formylcytosine", "CHEBI:76794"),
        ("5-carboxylcytosine", "EpiMod:5caC"),
        ("N6-methyladenine", "EpiMod:6mA"),
    ],
)
def test_dna_mod_lookup_resolves_known_variants(text, expected_curie):
    match = DNA_MOD_INDEX.lookup(text)
    assert match is not None, f"{text!r} should resolve"
    assert match["curie"] == expected_curie


def test_dna_mod_lookup_returns_none_for_unknown():
    assert DNA_MOD_INDEX.lookup("not-a-dna-modification") is None


def test_dna_mod_lookup_returns_none_for_non_string():
    assert DNA_MOD_INDEX.lookup(None) is None  # type: ignore[arg-type]


def test_dna_mod_lookup_payload_carries_verified_field():
    chebi_match = DNA_MOD_INDEX.lookup("5mC")
    assert chebi_match is not None
    assert chebi_match["verified"] == "chebi"
    local_match = DNA_MOD_INDEX.lookup("5caC")
    assert local_match is not None
    assert local_match["verified"] == "local"


# --- Regression guard: histone index no longer returns DNA mods -------------


@pytest.mark.parametrize("dna_mod_symbol", ["5mC", "5hmC", "5fC", "5caC", "6mA"])
def test_histone_index_no_longer_returns_dna_mods(dna_mod_symbol):
    """In v0.1.0 these rows lived in ``HISTONE_INDEX``. v0.2.0 §3 split
    them out; this is the regression guard that catches a future
    accidental re-introduction.
    """
    assert HISTONE_INDEX.lookup(dna_mod_symbol) is None


# --- End-to-end: methylation qualifiers route through DNA_MOD --------------


def test_bulk_methylated_in_resolves_to_dna_mod_curie():
    bulk = BulkEpigenomicsAdapter()
    enriched = bulk.enrich_claim(
        {
            "id": "c-mod-1",
            "subject": "TAP2_promoter",
            "qualifier": "methylated_in",
            "object": "5mC",
        }
    )
    assert enriched["ontology_terms"]["object"] == "CHEBI:27551"
    assert enriched["xrefs"]["object"]["source"] == "dna_mod"
    assert enriched["xrefs"]["object"]["verified"] == "chebi"


def test_singlecell_methylated_in_cell_type_resolves_to_dna_mod_curie():
    sc = SingleCellEpigenomicsAdapter()
    enriched = sc.enrich_claim(
        {
            "id": "c-mod-2",
            "subject": "CD4+ T cell",
            "qualifier": "methylated_in_cell_type",
            "object": "5hmC",
        }
    )
    assert enriched["ontology_terms"]["object"] == "CHEBI:76792"
    assert enriched["xrefs"]["object"]["source"] == "dna_mod"


def test_array_differentially_methylated_at_probe_resolves_to_dna_mod_curie():
    arr = MethylationArrayEpigenomicsAdapter()
    enriched = arr.enrich_claim(
        {
            "id": "c-mod-3",
            "subject": "cg00000029",
            "qualifier": "differentially_methylated_at_probe",
            "object": "5mC",
        }
    )
    assert enriched["ontology_terms"]["object"] == "CHEBI:27551"
    assert enriched["xrefs"]["object"]["source"] == "dna_mod"


def test_methylation_qualifier_falls_through_to_target_for_enzyme_symbol():
    """The methylation preference is ``(DNA_MOD, TARGET)`` — TARGET as
    the second hop so enzyme symbols like ``DNMT3A`` (a target, not a
    chemistry) still resolve when used as the object of a methylation
    qualifier.
    """
    bulk = BulkEpigenomicsAdapter()
    enriched = bulk.enrich_claim(
        {
            "id": "c-mod-4",
            "subject": "TAP2_promoter",
            "qualifier": "methylated_in",
            "object": "DNMT3A",
        }
    )
    assert enriched["ontology_terms"]["object"] == "UniProt:Q9Y6K1"
    assert enriched["xrefs"]["object"]["source"] == "target"


# --- TET / MBD retag — corrected functional class --------------------------


@pytest.mark.parametrize(
    "symbol, expected_category",
    [
        ("TET1", "dna_demethylase"),
        ("TET2", "dna_demethylase"),
        ("TET3", "dna_demethylase"),
        ("MBD1", "methyl_cpg_reader"),
        ("MBD2", "methyl_cpg_reader"),
        ("MBD3", "methyl_cpg_reader"),
        ("MBD4", "methyl_cpg_reader"),
    ],
)
def test_methylation_machinery_enzymes_have_corrected_category(
    symbol, expected_category
):
    """v0.1.0 tagged all TETs and MBDs as ``chromatin_remodeler`` which
    was domain-incorrect (TETs are demethylases, MBDs are readers).
    v0.2.0 §3 retags them to their actual functional class for downstream
    machine-reading clarity. Behaviour-pure: the resolver only uses the
    ``curie`` field of a TARGET row, so this is metadata cleanup, but
    it's an observable contract on the ``xrefs`` payload that ASB
    consumers may rely on.
    """
    row = TARGET_INDEX.lookup(symbol)
    assert row is not None
    assert row["category"] == expected_category
