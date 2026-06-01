"""IHEC core biosample index behaviour tests.

Covers:
- Direct ``IHEC_BIOSAMPLE_INDEX.lookup()`` for canonical symbols,
  aliases, case folding, and the broader dash / plus / underscore /
  whitespace normalisation needed for biosample literature (``K-562``,
  ``CD4+ T cell``, ``HeLa S3``).
- Cross-reference payload — populated entries return ``cl_curie``,
  ``uberon_curie`` and / or ``efo_curie`` in the lookup dict; entries
  without a stable upstream cross-ref intentionally omit those keys
  rather than emit empty strings.
- End-to-end ``enrich_claim`` resolution under accessibility / cell-type
  qualifiers (``accessible_in``, ``accessible_in_cell_type``,
  ``differentially_accessible_in``, ``marker_peak_of_cell_type``) —
  verifying the resolver wiring + biosample-first qualifier preferences.
- Fall-through behaviour — methylated_in is now ``(DNA_MOD, TARGET,
  IHEC_BIOSAMPLE)``: a chemistry object still wins DNA_MOD; an enzyme
  wins TARGET; a cell-type name correctly falls through to IHEC.
- CURIE provenance — every primary CURIE is ``IHECBio:<localpart>`` with
  no characters that break CURIE round-tripping (no dashes, no plus
  signs, no spaces).
- Vocabulary structure — known ``type`` categories, snake_case.
"""
from __future__ import annotations

import re

import pytest

from indicium_adapters_epigenomics import (
    BulkEpigenomicsAdapter,
    SingleCellEpigenomicsAdapter,
)
from indicium_adapters_epigenomics._biosamples import (
    IHEC_BIOSAMPLE_INDEX,
    IHECBiosampleIndex,
    _load_ihec_biosample_index_rows,
)

# --- Index integrity --------------------------------------------------------


def test_ihec_biosample_index_is_loaded_and_nonempty():
    assert isinstance(IHEC_BIOSAMPLE_INDEX, IHECBiosampleIndex)
    assert len(IHEC_BIOSAMPLE_INDEX) >= 30  # v0.2.0 ships ~33 curated rows


def test_every_row_has_ihecbio_curie_with_clean_localpart():
    """``IHECBio:<localpart>`` — local part must not contain whitespace,
    dash, plus, or other characters that break CURIE round-tripping.
    """
    snake_compat = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
    rows = _load_ihec_biosample_index_rows()
    bad: list[tuple[str, str]] = []
    for row in rows:
        curie = row["curie"]
        if not curie.startswith("IHECBio:"):
            bad.append((row["symbol"], curie))
            continue
        localpart = curie.split(":", 1)[1]
        if not snake_compat.match(localpart):
            bad.append((row["symbol"], curie))
    assert not bad, f"malformed IHECBio CURIE local parts: {bad}"


def test_every_row_type_is_a_known_biosample_category():
    """``type`` distinguishes cell_line / primary_cell / stem_cell /
    tissue — used downstream by ASB to decide which UBERON-vs-CL
    cross-ref to surface in display.
    """
    known = {"cell_line", "primary_cell", "stem_cell", "tissue"}
    offenders = [
        (row["symbol"], row["type"])
        for row in _load_ihec_biosample_index_rows()
        if row["type"] not in known
    ]
    assert not offenders, f"unknown biosample types: {offenders}"


def test_every_cl_curie_has_obo_form():
    """Cross-ref columns must hold valid OBO-style CURIEs when populated.
    Empty cells are allowed (intentionally — see _biosamples.py docstring).
    """
    cl_pat = re.compile(r"^CL:\d{7}$")
    uberon_pat = re.compile(r"^UBERON:\d{7}$")
    efo_pat = re.compile(r"^EFO:\d{7}$")
    bad: list[tuple[str, str, str]] = []
    for row in _load_ihec_biosample_index_rows():
        for col, pat in (
            ("cl_curie", cl_pat),
            ("uberon_curie", uberon_pat),
            ("efo_curie", efo_pat),
        ):
            value = row.get(col, "")
            if value and not pat.match(value):
                bad.append((row["symbol"], col, value))
    assert not bad, f"malformed cross-ref CURIEs: {bad}"


# --- Lookup variants --------------------------------------------------------


@pytest.mark.parametrize(
    "text, expected_curie",
    [
        # Canonical symbols.
        ("K562", "IHECBio:K562"),
        ("GM12878", "IHECBio:GM12878"),
        ("HepG2", "IHECBio:HepG2"),
        # Case-insensitive.
        ("k562", "IHECBio:K562"),
        ("GM12878", "IHECBio:GM12878"),
        # Dash / plus / underscore / whitespace normalisation.
        ("K-562", "IHECBio:K562"),
        ("K 562", "IHECBio:K562"),
        ("HeLa-S3", "IHECBio:HeLa_S3"),
        ("HeLa S3", "IHECBio:HeLa_S3"),
        ("HeLaS3", "IHECBio:HeLa_S3"),
        ("CD4+ T cell", "IHECBio:CD4_T_cell"),
        ("CD4-positive T cell", "IHECBio:CD4_T_cell"),
        ("CD8+ T cell", "IHECBio:CD8_T_cell"),
        # Long-form aliases.
        ("Human Umbilical Vein Endothelial Cells", "IHECBio:HUVEC"),
        ("hematopoietic stem cell", None),  # not in our alias set
    ],
)
def test_ihec_biosample_lookup_resolves_known_variants(text, expected_curie):
    match = IHEC_BIOSAMPLE_INDEX.lookup(text)
    if expected_curie is None:
        assert match is None
    else:
        assert match is not None, f"{text!r} should resolve"
        assert match["curie"] == expected_curie


def test_ihec_biosample_lookup_returns_none_for_unknown():
    assert IHEC_BIOSAMPLE_INDEX.lookup("not-a-biosample") is None


def test_ihec_biosample_lookup_returns_none_for_non_string():
    assert IHEC_BIOSAMPLE_INDEX.lookup(None) is None  # type: ignore[arg-type]


# --- Cross-reference payload ------------------------------------------------


def test_lookup_payload_includes_populated_cross_refs():
    """GM12878 has CL + EFO cross-refs but no UBERON — uberon_curie must
    be absent from the returned dict, not present-but-empty.
    """
    match = IHEC_BIOSAMPLE_INDEX.lookup("GM12878")
    assert match is not None
    assert match["cl_curie"] == "CL:0000236"
    assert match["efo_curie"] == "EFO:0002824"
    assert "uberon_curie" not in match


def test_lookup_payload_omits_empty_cross_refs():
    """K562 intentionally has no CL cross-ref — the K562 lineage in CL
    is contested. Empty cells must be omitted, not surfaced as empty
    strings.
    """
    match = IHEC_BIOSAMPLE_INDEX.lookup("K562")
    assert match is not None
    assert "cl_curie" not in match
    assert "uberon_curie" not in match
    assert match["efo_curie"] == "EFO:0002067"


def test_tissue_entry_has_uberon_but_no_cl():
    """Tissues are anatomical sites — UBERON applies but CL does not."""
    match = IHEC_BIOSAMPLE_INDEX.lookup("liver")
    assert match is not None
    assert match["type"] == "tissue"
    assert match["uberon_curie"] == "UBERON:0002107"
    assert "cl_curie" not in match


# --- End-to-end: biosample-aware qualifier preferences ----------------------


def test_bulk_accessible_in_resolves_to_ihecbio():
    bulk = BulkEpigenomicsAdapter()
    enriched = bulk.enrich_claim(
        {
            "id": "c-bio-1",
            "subject": "ATAC_peak_42",
            "qualifier": "accessible_in",
            "object": "K562",
        }
    )
    assert enriched["ontology_terms"]["object"] == "IHECBio:K562"
    assert enriched["xrefs"]["object"]["source"] == "ihec_biosample"


def test_singlecell_accessible_in_cell_type_resolves_to_ihecbio():
    sc = SingleCellEpigenomicsAdapter()
    enriched = sc.enrich_claim(
        {
            "id": "c-bio-2",
            "subject": "ATAC_peak_99",
            "qualifier": "accessible_in_cell_type",
            "object": "CD4+ T cell",
        }
    )
    assert enriched["ontology_terms"]["object"] == "IHECBio:CD4_T_cell"
    assert enriched["xrefs"]["object"]["source"] == "ihec_biosample"
    # Cross-refs are preserved into the xrefs payload via **match.
    assert enriched["xrefs"]["object"]["cl_curie"] == "CL:0000624"


def test_singlecell_marker_peak_of_cell_type_resolves_to_ihecbio_first():
    """``marker_peak_of_cell_type`` preference is
    ``(IHEC_BIOSAMPLE, HISTONE, TARGET)``. With a biosample object, IHEC
    wins; with a histone-mark object, the second hop fires.
    """
    sc = SingleCellEpigenomicsAdapter()
    enriched_biosample = sc.enrich_claim(
        {
            "id": "c-bio-3",
            "subject": "peak_X",
            "qualifier": "marker_peak_of_cell_type",
            "object": "CD8+ T cell",
        }
    )
    assert enriched_biosample["xrefs"]["object"]["source"] == "ihec_biosample"
    enriched_histone = sc.enrich_claim(
        {
            "id": "c-bio-3b",
            "subject": "peak_X",
            "qualifier": "marker_peak_of_cell_type",
            "object": "H3K4me3",
        }
    )
    assert enriched_histone["xrefs"]["object"]["source"] == "histone"


def test_methylated_in_extends_to_ihecbio_for_cell_type_object():
    """v0.2.0 §4 extends the methylation preference to ``(DNA_MOD,
    TARGET, IHEC_BIOSAMPLE)`` so a cell-type object now resolves too.
    Chemistry-object and enzyme-object regressions from §3 still hold.
    """
    bulk = BulkEpigenomicsAdapter()
    # Chemistry → DNA_MOD.
    chem = bulk.enrich_claim(
        {"id": "m1", "subject": "site", "qualifier": "methylated_in", "object": "5mC"}
    )
    assert chem["xrefs"]["object"]["source"] == "dna_mod"
    # Enzyme → TARGET (DNMT3A).
    enz = bulk.enrich_claim(
        {"id": "m2", "subject": "site", "qualifier": "methylated_in", "object": "DNMT3A"}
    )
    assert enz["xrefs"]["object"]["source"] == "target"
    # Cell type → IHEC_BIOSAMPLE (was missing before §4).
    bio = bulk.enrich_claim(
        {"id": "m3", "subject": "site", "qualifier": "methylated_in", "object": "K562"}
    )
    assert bio["xrefs"]["object"]["source"] == "ihec_biosample"
