"""SingleCellEpigenomicsAdapter — enrichment + qualifier tests.

The v0.1.0 single-cell adapter reuses the shared assay / histone / target /
QC / accession indices from ``_enrichment.py``. Cell-type (CL) and anatomy
(UBERON) enrichment is deferred to v0.2.0 per the adapter docstring — that
deferral is part of the public contract and is asserted here.
"""
from __future__ import annotations

import copy

import pytest

from indicium_adapters_epigenomics import SingleCellEpigenomicsAdapter

_EXPECTED_QUALIFIERS = frozenset({
    "accessible_in_cell_type",
    "methylated_in_cell_type",
    "bound_in_cell_type",
    "marker_peak_of_cell_type",
    "marker_dmr_of_cell_type",
    "enriched_in_cluster",
    "co_accessible_with",
    "localized_in_spatial_domain",
    # v0.2.0 §2 — chromatin-state markers per cell type.
    "marker_state_of_cell_type",
})


# --- qualifiers -------------------------------------------------------------


def test_all_documented_qualifiers_are_present():
    assert SingleCellEpigenomicsAdapter.qualifiers == _EXPECTED_QUALIFIERS


@pytest.mark.parametrize("qualifier", sorted(_EXPECTED_QUALIFIERS))
def test_each_qualifier_is_accepted_by_enrich_claim(qualifier):
    """An in-vocabulary qualifier must not cause ``enrich_claim`` to raise."""
    a = SingleCellEpigenomicsAdapter()
    # Use a known-resolvable subject so the enrichment branch fires too.
    out = a.enrich_claim({
        "subject": "H3K27ac",
        "qualifier": qualifier,
        "object": "CTCF",
    })
    assert isinstance(out, dict)


# --- shared enrichment indices still wire through --------------------------


@pytest.mark.parametrize(
    "raw, expected_curie",
    [
        ("scATAC-seq", "EFO:0010891"),
        ("snATAC-seq", "OBI:0002762"),
        ("snmC-seq", "EFO:0008939"),
        ("snmCAT-seq", "EpiAssay:snmCAT-seq"),
        ("sci-MET", "EpiAssay:sci-MET"),
        ("scChIP-seq", "EpiAssay:scChIP-seq"),
        ("scNOMe-seq", "EpiAssay:scNOMe-seq"),
        ("scHi-C", "EFO:0009979"),
    ],
)
def test_enrich_claim_resolves_single_cell_assay_names(raw, expected_curie):
    a = SingleCellEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": raw,
        "qualifier": "marker_peak_of_cell_type",
        "object": "study",
    })
    assert out["ontology_terms"]["subject"] == expected_curie
    assert out["xrefs"]["subject"]["source"] == "assay"
    assert out["xrefs"]["subject"]["assay_class"] == "single_cell"


def test_enrich_claim_resolves_histone_and_target_slots():
    a = SingleCellEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "H3K27ac",
        "qualifier": "marker_peak_of_cell_type",
        "object": "CTCF",
    })
    assert out["ontology_terms"]["subject"] == "EpiHist:H3K27ac"
    assert out["ontology_terms"]["object"] == "UniProt:P49711"
    assert out["xrefs"]["subject"]["source"] == "histone"
    assert out["xrefs"]["object"]["source"] == "target"


# --- cell-type slot — deferred to v0.2.0 -----------------------------------


@pytest.mark.parametrize(
    "cell_type",
    # Labels NOT in the v0.2.0 IHEC core biosample index — these still
    # pass through un-CURIE'd. The index ships a curated seed (~33 of
    # the most-cited reference biosamples); broader cell-type coverage
    # tracks the CL / UBERON cascade planned in a later v0.2.0 §.
    ["neuron", "hepatocyte", "B cell"],
)
def test_cell_type_label_not_in_ihec_index_passes_through_uncuried(cell_type):
    a = SingleCellEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "some_peak",
        "qualifier": "accessible_in_cell_type",
        "object": cell_type,
    })
    assert "object" not in out["ontology_terms"]


@pytest.mark.parametrize(
    "cell_type, expected_curie",
    [
        # v0.2.0 §4 — these now resolve via IHEC_BIOSAMPLE.
        ("CD8+ T cell", "IHECBio:CD8_T_cell"),
        ("CD4+ T cell", "IHECBio:CD4_T_cell"),
        ("K562", "IHECBio:K562"),
        ("GM12878", "IHECBio:GM12878"),
    ],
)
def test_cell_type_label_in_ihec_index_resolves_to_ihecbio(
    cell_type, expected_curie
):
    """v0.2.0 §4 added the IHEC biosample matrix; cell-type labels
    matching the curated seed now resolve to ``IHECBio:`` CURIEs under
    ``accessible_in_cell_type`` (and the rest of the §4 qualifier
    preferences). This pins the new behaviour.
    """
    a = SingleCellEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "some_peak",
        "qualifier": "accessible_in_cell_type",
        "object": cell_type,
    })
    assert out["ontology_terms"]["object"] == expected_curie


# --- xrefs payload shape ----------------------------------------------------


def test_every_xrefs_entry_has_a_source_discriminator_key():
    a = SingleCellEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "scATAC-seq",
        "qualifier": "marker_peak_of_cell_type",
        "object": "CTCF",
    })
    for slot, payload in out["xrefs"].items():
        assert "source" in payload, (
            f"xrefs[{slot!r}] missing 'source' discriminator: {payload!r}"
        )


# --- pure-function invariants ----------------------------------------------


def test_enrich_claim_does_not_mutate_input_dict():
    a = SingleCellEpigenomicsAdapter()
    original = {
        "subject": "H3K27ac",
        "qualifier": "marker_peak_of_cell_type",
        "object": "CTCF",
    }
    snapshot = copy.deepcopy(original)
    a.enrich_claim(original)
    assert original == snapshot


def test_enrich_claim_returns_a_new_dict_object():
    a = SingleCellEpigenomicsAdapter()
    claim = {
        "subject": "scATAC-seq",
        "qualifier": "marker_peak_of_cell_type",
        "object": "study",
    }
    assert a.enrich_claim(claim) is not claim


def test_enrich_claim_is_deterministic_across_repeated_calls():
    a = SingleCellEpigenomicsAdapter()
    claim = {
        "subject": "H3K27ac",
        "qualifier": "marker_peak_of_cell_type",
        "object": "CTCF",
    }
    assert a.enrich_claim(claim) == a.enrich_claim(claim)
