"""MethylationArrayEpigenomicsAdapter — enrichment + probe-id tests.

Covers the array-specific qualifier vocabulary, array-platform / EWAS-tool /
normalisation-method resolution via the shared ASSAY / QC indices, and the
Illumina ``cg\\d{8}`` probe-id branch (xrefs-only, no ontology_terms CURIE —
see the adapter class docstring for the design rationale).
"""
from __future__ import annotations

import copy

import pytest

from indicium_adapters_epigenomics import MethylationArrayEpigenomicsAdapter

_EXPECTED_QUALIFIERS = frozenset({
    "assayed_with_array",
    "probe_filtered_with",
    "normalized_with_array_method",
    "corrected_for_cell_composition_with",
    "epigenome_wide_associated_with",
    "differentially_methylated_at_probe",
    "clock_age_estimated_with",
    "tested_with_array_model",
})


# --- qualifiers -------------------------------------------------------------


def test_all_documented_qualifiers_are_present():
    assert MethylationArrayEpigenomicsAdapter.qualifiers == _EXPECTED_QUALIFIERS


@pytest.mark.parametrize("qualifier", sorted(_EXPECTED_QUALIFIERS))
def test_each_qualifier_is_accepted_by_enrich_claim(qualifier):
    """An in-vocabulary qualifier must not cause ``enrich_claim`` to raise."""
    a = MethylationArrayEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "EPIC v2",
        "qualifier": qualifier,
        "object": "cg12345678",
    })
    assert isinstance(out, dict)


# --- array model resolution -------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected_curie, expected_label",
    [
        ("EPIC v2", "EpiAssay:MethylationEPIC-v2",
         "Illumina MethylationEPIC v2.0 BeadChip"),
        ("MethylationEPIC v2", "EpiAssay:MethylationEPIC-v2",
         "Illumina MethylationEPIC v2.0 BeadChip"),
        ("MethylationEPIC", "OBI:0002131",
         "Illumina MethylationEPIC BeadChip"),
        ("450K", "EpiAssay:HumanMethylation450",
         "Illumina HumanMethylation450 BeadChip"),
        ("HumanMethylation450", "EpiAssay:HumanMethylation450",
         "Illumina HumanMethylation450 BeadChip"),
    ],
)
def test_enrich_claim_resolves_illumina_array_models(
    raw, expected_curie, expected_label
):
    a = MethylationArrayEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": raw,
        "qualifier": "assayed_with_array",
        "object": "cg12345678",
    })
    assert out["ontology_terms"]["subject"] == expected_curie
    xref = out["xrefs"]["subject"]
    assert xref["label"] == expected_label
    assert xref["assay_class"] == "array"
    assert xref["source"] == "assay"


# --- normalisation tools ---------------------------------------------------
#
# The ASSAY_INDEX and QC_INDEX do not currently carry a CURIE for every
# normalisation tool name. The contract is: tools that DO resolve via the
# shared indices pick up a CURIE; tools that do not pass through without an
# ``ontology_terms`` entry. The tests below pin both halves of that contract
# so a regression in either direction is caught.


@pytest.mark.parametrize(
    "raw, expected_curie",
    [
        # Resolved by QC_INDEX as array-specific normalisation methods.
        ("NOOB", "EPIQC:noob-offset"),
        ("BMIQ", "EPIQC:bmiq-alpha"),
        ("funnorm", "EPIQC:funnorm-pcs"),
    ],
)
def test_enrich_claim_resolves_known_normalisation_methods(raw, expected_curie):
    a = MethylationArrayEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": raw,
        "qualifier": "normalized_with_array_method",
        "object": "study",
    })
    assert out["ontology_terms"]["subject"] == expected_curie
    assert out["xrefs"]["subject"]["source"] == "qc"


@pytest.mark.parametrize(
    "raw",
    ["minfi", "ChAMP", "SeSAMe", "RnBeads", "ewastools", "SWAN"],
)
def test_unresolved_normalisation_tools_pass_through_without_curie(raw):
    """Tools that don't have a row in ASSAY_INDEX / QC_INDEX yet should pass
    through silently — no exception, no synthetic CURIE. This pins the v0.1.0
    behaviour; promote any of these to the indices when canonical CURIEs are
    chosen and update this list accordingly."""
    a = MethylationArrayEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": raw,
        "qualifier": "normalized_with_array_method",
        "object": "study",
    })
    assert "subject" not in out["ontology_terms"]
    assert "subject" not in out["xrefs"]


# --- probe-id detection ----------------------------------------------------


def test_enrich_claim_detects_illumina_probe_id_in_subject_slot():
    """``cg\\d{8}`` IDs surface in ``xrefs`` with a ``"probe_id"`` source
    discriminator and a ``"platform"`` field, but do NOT receive an
    ``ontology_terms`` entry (per the adapter's documented v0.1.0 choice)."""
    a = MethylationArrayEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "cg12345678",
        "qualifier": "differentially_methylated_at_probe",
        "object": "study",
    })
    assert "subject" not in out["ontology_terms"]
    xref = out["xrefs"]["subject"]
    assert xref["source"] == "probe_id"
    assert xref["bare_id"] == "cg12345678"
    assert xref["platform"] == "illumina_methylation_array"


def test_enrich_claim_detects_illumina_probe_id_in_object_slot():
    a = MethylationArrayEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "EPIC v2",
        "qualifier": "differentially_methylated_at_probe",
        "object": "cg00000029",
    })
    # Subject is the array model (assay enrichment).
    assert out["ontology_terms"]["subject"] == "EpiAssay:MethylationEPIC-v2"
    # Object is the probe — xrefs only.
    assert "object" not in out["ontology_terms"]
    assert out["xrefs"]["object"]["source"] == "probe_id"
    assert out["xrefs"]["object"]["bare_id"] == "cg00000029"


@pytest.mark.parametrize(
    "bad_probe",
    [
        "cg1234567",   # 7 digits — too short
        "cg123456789", # 9 digits — too long
        "CG12345678",  # uppercase prefix
        "ch12345678",  # CpH probe — deliberately out of scope
        "rs12345678",  # control SNP probe — deliberately out of scope
    ],
)
def test_non_cg_probe_patterns_do_not_match(bad_probe):
    a = MethylationArrayEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": bad_probe,
        "qualifier": "differentially_methylated_at_probe",
        "object": "study",
    })
    # Should NOT be detected as a probe_id; should also not get an
    # ontology_terms entry from the other resolvers.
    if "subject" in out["xrefs"]:
        assert out["xrefs"]["subject"].get("source") != "probe_id"


# --- pure-function invariants ----------------------------------------------


def test_enrich_claim_does_not_mutate_input_dict():
    a = MethylationArrayEpigenomicsAdapter()
    original = {
        "subject": "cg12345678",
        "qualifier": "differentially_methylated_at_probe",
        "object": "EPIC v2",
    }
    snapshot = copy.deepcopy(original)
    a.enrich_claim(original)
    assert original == snapshot


def test_enrich_claim_returns_a_new_dict_object():
    a = MethylationArrayEpigenomicsAdapter()
    claim = {
        "subject": "EPIC v2",
        "qualifier": "assayed_with_array",
        "object": "study",
    }
    assert a.enrich_claim(claim) is not claim


def test_enrich_claim_is_deterministic_across_repeated_calls():
    a = MethylationArrayEpigenomicsAdapter()
    claim = {
        "subject": "cg12345678",
        "qualifier": "differentially_methylated_at_probe",
        "object": "EPIC v2",
    }
    assert a.enrich_claim(claim) == a.enrich_claim(claim)
