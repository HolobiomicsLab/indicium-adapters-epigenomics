"""BulkEpigenomicsAdapter — enrichment pipeline tests.

Exercises every enrichment branch in ``BulkEpigenomicsAdapter.enrich_claim``:
repository-accession resolution, assay / assembly / histone-mark / chromatin-
target / QC-metric lookups, multi-slot population, xrefs payload shape, and
pure-function invariants (no input mutation, deterministic, no I/O).
"""
from __future__ import annotations

import copy

import pytest

from indicium_adapters_epigenomics import BulkEpigenomicsAdapter
from indicium_adapters_epigenomics._registries import accession_to_curie

# --- constants --------------------------------------------------------------


_QUALIFIER_OK = "marked_by"  # arbitrary valid qualifier for the test claims


# --- helpers ----------------------------------------------------------------


def _make_claim(**slots):
    claim = {"id": "c1", "qualifier": _QUALIFIER_OK}
    claim.update(slots)
    return claim


# --- accession enrichment (per-repo parametrised) ---------------------------


@pytest.mark.parametrize(
    "raw, expected_curie, expected_repository",
    [
        ("GSE123456", "GEO:GSE123456", "GEO"),
        ("GSM4567890", "GEO:GSM4567890", "GEO"),
        ("SRR1234567", "SRA:SRR1234567", "SRA"),
        ("SRX1234567", "SRA:SRX1234567", "SRA"),
        ("ERR1234567", "ENA:ERR1234567", "ENA"),
        ("ENCSR000ABC", "ENCODE:ENCSR000ABC", "ENCODE"),
        ("ENCFF123XYZ", "ENCODE:ENCFF123XYZ", "ENCODE"),
        ("IHECRE00000001.1", "EpiRR:IHECRE00000001.1", "IHEC"),
        ("EGAS00000000001", "EGA:EGAS00000000001", "EGA"),
        ("EGAD00000000001", "EGA:EGAD00000000001", "EGA"),
        ("PRJNA12345", "BioProject:PRJNA12345", "BioProject"),
        ("PRJEB104129", "BioProject:PRJEB104129", "BioProject"),
        ("E-MTAB-1234", "ArrayExpress:E-MTAB-1234", "ArrayExpress"),
        ("ChIP-Atlas:12345", "ChIPAtlas:12345", "ChIP-Atlas"),
        ("CA98765", "ChIPAtlas:98765", "ChIP-Atlas"),
        ("Cistrome:45123", "Cistrome:45123", "Cistrome"),
    ],
)
def test_enrich_claim_resolves_repository_accession(
    raw, expected_curie, expected_repository
):
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim(_make_claim(subject=raw))
    assert out["ontology_terms"]["subject"] == expected_curie
    xref = out["xrefs"]["subject"]
    assert xref["curie"] == expected_curie
    assert xref["repository"] == expected_repository
    assert xref["source"] == "accession"


def test_accession_to_curie_returns_none_for_unrecognised_token():
    assert accession_to_curie("not-an-accession") is None
    assert accession_to_curie("") is None


# --- SEMAPV mapping-justification (indicium-adapters A-convention) -----------


def test_enrich_claim_sets_composite_justification_for_accession():
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim(_make_claim(subject="ENCSR000ABC"))
    assert out["ontology_terms"]["subject"] == "ENCODE:ENCSR000ABC"
    assert out["ontology_term_justifications"]["subject"] == "semapv:CompositeMatching"


def test_enrich_claim_sets_lexical_justification_for_curated_lookup():
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim(_make_claim(subject="ChIP-seq"))
    assert out["xrefs"]["subject"]["source"] == "assay"
    assert out["ontology_term_justifications"]["subject"] == "semapv:LexicalMatching"


def test_justification_keys_mirror_ontology_terms_keys():
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim(_make_claim(subject="ENCSR000ABC", object="ChIP-seq"))
    assert set(out["ontology_term_justifications"]) == set(out["ontology_terms"])


def test_no_justification_when_nothing_resolves():
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim(_make_claim(subject="an unresolvable free-text value"))
    assert out["ontology_terms"] == {}
    assert out["ontology_term_justifications"] == {}


# --- assay enrichment -------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected_curie",
    [
        ("ChIP-seq", "OBI:0000716"),
        ("ChIPseq", "OBI:0000716"),
        ("ATAC-seq", "EFO:0007045"),
        ("WGBS", "OBI:0001863"),
        ("Hi-C", "EFO:0007693"),
        ("Micro-C", "EpiAssay:Micro-C"),
        ("CUT&RUN", "EFO:0009973"),
        ("CUT&Tag", "EFO:0920004"),
        ("MNase-seq", "OBI:0001924"),
        ("DNase-seq", "OBI:0001853"),
        ("RRBS", "EFO:0008900"),
        ("EM-seq", "EpiAssay:EM-seq"),
    ],
)
def test_enrich_claim_resolves_assay_name(raw, expected_curie):
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim(_make_claim(subject=raw))
    assert out["ontology_terms"]["subject"] == expected_curie
    assert out["xrefs"]["subject"]["source"] == "assay"


# --- assembly enrichment ----------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected_curie, expected_organism",
    [
        ("GRCh38", "RefSeq:GCF_000001405.40", "NCBITaxon:9606"),
        ("hg38", "RefSeq:GCF_000001405.40", "NCBITaxon:9606"),
        ("GRCh37", "RefSeq:GCF_000001405.25", "NCBITaxon:9606"),
        ("hg19", "RefSeq:GCF_000001405.25", "NCBITaxon:9606"),
        ("mm10", "RefSeq:GCF_000001635.26", "NCBITaxon:10090"),
        ("mm39", "RefSeq:GCF_000001635.27", "NCBITaxon:10090"),
        ("T2T-CHM13v2", "RefSeq:GCF_009914755.1", "NCBITaxon:9606"),
        ("hs1", "RefSeq:GCF_009914755.1", "NCBITaxon:9606"),
    ],
)
def test_enrich_claim_resolves_assembly_name(raw, expected_curie, expected_organism):
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim(_make_claim(subject=raw))
    assert out["ontology_terms"]["subject"] == expected_curie
    xref = out["xrefs"]["subject"]
    assert xref["organism_curie"] == expected_organism
    assert xref["source"] == "assembly"


# --- histone enrichment (notation variants) ---------------------------------


@pytest.mark.parametrize(
    "raw, expected_curie",
    [
        ("H3K27ac", "EpiHist:H3K27ac"),
        ("H3K27Ac", "EpiHist:H3K27ac"),
        ("H3 K27 ac", "EpiHist:H3K27ac"),
        ("H3K27-ac", "EpiHist:H3K27ac"),
        ("H3K4me3", "EpiHist:H3K4me3"),
        ("H3K4Me3", "EpiHist:H3K4me3"),
        ("H3K9me3", "EpiHist:H3K9me3"),
        ("H3K36me3", "EpiHist:H3K36me3"),
    ],
)
def test_enrich_claim_resolves_histone_modification(raw, expected_curie):
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim(_make_claim(subject=raw))
    assert out["ontology_terms"]["subject"] == expected_curie
    assert out["xrefs"]["subject"]["source"] == "histone"


# --- target enrichment ------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected_curie, expected_organism",
    [
        ("CTCF", "UniProt:P49711", "NCBITaxon:9606"),
        ("EZH2", "UniProt:Q15910", "NCBITaxon:9606"),
        ("EP300", "UniProt:Q09472", "NCBITaxon:9606"),
        # ``P300`` is a canonical alias of EP300; the lookup is alias-aware.
        ("P300", "UniProt:Q09472", "NCBITaxon:9606"),
        ("SPI1", "UniProt:P17947", "NCBITaxon:9606"),
    ],
)
def test_enrich_claim_resolves_chromatin_target_to_canonical_human_ortholog(
    raw, expected_curie, expected_organism
):
    """Without an explicit organism hint, the lookup must pick the canonical
    (human) ortholog per TARGET_INDEX's documented row-order policy."""
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim(_make_claim(object=raw))
    assert out["ontology_terms"]["object"] == expected_curie
    xref = out["xrefs"]["object"]
    assert xref["organism_curie"] == expected_organism
    assert xref["source"] == "target"


# --- QC enrichment ----------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected_curie",
    [
        ("FRiP", "EPIQC:frip"),
        ("TSS-enrichment", "EPIQC:tss-enrichment"),
        ("TSS enrichment", "EPIQC:tss-enrichment"),
        ("bisulfite-conversion-rate", "EPIQC:bisulfite-conversion-rate"),
        ("NSC", "EPIQC:nsc"),
        ("RSC", "EPIQC:rsc"),
        ("PBC1", "EPIQC:pbc1"),
    ],
)
def test_enrich_claim_resolves_qc_metric(raw, expected_curie):
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim(_make_claim(subject=raw))
    assert out["ontology_terms"]["subject"] == expected_curie
    assert out["xrefs"]["subject"]["source"] == "qc"


# --- multi-slot population --------------------------------------------------


def test_enrich_claim_populates_both_subject_and_object_slots():
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "H3K27ac",
        "qualifier": "marked_by",
        "object": "CTCF",
    })
    assert out["ontology_terms"]["subject"] == "EpiHist:H3K27ac"
    assert out["ontology_terms"]["object"] == "UniProt:P49711"
    # Both xrefs carry a ``source`` discriminator that names the resolver that
    # produced them.
    assert out["xrefs"]["subject"]["source"] == "histone"
    assert out["xrefs"]["object"]["source"] == "target"


# --- xrefs payload shape ----------------------------------------------------


def test_every_xrefs_entry_has_a_source_discriminator_key():
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "GSE123456",
        "qualifier": "consolidates_dataset",
        "object": "H3K27ac",
    })
    for slot, payload in out["xrefs"].items():
        assert "source" in payload, (
            f"xrefs[{slot!r}] missing 'source' discriminator: {payload!r}"
        )
        assert isinstance(payload["source"], str) and payload["source"]


# --- pure-function invariants -----------------------------------------------


def test_enrich_claim_does_not_mutate_input_dict():
    a = BulkEpigenomicsAdapter()
    original = {
        "id": "c1",
        "subject": "H3K27ac",
        "qualifier": "marked_by",
        "object": "CTCF",
    }
    snapshot = copy.deepcopy(original)
    a.enrich_claim(original)
    assert original == snapshot


def test_enrich_claim_returns_a_new_dict_object():
    a = BulkEpigenomicsAdapter()
    claim = {"subject": "CTCF", "qualifier": "bound_by", "object": "study"}
    assert a.enrich_claim(claim) is not claim


def test_enrich_claim_is_deterministic_across_repeated_calls():
    a = BulkEpigenomicsAdapter()
    claim = {
        "subject": "H3K27ac",
        "qualifier": "marked_by",
        "object": "CTCF",
    }
    assert a.enrich_claim(claim) == a.enrich_claim(claim)


# --- unrecognised input handling -------------------------------------------


def test_enrich_claim_leaves_unrecognised_slot_value_un_curie_d():
    """An unrecognised free-text value must not produce an ``ontology_terms``
    entry for that slot (and must not raise)."""
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "random_nonsense_xyz",
        "qualifier": "marked_by",
        "object": "another_random_token",
    })
    # The adapter pre-creates the ``ontology_terms`` dict — assert it's
    # empty (no per-slot entries) rather than asserting the dict is absent.
    assert out["ontology_terms"] == {}
    assert out["xrefs"] == {}


@pytest.mark.parametrize("bad_value", [None, 42, 3.14, ["nested"], {"k": "v"}])
def test_enrich_claim_skips_non_string_slot_values_silently(bad_value):
    """Non-string slot values must be skipped, not raise."""
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": bad_value,
        "qualifier": "marked_by",
        "object": "CTCF",
    })
    # Only the well-typed object slot is enriched.
    assert "subject" not in out["ontology_terms"]
    assert out["ontology_terms"]["object"] == "UniProt:P49711"


def test_enrich_claim_handles_empty_dict_without_raising():
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim({})
    assert out["ontology_terms"] == {}
    assert out["xrefs"] == {}


def test_enrich_claim_handles_whitespace_only_slot_values():
    a = BulkEpigenomicsAdapter()
    out = a.enrich_claim({
        "subject": "   ",
        "qualifier": "marked_by",
        "object": "",
    })
    assert out["ontology_terms"] == {}
