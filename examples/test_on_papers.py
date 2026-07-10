"""Smoke-test the three epigenomics adapters on claims drawn from two real papers.

Papers:
  P1: Roadmap Epigenomics Consortium (2015), "Integrative analysis of 111
      reference human epigenomes", Nature 518, 317-330.
      doi:10.1038/nature14248
  P2: Hall et al. (2014), "Sex differences in the genome-wide DNA methylation
      pattern and impact on gene expression, microRNA levels and insulin
      secretion in human pancreatic islets", Genome Biology 15:522.

Each claim was hand-extracted from the abstract / methods section to mimic
what an LLM extractor in the ASB pipeline would emit. We feed each claim
through the relevant DomainAdapter and print the enriched output so we can
eyeball whether the adapter is resolving the expected CURIEs.
"""

from __future__ import annotations

from pprint import pprint

import rdflib
from pyshacl import validate

from indicium_adapters_epigenomics import (
    BulkEpigenomicsAdapter,
    MethylationArrayEpigenomicsAdapter,
    SingleCellEpigenomicsAdapter,
)

# ---------------------------------------------------------------------------
# Paper 1 — Roadmap Epigenomics: bulk multi-assay (ChIP-seq, DNase, WGBS, …)
# ---------------------------------------------------------------------------
ROADMAP_CLAIMS = [
    # "profiled for histone modification patterns" — H3K4me3 ChIP-seq is a
    # canonical Roadmap mark; here we encode the mark→assay relationship.
    {
        "id": "roadmap-c1",
        "subject": "H3K4me3",
        "qualifier": "assayed_with",
        "object": "ChIP-seq",
        "source_doi": "10.1038/nature14248",
    },
    # H3K27ac marks active enhancers — same assay relation.
    {
        "id": "roadmap-c2",
        "subject": "H3K27ac",
        "qualifier": "assayed_with",
        "object": "ChIP-seq",
        "source_doi": "10.1038/nature14248",
    },
    # "DNA accessibility … DNase I" — DNase-seq protocol.
    {
        "id": "roadmap-c3",
        "subject": "chromatin_accessibility",
        "qualifier": "assayed_with",
        "object": "DNase-seq",
        "source_doi": "10.1038/nature14248",
    },
    # "bisulfite treatment" — WGBS for DNA methylation profiling.
    {
        "id": "roadmap-c4",
        "subject": "DNA_methylation",
        "qualifier": "quantified_methylation_with",
        "object": "WGBS",
        "source_doi": "10.1038/nature14248",
    },
    # "an additional 16 epigenomes previously reported by … (ENCODE)"
    # — the Roadmap analysis consolidates ENCODE epigenome data.
    {
        "id": "roadmap-c5",
        "subject": "reference_epigenomes_127",
        "qualifier": "consolidates_dataset",
        "object": "ENCSR000DZQ",  # representative ENCODE Series accession
        "source_doi": "10.1038/nature14248",
    },
    # The Roadmap consortium itself deposited data to GEO under GSE16256.
    {
        "id": "roadmap-c6",
        "subject": "roadmap_consortium_data",
        "qualifier": "uses_database",
        "object": "GSE16256",
        "source_doi": "10.1038/nature14248",
    },
    # The IHEC reference-epigenome accession family (Roadmap rolls up into IHEC).
    {
        "id": "roadmap-c7",
        "subject": "roadmap_consortium_data",
        "qualifier": "consolidates_dataset",
        "object": "IHECRE00000001.1",
        "source_doi": "10.1038/nature14248",
    },
    # Marker-mark biology: H3K27ac "marked_by" assertion at a peak.
    {
        "id": "roadmap-c8",
        "subject": "active_enhancer_peak_chr1",
        "qualifier": "marked_by",
        "object": "H3K27ac",
        "source_doi": "10.1038/nature14248",
    },
    # CTCF binding (Roadmap profiles transcription factor occupancy via ChIP-seq).
    {
        "id": "roadmap-c9",
        "subject": "insulator_chr1_1000000",
        "qualifier": "bound_by",
        "object": "CTCF",
        "source_doi": "10.1038/nature14248",
    },
]


# ---------------------------------------------------------------------------
# Paper 2 — Hall et al. 2014: Illumina HumanMethylation450 BeadChip (EWAS)
# ---------------------------------------------------------------------------
HALL_CLAIMS = [
    # "analyzed with the Infinium HumanMethylation450 BeadChip array"
    {
        "id": "hall-c1",
        "subject": "pancreatic_islet_methylome_n87",
        "qualifier": "assayed_with_array",
        "object": "HumanMethylation450",
        "source_doi": "10.1186/s13059-014-0522-z",
    },
    # The same in alias form — the adapter should still recognise it.
    {
        "id": "hall-c2",
        "subject": "pancreatic_islet_methylome_n87",
        "qualifier": "assayed_with_array",
        "object": "450K",
        "source_doi": "10.1186/s13059-014-0522-z",
    },
    # "8,140 individual X-chromosome sites … shows sex-specific differences"
    # — DMP at the AR-locus probe (a real probe on the 450K array).
    {
        "id": "hall-c3",
        "subject": "cg10001023",
        "qualifier": "differentially_methylated_at_probe",
        "object": "AR_locus_x_chromosome",
        "source_doi": "10.1186/s13059-014-0522-z",
    },
    # Genome-wide EWAS association with sex covariate.
    {
        "id": "hall-c4",
        "subject": "cg23577020",
        "qualifier": "epigenome_wide_associated_with",
        "object": "sex",
        "source_doi": "10.1186/s13059-014-0522-z",
    },
]


# ---------------------------------------------------------------------------
# Adapter registry — pick the right adapter per claim by its qualifier.
# ---------------------------------------------------------------------------
ADAPTERS = [
    BulkEpigenomicsAdapter(),
    SingleCellEpigenomicsAdapter(),
    MethylationArrayEpigenomicsAdapter(),
]

QUALIFIER_OWNER = {q: a for a in ADAPTERS for q in a.qualifiers}


def route(claim: dict):
    """Pick the adapter whose qualifiers contain this claim's qualifier."""
    return QUALIFIER_OWNER.get(claim["qualifier"])


def run_section(title: str, claims: list[dict]) -> None:
    bar = "=" * 78
    print(f"\n{bar}\n{title}\n{bar}")
    for claim in claims:
        adapter = route(claim)
        if adapter is None:
            print(f"\n[{claim['id']}] no adapter owns qualifier '{claim['qualifier']}'")
            continue
        enriched = adapter.enrich_claim(claim)
        print(f"\n[{claim['id']}] via {adapter.domain_id}")
        print(f"  subject:   {claim['subject']!r}")
        print(f"  qualifier: {claim['qualifier']!r}")
        print(f"  object:    {claim['object']!r}")
        ot = enriched.get("ontology_terms", {})
        xr = enriched.get("xrefs", {})
        if ot:
            print("  ontology_terms:")
            for slot, curie in ot.items():
                src = xr.get(slot, {}).get("source", "?")
                label = xr.get(slot, {}).get("label", "")
                extra = f"  ({label})" if label else ""
                print(f"    {slot:9s} -> {curie:25s}  [via {src}]{extra}")
        else:
            print("  ontology_terms: <none — no slot resolved>")
        if xr and not ot:
            # E.g. probe IDs which only land in xrefs.
            print("  xrefs:")
            pprint(xr, indent=4, width=88, sort_dicts=False)


def coverage_summary(claims: list[dict], label: str) -> None:
    routed = sum(1 for c in claims if route(c) is not None)
    enriched_slots = 0
    total_str_slots = 0
    for c in claims:
        a = route(c)
        if a is None:
            continue
        out = a.enrich_claim(c)
        ot = out.get("ontology_terms", {})
        for slot in ("subject", "object"):
            if isinstance(c.get(slot), str):
                total_str_slots += 1
                if slot in ot:
                    enriched_slots += 1
    pct = (enriched_slots / total_str_slots * 100) if total_str_slots else 0
    print(
        f"  {label}: routed {routed}/{len(claims)} claims, "
        f"resolved {enriched_slots}/{total_str_slots} slots ({pct:.0f}%)"
    )


def shacl_check(claim: dict, adapter) -> None:
    """Materialise a single claim as RDF and run pyshacl against the adapter."""
    asb = rdflib.Namespace("https://w3id.org/asb#")
    epi = rdflib.Namespace("http://holobiomicslab.org/adapters/epigenomics/")
    enriched = adapter.enrich_claim(claim)
    g = rdflib.Graph()
    node = rdflib.URIRef(f"urn:claim:{claim['id']}")
    g.add((node, rdflib.RDF.type, asb.Claim))
    # Predicate name = the qualifier; object = the CURIE if we have one,
    # else the literal value (which the shape may flag as malformed).
    pred = epi[claim["qualifier"]]
    obj_val = enriched.get("ontology_terms", {}).get("object", claim["object"])
    g.add((node, pred, rdflib.Literal(obj_val)))
    shapes = adapter.shacl_shapes()
    conforms, _, report = validate(data_graph=g, shacl_graph=shapes, inference=None)
    status = "PASS" if conforms else "FAIL"
    print(f"  [{claim['id']}] SHACL via {adapter.domain_id}: {status}")
    if not conforms:
        # Just the first reported violation, trimmed.
        if "Constraint Violation" in report:
            head = report.strip().split("Constraint Violation")[1][:240]
        else:
            head = report[:240]
        print(f"    -> {head.strip()}")


def main() -> None:
    run_section("Paper 1 — Roadmap Epigenomics (bulk multi-assay)", ROADMAP_CLAIMS)
    run_section("Paper 2 — Hall et al. 2014 (Illumina 450K EWAS)", HALL_CLAIMS)

    print("\n" + "=" * 78)
    print("Resolution summary")
    print("=" * 78)
    coverage_summary(ROADMAP_CLAIMS, "Roadmap")
    coverage_summary(HALL_CLAIMS, "Hall")

    print("\n" + "=" * 78)
    print("SHACL spot-checks")
    print("=" * 78)
    # Take one well-formed enriched claim per adapter and validate it.
    for claim in [ROADMAP_CLAIMS[5], ROADMAP_CLAIMS[6], HALL_CLAIMS[0]]:
        adapter = route(claim)
        if adapter is not None:
            shacl_check(claim, adapter)


if __name__ == "__main__":
    main()
