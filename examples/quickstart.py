"""indicium-adapters-epigenomics — quickstart example.

Demonstrates the three pieces of public surface that downstream consumers
typically wire up:

1. ``indicium_adapters.discover_adapters()`` — entry-point discovery returns
   every ``DomainAdapter`` registered under ``indicium.adapters`` across all
   installed adapter wheels, including the three from this package.
2. ``adapter.enrich_claim(claim)`` — pure-function enrichment of a draft
   claim dict; populates ``ontology_terms`` (per-slot CURIE) and ``xrefs``
   (per-slot resolver metadata + ``source`` discriminator).
3. ``adapter.shacl_shapes()`` + ``pyshacl.validate(...)`` — SHACL validation
   of the enriched claim against the bundled per-qualifier CURIE-prefix
   constraints.

Run it from the package root:

    python examples/quickstart.py
"""
from __future__ import annotations

from pprint import pprint

import pyshacl
import rdflib
from indicium_adapters import discover_adapters

from indicium_adapters_epigenomics import (
    BulkEpigenomicsAdapter,
    MethylationArrayEpigenomicsAdapter,
    SingleCellEpigenomicsAdapter,
)


def _section(title: str) -> None:
    bar = "=" * 72
    print(f"\n{bar}\n{title}\n{bar}")


def main() -> None:
    # --- 1. Entry-point discovery ------------------------------------------
    _section("1. discover_adapters() — every registered DomainAdapter")
    discovered = discover_adapters()
    for domain_id in sorted(discovered):
        print(f"  - {domain_id}")

    # --- 2. Per-adapter enrichment -----------------------------------------
    bulk = BulkEpigenomicsAdapter()
    sc = SingleCellEpigenomicsAdapter()
    array = MethylationArrayEpigenomicsAdapter()

    bulk_claim = {
        "id": "claim-bulk-1",
        "subject": "H3K27ac",
        "qualifier": "marked_by",
        "object": "CTCF",
    }
    sc_claim = {
        "id": "claim-sc-1",
        "subject": "scATAC-seq",
        "qualifier": "marker_peak_of_cell_type",
        "object": "CTCF",
    }
    array_claim = {
        "id": "claim-array-1",
        "subject": "cg00000029",
        "qualifier": "differentially_methylated_at_probe",
        "object": "EPIC v2",
    }

    enriched_bulk = bulk.enrich_claim(bulk_claim)
    enriched_sc = sc.enrich_claim(sc_claim)
    enriched_array = array.enrich_claim(array_claim)

    _section("2a. BulkEpigenomicsAdapter.enrich_claim")
    pprint(enriched_bulk)
    _section("2b. SingleCellEpigenomicsAdapter.enrich_claim")
    pprint(enriched_sc)
    _section("2c. MethylationArrayEpigenomicsAdapter.enrich_claim")
    pprint(enriched_array)

    # --- 3. SHACL validation against the bundled shapes --------------------
    _section("3. SHACL validation of a well-formed bulk-claim")
    # Build a tiny RDF graph that models a single ``asb:Claim`` whose
    # ``epi:marked_by`` property carries the histone MOD: CURIE the bulk
    # adapter just resolved for H3K27ac.
    marked_by_curie = enriched_bulk["ontology_terms"]["subject"]
    data_ttl = f"""
    @prefix asb: <http://holobiomicslab.org/asb/> .
    @prefix epi: <http://holobiomicslab.org/adapters/epigenomics/> .
    <urn:{enriched_bulk["id"]}> a asb:Claim ; epi:marked_by "{marked_by_curie}" .
    """
    data_graph = rdflib.Graph()
    data_graph.parse(data=data_ttl, format="turtle")

    conforms, _report_graph, report_text = pyshacl.validate(
        data_graph=data_graph,
        shacl_graph=bulk.shacl_shapes(),
    )
    print(f"  marked_by CURIE: {marked_by_curie}")
    print(f"  conforms: {conforms}")
    if not conforms:
        print(report_text)


if __name__ == "__main__":
    main()
