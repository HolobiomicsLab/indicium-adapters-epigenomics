"""SingleCellEpigenomicsAdapter — single-cell and single-nucleus epigenomics adapter.

Covers cell-resolution epigenomic claims for chromatin accessibility, DNA
methylation, histone modification, chromatin binding, and 3D conformation at
single-cell or single-nucleus level:

- scATAC-seq, snATAC-seq
- snmC-seq, sci-MET, snmCAT-seq
- scChIP-seq, single-cell CUT&Tag
- scNOMe-seq
- scHi-C

Disjoint from ``BulkEpigenomicsAdapter`` and
``MethylationArrayEpigenomicsAdapter``: bulk biology/assay qualifiers and
methylation-array-specific qualifiers live in their respective siblings.
This adapter holds only the qualifiers that require per-cell, per-nucleus,
per-cluster or per-spatial-domain resolution.
"""
from __future__ import annotations

import copy
from importlib.resources import files
from typing import Final

import rdflib

from indicium_adapters_epigenomics._enrichment import (
    ASSAY_INDEX,
    ASSEMBLY_INDEX,
    HISTONE_INDEX,
    QC_INDEX,
    TARGET_INDEX,
)
from indicium_adapters_epigenomics._namespaces import EPIGENOMICS_PREFIXES
from indicium_adapters_epigenomics._registries import accession_to_curie

__all__ = ["SingleCellEpigenomicsAdapter"]


_ENRICHABLE_SLOTS: Final[tuple[str, ...]] = ("subject", "object", "qualifier")


class SingleCellEpigenomicsAdapter:
    """Single-cell / single-nucleus epigenomics DomainAdapter (api_version=1).

    Qualifiers are deliberately scoped to assertions that require cell-,
    nucleus-, cluster- or spatial-domain resolution. Bulk biology qualifiers
    (``bound_by``, ``marked_by``, ``accessible_in``, ``methylated_in``, ...)
    stay in ``BulkEpigenomicsAdapter`` to preserve disjointness across
    siblings.

    Enrichment in v0.1.0 reuses the shared assay / assembly / histone-mark /
    chromatin-target / QC / repository-accession indices. Cell-type (CL) and
    anatomy (UBERON) enrichment for the ``*_in_cell_type`` qualifiers is
    deferred to v0.2.0 — see ``ROADMAP.md`` — once a reference paper anchors
    the bundled term set.
    """

    domain_id = "single-cell-epigenomics"
    api_version = 1

    qualifiers: frozenset[str] = frozenset({
        # Cell-resolution qualifiers — disjoint from the bulk and array sets.
        "accessible_in_cell_type",
        "methylated_in_cell_type",
        "bound_in_cell_type",
        "marker_peak_of_cell_type",
        "marker_dmr_of_cell_type",
        "enriched_in_cluster",
        "co_accessible_with",
        "localized_in_spatial_domain",
    })

    ontology_prefixes: dict[str, str] = EPIGENOMICS_PREFIXES

    def extraction_context(self) -> str:
        return (
            "Domain: single-cell epigenomics (scATAC-seq, snATAC-seq, "
            "snmC-seq, sci-MET, snmCAT-seq, scChIP-seq, scNOMe-seq, scHi-C, "
            "single-cell CUT&Tag). Use one of the following qualifiers: "
            "accessible_in_cell_type, methylated_in_cell_type, "
            "bound_in_cell_type, marker_peak_of_cell_type, "
            "marker_dmr_of_cell_type, enriched_in_cluster, "
            "co_accessible_with, localized_in_spatial_domain. "
            "Subjects/objects should typically be CL cell-type labels, "
            "UBERON tissue labels, histone marks, chromatin-target gene "
            "symbols, assay/tool names, or cluster identifiers like "
            "'cluster_4'."
        )

    def enrich_claim(self, claim: dict) -> dict:
        """Enrich a claim with epigenomics CURIEs.

        Per-slot resolution order (first non-None match wins):

        1. ``accession_to_curie`` — repository accessions.
        2. ``ASSAY_INDEX`` — single-cell epigenomics assay names.
        3. ``ASSEMBLY_INDEX`` — reference-genome assembly names.
        4. ``HISTONE_INDEX`` — histone PTMs / variants.
        5. ``TARGET_INDEX`` — TFs and chromatin modifiers.
        6. ``QC_INDEX`` — QC-metric names.

        TODO: cell-type (CL) and anatomy (UBERON) lookup land in v0.2.0; no
        ``CL_INDEX`` / ``UBERON_INDEX`` exists in ``_enrichment.py`` yet, so
        the cell-type-biased branches are intentionally absent. Subjects /
        objects carrying a free-text cell-type label will pass through
        un-CURIE'd until then.

        Pure function — never mutates the input dict.
        """
        result = copy.deepcopy(claim)
        result.setdefault("ontology_terms", {})
        result.setdefault("xrefs", {})

        for slot in _ENRICHABLE_SLOTS:
            raw = result.get(slot)
            if not isinstance(raw, str):
                continue
            value = raw.strip()
            if not value:
                continue

            match, source = _resolve(value)
            if match is None:
                continue

            result["ontology_terms"][slot] = match["curie"]
            result["xrefs"][slot] = {**match, "source": source}

        return result

    def shacl_shapes(self) -> rdflib.Graph:
        """Return ASB-claim SHACL shapes (CURIE format validation only).

        Shares the epigenomics shapes TTL with the bulk and array adapters.
        Returns an empty graph until the TTL is added by the follow-up task.
        """
        return _load_shapes()


def _resolve(value: str) -> tuple[dict[str, str] | None, str]:
    """Same precedence as bulk; cell-type / anatomy branches are deferred."""
    if match := accession_to_curie(value):
        return match, "accession"
    if match := ASSAY_INDEX.lookup(value):
        return match, "assay"
    if match := ASSEMBLY_INDEX.lookup(value):
        return match, "assembly"
    if match := HISTONE_INDEX.lookup(value):
        return match, "histone"
    if match := TARGET_INDEX.lookup(value):
        return match, "target"
    if match := QC_INDEX.lookup(value):
        return match, "qc"
    return None, ""


def _load_shapes() -> rdflib.Graph:
    """Shared shapes loader. Returns an empty graph if the TTL is absent."""
    g = rdflib.Graph()
    ttl = files("indicium_adapters_epigenomics.shapes").joinpath(
        "epigenomics.shapes.ttl"
    )
    # TODO: load shapes — TTL is produced by the follow-up shapes task.
    if ttl.is_file():
        g.parse(ttl.open("rb"), format="turtle")
    return g
