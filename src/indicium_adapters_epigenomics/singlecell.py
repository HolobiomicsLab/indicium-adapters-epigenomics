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

from indicium_adapters_epigenomics._assays import techniques_for_class
from indicium_adapters_epigenomics._namespaces import EPIGENOMICS_PREFIXES
from indicium_adapters_epigenomics._resolver import justification_for, resolve

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

    Every qualifier in the ``qualifiers`` frozenset MUST match
    ``^[a-z][a-z0-9_]*$`` (snake_case). Required by indicium >= 1.7.0:
    Perspicacite routes adapter qualifiers to the open
    ``asb:domainQualifier`` slot, whose SHACL shape enforces exactly this
    pattern. Non-conforming values would be silently rejected by
    ``indicium.validate_graph`` at the RDF boundary. The constraint is
    enforced by ``tests/test_protocol_conformance.py``.
    """

    domain_id = "single-cell-epigenomics"
    api_version = 1

    # ASB's classification layer reads ``primary_domain`` in preference to
    # ``domain_id`` and unions the declared vocabulary into its controlled
    # terms (agentic_science_builder.classification._adapter_vocab). Without
    # it an epigenomics paper classifies as "other" and carries no technique.
    # ``techniques`` is derived from the assay table's own rows for this
    # adapter's ``assay_class``, so the vocabulary moves with the resolver
    # rather than with a second hand-maintained list.
    primary_domain = "epigenomics"
    subdomains: tuple[str, ...] = (domain_id,)
    techniques: tuple[str, ...] = techniques_for_class("single_cell")

    # Every qualifier MUST match ``^[a-z][a-z0-9_]*$`` (snake_case).
    # Required by indicium >= 1.7.0: at serialisation time, Perspicacite
    # routes any qualifier from this set to ``asb:domainQualifier`` (the
    # open slot for domain vocabularies), which enforces the snake_case
    # pattern at the SHACL boundary. Non-conforming values would be
    # rejected by ``indicium.validate_graph``. Enforced for the package
    # by ``tests/test_protocol_conformance.py``.
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
        # v0.2.0 §2 — cell-type-resolved segmentation-state markers.
        # Subject is a cell type / cluster, object is a ChromHMM or
        # Segway state label characteristic of that population.
        "marker_state_of_cell_type",
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
            "co_accessible_with, localized_in_spatial_domain, "
            "marker_state_of_cell_type. Subjects/objects should "
            "typically be CL cell-type labels, UBERON tissue labels, "
            "histone marks, chromatin-target gene symbols, assay/tool "
            "names, cluster identifiers like 'cluster_4', or ChromHMM / "
            "Segway segmentation labels like TssA, EnhA1, ReprPC, Quies "
            "(for marker_state_of_cell_type)."
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
        result.setdefault("ontology_term_justifications", {})
        result.setdefault("xrefs", {})

        qualifier_hint = (
            result["qualifier"] if isinstance(result.get("qualifier"), str) else None
        )

        for slot in _ENRICHABLE_SLOTS:
            raw = result.get(slot)
            if not isinstance(raw, str):
                continue
            value = raw.strip()
            if not value:
                continue

            match, source = resolve(value, qualifier_hint)
            if match is None:
                continue

            result["ontology_terms"][slot] = match["curie"]
            result["ontology_term_justifications"][slot] = justification_for(source)
            result["xrefs"][slot] = {**match, "source": source}

        return result

    def shacl_shapes(self) -> rdflib.Graph:
        """Return ASB-claim SHACL shapes (CURIE format validation only).

        Shares the epigenomics shapes TTL with the bulk and array adapters.
        Returns an empty graph until the TTL is added by the follow-up task.
        """
        return _load_shapes()


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
