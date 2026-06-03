"""MethylationArrayEpigenomicsAdapter — Illumina methylation array adapter.

Covers methylation-array claims for the Illumina platforms:

- HumanMethylation450 (450K)
- MethylationEPIC v1
- MethylationEPIC v2
- MethylationEPIC mouse

and the EWAS workflow vocabulary built around them (probe filtering,
normalisation, cell-composition correction, EWAS testing, differentially
methylated probe / region calling, epigenetic-age clock estimation).

Disjoint from ``BulkEpigenomicsAdapter`` and ``SingleCellEpigenomicsAdapter``
on purpose — array-specific tooling (450K / EPIC array models, probe
filtering, Houseman cell-composition correction, Horvath / Hannum / GrimAge /
PhenoAge clocks) does not generalise to sequencing-based assays and so lives
here.
"""
from __future__ import annotations

import copy
import re
from importlib.resources import files
from typing import Final

import rdflib

from indicium_adapters_epigenomics._namespaces import EPIGENOMICS_PREFIXES
from indicium_adapters_epigenomics._resolver import justification_for, resolve

__all__ = ["MethylationArrayEpigenomicsAdapter"]


_ENRICHABLE_SLOTS: Final[tuple[str, ...]] = ("subject", "object", "qualifier")

# Illumina CpG-probe identifier pattern. The standard "cg" probes are the
# bulk of the array; "ch" (CpH) and "rs" (control SNP) probes are deliberately
# out of scope for this first pass.
_PROBE_ID_RE: Final[re.Pattern[str]] = re.compile(r"^cg\d{8}$")


class MethylationArrayEpigenomicsAdapter:
    """Methylation-array epigenomics DomainAdapter (api_version=1).

    Qualifiers cover the array-specific workflow vocabulary: which array
    model was used, what probe-filtering / normalisation / cell-composition-
    correction / EWAS-test / DMP-calling / clock-estimator method was
    applied. Biology-level qualifiers (binding, accessibility, 3D
    interaction) and cell-resolution qualifiers (cluster enrichment, spatial
    localisation) live in the sibling adapters to keep the three qualifier
    sets pairwise disjoint.

    Enrichment recognises repository accessions, assay names, reference-
    genome assemblies, histone marks, chromatin targets, QC-metric names, and
    Illumina CpG-probe identifiers (``cg\\d{8}``).

    Design decision (v0.1.0): Illumina ``cg`` probes are surfaced via
    ``xrefs[slot]`` with a ``"source": "probe_id"`` discriminator but are
    *not* given an ``ontology_terms`` CURIE — minting a new ``EPICpg`` prefix
    in ``_namespaces.py`` is deferred until a stable identifiers.org or
    upstream registry endpoint is chosen. Consumers that need a CURIE can
    construct one from the bundled ``bare_id`` + ``platform`` fields.

    Every qualifier in the ``qualifiers`` frozenset MUST match
    ``^[a-z][a-z0-9_]*$`` (snake_case). Required by indicium >= 1.7.0:
    Perspicacite routes adapter qualifiers to the open
    ``asb:domainQualifier`` slot, whose SHACL shape enforces exactly this
    pattern. Non-conforming values would be silently rejected by
    ``indicium.validate_graph`` at the RDF boundary. The constraint is
    enforced by ``tests/test_protocol_conformance.py``.
    """

    domain_id = "methylation-array-epigenomics"
    api_version = 1

    # Every qualifier MUST match ``^[a-z][a-z0-9_]*$`` (snake_case).
    # Required by indicium >= 1.7.0: at serialisation time, Perspicacite
    # routes any qualifier from this set to ``asb:domainQualifier`` (the
    # open slot for domain vocabularies), which enforces the snake_case
    # pattern at the SHACL boundary. Non-conforming values would be
    # rejected by ``indicium.validate_graph``. Enforced for the package
    # by ``tests/test_protocol_conformance.py``.
    qualifiers: frozenset[str] = frozenset({
        "assayed_with_array",
        "probe_filtered_with",
        "normalized_with_array_method",
        "corrected_for_cell_composition_with",
        "epigenome_wide_associated_with",
        "differentially_methylated_at_probe",
        "clock_age_estimated_with",
        "tested_with_array_model",
    })

    ontology_prefixes: dict[str, str] = EPIGENOMICS_PREFIXES

    def extraction_context(self) -> str:
        return (
            "Domain: methylation arrays (Illumina 450K, EPIC v1, EPIC v2, "
            "MethylationEPIC mouse; EWAS workflows). Use one of the "
            "following qualifiers: assayed_with_array, probe_filtered_with, "
            "normalized_with_array_method, "
            "corrected_for_cell_composition_with, "
            "epigenome_wide_associated_with, "
            "differentially_methylated_at_probe, clock_age_estimated_with, "
            "tested_with_array_model. Subjects/objects: array model names "
            "(450K, EPIC), normalization methods (minfi/SeSAMe/ChAMP/"
            "RnBeads/ewastools/NOOB/BMIQ/SWAN/funnorm), probe identifiers "
            "like cg12345678, CpG/DMP sets, reference cell-composition "
            "methods (Houseman), epigenetic age estimators (Horvath, "
            "Hannum, GrimAge, PhenoAge). Cell types and tissues should be "
            "CL/UBERON labels."
        )

    def enrich_claim(self, claim: dict) -> dict:
        """Enrich a claim with methylation-array CURIEs and probe xrefs.

        Per-slot resolution order (first non-None match wins):

        1. ``accession_to_curie`` — repository accessions.
        2. ``ASSAY_INDEX`` — array-platform / EWAS-tool names.
        3. ``ASSEMBLY_INDEX`` — reference-genome assembly names.
        4. ``HISTONE_INDEX`` — histone PTMs / variants.
        5. ``TARGET_INDEX`` — TFs and chromatin modifiers.
        6. ``QC_INDEX`` — QC-metric names.
        7. Illumina probe-id regex (``cg\\d{8}``) — xrefs only, no CURIE.

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
            if match is not None:
                result["ontology_terms"][slot] = match["curie"]
                result["ontology_term_justifications"][slot] = justification_for(source)
                result["xrefs"][slot] = {**match, "source": source}
                continue

            # Probe-id detection — xrefs only, no ontology_terms CURIE. See
            # class docstring for the rationale.
            if _PROBE_ID_RE.match(value):
                result["xrefs"][slot] = {
                    "source": "probe_id",
                    "bare_id": value,
                    "platform": "illumina_methylation_array",
                }

        return result

    def shacl_shapes(self) -> rdflib.Graph:
        """Return ASB-claim SHACL shapes (CURIE format validation only).

        Shares the epigenomics shapes TTL with the bulk and single-cell
        adapters. Returns an empty graph until the TTL is added by the
        follow-up task.
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
