"""Shared qualifier-aware resolver for the three epigenomics adapters.

Before v0.2.0 each of ``bulk.py``, ``singlecell.py`` and ``array.py`` carried
its own byte-identical ``_resolve(value)`` helper that walked a *fixed*
six-step priority (accession → assay → assembly → histone → target → qc)
regardless of what the claim's qualifier was. That worked for v0.1.0 because
no token in the bundled indices cross-hits two domains, but several v0.2.0
additions (ChromHMM states, IHEC biosamples, ENCODE cCRE) introduce
overlap with the existing histone / target / assay indices, so resolution
order has to be qualifier-aware.

The pattern (qualifier → preferred-index-tuple map with the fixed list as
fallback) is the reusable structural convention for any future
``DomainAdapter`` whose index family grows past three or four members —
which is why it lives in its own module rather than inline in each
adapter.

Design properties:

- **Accession always wins.** Repository accessions are regex-deterministic
  and outrank every free-text resolution, regardless of qualifier — this
  matches the explicit v0.1.0 contract documented in
  ``BulkEpigenomicsAdapter.enrich_claim()``.
- **Pure function, no I/O.** Mirrors the existing ``_resolve`` discipline.
- **Backwards-compatible default.** When ``qualifier`` is ``None`` or no
  preference is registered for it, the fixed six-step order is preserved,
  so every v0.1.0 behavioural test continues to pass unchanged.
- **Preference shadows the default; it does not replace it.** Sources
  listed in a qualifier preference are tried first in their declared
  order; remaining sources are then tried in the default order. A
  preferred source that doesn't match still falls through to the rest of
  the cascade — preference is a hint, not a hard route.
"""
from __future__ import annotations

from typing import Final

from indicium_adapters_epigenomics._enrichment import (
    ASSAY_INDEX,
    ASSEMBLY_INDEX,
    CCRE_CLASS_INDEX,
    CHROM_STATE_INDEX,
    CONFORMATION_INDEX,
    DNA_MOD_INDEX,
    HISTONE_INDEX,
    IHEC_BIOSAMPLE_INDEX,
    QC_INDEX,
    TARGET_INDEX,
)
from indicium_adapters_epigenomics._registries import accession_to_curie

__all__ = ["resolve", "Source", "SOURCE_TO_SEMAPV", "justification_for"]


# Source labels — also used as keys in `_QUALIFIER_PREFERENCES`. Kept as
# plain strings (not an Enum) to match the existing ``_resolve`` return
# contract used by callers and exercised in test_bulk / test_singlecell /
# test_array.
ACCESSION: Final = "accession"
ASSAY: Final = "assay"
ASSEMBLY: Final = "assembly"
HISTONE: Final = "histone"
TARGET: Final = "target"
QC: Final = "qc"
CHROM_STATE: Final = "chrom_state"
DNA_MOD: Final = "dna_mod"
IHEC_BIOSAMPLE: Final = "ihec_biosample"
CCRE: Final = "ccre"
CONFORMATION: Final = "conformation"

Source = str  # alias for the literal source label


# SEMAPV (https://w3id.org/semapv/vocab/) mapping-justification per source —
# honest provenance for *how* each ontology_terms CURIE was obtained. Accession
# → CURIE registry resolution is identifier (composite) matching; the curated
# controlled-vocabulary lookups match the input name/label string against a
# bundled vocabulary, i.e. lexical matching. Adapters expose this as the
# optional ``ontology_term_justifications`` companion-map documented by the
# indicium-adapters core DomainAdapter Protocol; downstream consumers (e.g.
# Perspicacite) promote ``ontology_terms`` + this map to first-class
# ``sssom:Mapping`` nodes carrying the SEMAPV justification.
SOURCE_TO_SEMAPV: Final[dict[str, str]] = {
    ACCESSION: "semapv:CompositeMatching",
    ASSAY: "semapv:LexicalMatching",
    ASSEMBLY: "semapv:LexicalMatching",
    HISTONE: "semapv:LexicalMatching",
    TARGET: "semapv:LexicalMatching",
    QC: "semapv:LexicalMatching",
    CHROM_STATE: "semapv:LexicalMatching",
    DNA_MOD: "semapv:LexicalMatching",
    IHEC_BIOSAMPLE: "semapv:LexicalMatching",
    CCRE: "semapv:LexicalMatching",
    CONFORMATION: "semapv:LexicalMatching",
}


def justification_for(source: str) -> str:
    """Return the SEMAPV mapping-justification CURIE for a resolver *source*
    label. Falls back to lexical matching for any unmapped source."""
    return SOURCE_TO_SEMAPV.get(source, "semapv:LexicalMatching")


# Fixed default order — extends the v0.1.0 cascade. New indices are
# appended to the END so any value that resolved in v0.1.0 still resolves
# to the same source (regression-protected by test_resolver). The new
# CHROM_STATE / DNA_MOD / IHEC_BIOSAMPLE / CCRE / CONFORMATION sources
# only ever win through a qualifier preference; they never overtake the
# v0.1.0 indices in the default cascade.
_DEFAULT_ORDER: Final[tuple[Source, ...]] = (
    ACCESSION,
    ASSAY,
    ASSEMBLY,
    HISTONE,
    TARGET,
    QC,
    CHROM_STATE,
    DNA_MOD,
    IHEC_BIOSAMPLE,
    CCRE,
    CONFORMATION,
)


# Qualifier → preferred lookup-source tuple. Accession is always tried
# first (regex-deterministic) before the qualifier-aware order kicks in; a
# source not listed here for a given qualifier falls back to the default
# cascade, so the preference is purely a hint about which index to consult
# first when more than one could plausibly match a free-text value.
#
# Mapping rationale (epigenomics-specific):
#
# - Chromatin-binding qualifiers prefer TARGET over HISTONE: an object like
#   ``"EZH2"`` is a Polycomb repressor (target), not a histone PTM, and
#   should never be misresolved to a histone CURIE under ``bound_by``.
# - Histone-marking qualifiers do the inverse: ``"H3K27ac"`` under
#   ``marked_by`` must resolve to the histone PTM, not be mistaken for a
#   protein symbol.
# - Assay-tool qualifiers prefer ASSAY so that ambiguous acronyms like
#   ``"MNase"`` (assay + an enzyme that could appear in a target index)
#   resolve to the assay first.
# - Assembly qualifier prefers ASSEMBLY so version strings like
#   ``"GRCh38"`` go straight to the assembly index.
#
# The map grows in v0.2.0 §2-§6 as new indices ship; the structural
# contract (preference is a partial order, not a hard route) makes those
# additions a straight superset.
_QUALIFIER_PREFERENCES: Final[dict[str, tuple[Source, ...]]] = {
    # Chromatin-binding family. The single-cell ``bound_in_cell_type``
    # qualifier additionally tries IHEC_BIOSAMPLE on the fall-through so
    # a cell-type object (``K562``, ``CD4_T_cell``) still resolves.
    "bound_by": (TARGET, HISTONE),
    "co_occupies_with": (TARGET, HISTONE),
    "differentially_bound_in": (TARGET, HISTONE, IHEC_BIOSAMPLE),
    "bound_in_cell_type": (TARGET, HISTONE, IHEC_BIOSAMPLE),
    # Histone-marking family. ``marker_peak_of_cell_type`` was here in
    # §1; v0.2.0 §4 moves it down to the accessibility / cell-type
    # section because IHEC_BIOSAMPLE is the most likely target for the
    # object slot once the biosample index lands.
    "marked_by": (HISTONE, TARGET),
    "enriched_for": (HISTONE, TARGET),
    # Assay-tool family.
    "assayed_with": (ASSAY,),
    "sequenced_on": (ASSAY,),
    "called_peaks_with": (ASSAY,),
    "quantified_methylation_with": (ASSAY,),
    "normalized_with": (ASSAY,),
    "tested_with": (ASSAY,),
    # Methylation-array workflow qualifiers — same assay-first preference.
    # ``differentially_methylated_at_probe`` is intentionally absent here:
    # it's a methylation qualifier and its preference is set below
    # (DNA_MOD, TARGET) to match the other methylated_* qualifiers,
    # rather than ASSAY-first like the other array-workflow qualifiers.
    "assayed_with_array": (ASSAY,),
    "probe_filtered_with": (ASSAY,),
    "normalized_with_array_method": (ASSAY,),
    "corrected_for_cell_composition_with": (ASSAY,),
    "epigenome_wide_associated_with": (ASSAY,),
    "clock_age_estimated_with": (ASSAY,),
    "tested_with_array_model": (ASSAY,),
    # Assembly qualifier.
    "aligned_with": (ASSEMBLY,),
    # Chromatin-state family (v0.2.0 §2). Both qualifiers route through
    # CHROM_STATE first; any value that doesn't match a known state
    # (e.g. an LLM-extracted descriptor like "open chromatin region")
    # falls through to the rest of the cascade.
    "chromatin_state_of": (CHROM_STATE,),
    "marker_state_of_cell_type": (CHROM_STATE,),
    # Methylation family (v0.2.0 §3). DNA_MOD goes first so that an
    # object like "5mC" or "5hmC" resolves to the CHEBI / EpiMod
    # chemistry CURIE before the histone index gets a look. TARGET is
    # the next hop so that methylation-machinery enzymes (DNMTs, TETs,
    # MBDs) — which live in TARGET — still resolve when the object is
    # an enzyme symbol rather than a chemistry name. IHEC_BIOSAMPLE
    # follows so a biosample object (``methylated_in K562``,
    # ``methylated_in_cell_type CD4+ T cell``) resolves to the IHEC
    # term too. v0.2.0 §4 extension.
    "methylated_in": (DNA_MOD, TARGET, IHEC_BIOSAMPLE),
    "methylated_in_cell_type": (DNA_MOD, TARGET, IHEC_BIOSAMPLE),
    "differentially_methylated_in": (DNA_MOD, TARGET, IHEC_BIOSAMPLE),
    "methylated_at_locus": (DNA_MOD, TARGET, IHEC_BIOSAMPLE),
    "differentially_methylated_at_probe": (DNA_MOD, TARGET, IHEC_BIOSAMPLE),
    # Accessibility / cell-type family (v0.2.0 §4). The object is most
    # commonly an IHEC biosample (cell type or tissue); fall through to
    # the rest of the cascade for the rarer case where the object is an
    # assay name or histone mark.
    "accessible_in": (IHEC_BIOSAMPLE,),
    "accessible_in_cell_type": (IHEC_BIOSAMPLE,),
    "differentially_accessible_in": (IHEC_BIOSAMPLE,),
    "marker_peak_of_cell_type": (IHEC_BIOSAMPLE, HISTONE, TARGET),
    "marker_dmr_of_cell_type": (IHEC_BIOSAMPLE, DNA_MOD),
    "enriched_in_cluster": (IHEC_BIOSAMPLE,),
    "localized_in_spatial_domain": (IHEC_BIOSAMPLE,),
    "co_accessible_with": (IHEC_BIOSAMPLE,),
    # ENCODE SCREEN cCRE classification (v0.2.0 §5). Object is one of
    # PLS / pELS / dELS / CA-CTCF / CA-H3K4me3 / CA-TF / CA-only.
    # ACCESSION still always wins first (an ``EH38E0000001`` accession
    # outranks anything), and the CCRE source acts only as the
    # qualifier-preferred fall-through for the bare class label.
    "classified_as_ccre": (CCRE,),
    # 3D-genome semantics (v0.2.0 §6 split of ``interacts_3d_with``).
    # Only ``compartment_of`` has a controlled label vocabulary; the
    # other 3D qualifiers (``contacts_with``, ``loops_with``,
    # ``anchored_to``, ``tad_boundary_with``, ``subtad_of``) take
    # coordinate / anchor / TAD-ID objects and rely on the default
    # cascade — preference here would only help once those objects
    # have a CURIE-resolving registry.
    "compartment_of": (CONFORMATION,),
}


def _lookup_source(source: Source, value: str) -> dict[str, str] | None:
    """Single dispatch from a source label to the underlying index lookup.

    Kept as a function (rather than a dict of bound methods) so test code
    can monkey-patch ``HISTONE_INDEX.lookup`` / ``TARGET_INDEX.lookup``
    individually without rebuilding the dispatch table.
    """
    if source == ACCESSION:
        return accession_to_curie(value)
    if source == ASSAY:
        return ASSAY_INDEX.lookup(value)
    if source == ASSEMBLY:
        return ASSEMBLY_INDEX.lookup(value)
    if source == HISTONE:
        return HISTONE_INDEX.lookup(value)
    if source == TARGET:
        return TARGET_INDEX.lookup(value)
    if source == QC:
        return QC_INDEX.lookup(value)
    if source == CHROM_STATE:
        return CHROM_STATE_INDEX.lookup(value)
    if source == DNA_MOD:
        return DNA_MOD_INDEX.lookup(value)
    if source == IHEC_BIOSAMPLE:
        return IHEC_BIOSAMPLE_INDEX.lookup(value)
    if source == CCRE:
        return CCRE_CLASS_INDEX.lookup(value)
    if source == CONFORMATION:
        return CONFORMATION_INDEX.lookup(value)
    return None


def resolve(
    value: str,
    qualifier: str | None = None,
) -> tuple[dict[str, str] | None, Source]:
    """Resolve ``value`` to a ``(match_dict, source_label)`` pair.

    Lookup order:
      1. ACCESSION — always first; regex-deterministic and outranks
         free-text matches under any qualifier.
      2. Qualifier-specific preferences (in declared order).
      3. Remaining sources in the default order — any source not already
         tried in (1) or (2).

    Returns ``(None, "")`` when nothing matches.

    Pure function — does not mutate the input string or any index.
    """
    # Step 1: accession always first.
    if match := _lookup_source(ACCESSION, value):
        return match, ACCESSION

    tried: set[Source] = {ACCESSION}

    # Step 2: qualifier-specific preferences.
    preferred = _QUALIFIER_PREFERENCES.get(qualifier or "", ())
    for source in preferred:
        if source in tried:
            continue
        tried.add(source)
        if match := _lookup_source(source, value):
            return match, source

    # Step 3: remaining default order.
    for source in _DEFAULT_ORDER:
        if source in tried:
            continue
        tried.add(source)
        if match := _lookup_source(source, value):
            return match, source

    return None, ""
