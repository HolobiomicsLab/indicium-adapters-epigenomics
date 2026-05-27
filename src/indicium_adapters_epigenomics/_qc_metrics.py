"""Epigenomics QC-metric name → CURIE resolver.

Pure functions, no I/O. Kept inline as a Python literal (the QC vocabulary
is small and curated, mirroring transcriptomics' inline style for small
registries).

CURIE choice policy: EDAM operation IRIs where a precise term exists,
otherwise a local ``EPIQC:`` vocabulary term — these are marked in
``_namespaces.py`` for later promotion to STATO / EDAM when upstream terms
land.
"""
from __future__ import annotations

import re
from types import MappingProxyType
from typing import Final, NamedTuple


class _QCSpec(NamedTuple):
    """One QC-metric table row."""

    pattern: re.Pattern[str]
    curie: str
    label: str
    applies_to_assay_class: str  # bulk_sequencing | single_cell | array | any


_QC_TABLE: Final[tuple[_QCSpec, ...]] = (
    # --- ChIP-seq / CUT&RUN / CUT&Tag --------------------------------------
    _QCSpec(
        re.compile(r"\bFRiP\b", re.I),
        "EPIQC:frip", "fraction of reads in peaks", "bulk_sequencing",
    ),
    _QCSpec(
        re.compile(
            r"\b(?:NSC|normalized\s+strand\s+cross[- ]?correlation)\b",
            re.I,
        ),
        "EPIQC:nsc", "normalized strand cross-correlation", "bulk_sequencing",
    ),
    _QCSpec(
        re.compile(
            r"\b(?:RSC|relative\s+strand\s+cross[- ]?correlation)\b",
            re.I,
        ),
        "EPIQC:rsc", "relative strand cross-correlation", "bulk_sequencing",
    ),
    _QCSpec(
        re.compile(r"\bIDR(?:[- ]threshold)?\b", re.I),
        "EPIQC:idr-threshold", "irreproducible discovery rate threshold",
        "bulk_sequencing",
    ),
    # --- ATAC / DNase / accessibility --------------------------------------
    _QCSpec(
        re.compile(r"\bTSS[- ]?enrichment\b", re.I),
        "EPIQC:tss-enrichment", "TSS enrichment", "bulk_sequencing",
    ),
    _QCSpec(
        re.compile(r"\bfragment[- ]size\s+periodicity\b", re.I),
        "EPIQC:fragment-size-periodicity", "fragment size periodicity",
        "bulk_sequencing",
    ),
    # --- Library complexity (generic sequencing) ---------------------------
    _QCSpec(
        re.compile(r"\bPBC1\b", re.I),
        "EPIQC:pbc1", "PCR bottleneck coefficient 1", "bulk_sequencing",
    ),
    _QCSpec(
        re.compile(r"\bPBC2\b", re.I),
        "EPIQC:pbc2", "PCR bottleneck coefficient 2", "bulk_sequencing",
    ),
    _QCSpec(
        re.compile(r"\bNRF\b", re.I),
        "EPIQC:nrf", "non-redundant fraction", "bulk_sequencing",
    ),
    _QCSpec(
        re.compile(r"\b(?:library[- ]complexity)\b", re.I),
        "EPIQC:library-complexity", "library complexity", "any",
    ),
    _QCSpec(
        re.compile(r"\bduplication[- ]rate\b", re.I),
        "EPIQC:duplication-rate", "duplication rate", "any",
    ),
    _QCSpec(
        re.compile(r"\bmapping[- ]rate\b", re.I),
        "EDAM:operation_3198", "mapping rate", "any",
    ),
    # --- Peak / region metrics --------------------------------------------
    _QCSpec(
        re.compile(r"\bpeak[- ]count\b", re.I),
        "EPIQC:peak-count", "peak count", "bulk_sequencing",
    ),
    _QCSpec(
        re.compile(r"\bpeak[- ]width(?:[- ]distribution)?\b", re.I),
        "EPIQC:peak-width-distribution", "peak width distribution",
        "bulk_sequencing",
    ),
    # --- WGBS / methylation ------------------------------------------------
    _QCSpec(
        re.compile(
            r"\bbisulfite[- ]conversion(?:[- ](?:rate|efficiency))?\b",
            re.I,
        ),
        "EPIQC:bisulfite-conversion-rate", "bisulfite conversion rate",
        "bulk_sequencing",
    ),
    # --- Methylation array -------------------------------------------------
    _QCSpec(
        re.compile(r"\bdetection[ -]?p[- ]value\b", re.I),
        "EPIQC:detection-p-value", "detection p-value", "array",
    ),
    _QCSpec(
        re.compile(r"\bbead[ -]?count\b", re.I),
        "EPIQC:bead-count", "bead count", "array",
    ),
    _QCSpec(
        re.compile(r"\bbeta[ -]?value(?:[ -]distribution)?\b", re.I),
        "EPIQC:beta-value-distribution", "beta-value distribution", "array",
    ),
    _QCSpec(
        re.compile(r"\bcross[- ]reactive(?:\s+probe(?:s)?)?(?:[- ]filter)?\b", re.I),
        "EPIQC:cross-reactive-probe-filter", "cross-reactive probe filter",
        "array",
    ),
    _QCSpec(
        re.compile(
            r"\b(?:cell\s+composition|cell[- ]type\s+composition|"
            r"Houseman(?:\s+reference)?)\b",
            re.I,
        ),
        "EPIQC:houseman-reference", "Houseman cell-composition reference",
        "array",
    ),
    _QCSpec(
        re.compile(r"\bfunnorm(?:[ -]PCs)?\b", re.I),
        "EPIQC:funnorm-pcs", "functional normalization PCs", "array",
    ),
    _QCSpec(
        re.compile(r"\bBMIQ(?:[ -]alpha)?\b", re.I),
        "EPIQC:bmiq-alpha", "BMIQ alpha parameter", "array",
    ),
    _QCSpec(
        re.compile(r"\bNOOB(?:[ -]offset)?\b", re.I),
        "EPIQC:noob-offset", "NOOB background offset", "array",
    ),
)


class QCIndex:
    """Case-insensitive QC-metric name → CURIE lookup."""

    def __init__(self, table: tuple[_QCSpec, ...]):
        self._table = table

    def lookup(self, text: str) -> dict[str, str] | None:
        """Return the matching QC-metric row or ``None``.

        The first regex in the precedence-ordered table to match wins. Match
        is a ``re.search`` (not anchored) — pass a bare metric name or a
        snippet of methods text.
        """
        if not isinstance(text, str) or not text:
            return None
        for spec in self._table:
            if spec.pattern.search(text):
                return {
                    "curie": spec.curie,
                    "label": spec.label,
                    "applies_to_assay_class": spec.applies_to_assay_class,
                }
        return None

    def __len__(self) -> int:
        return len(self._table)


QC_INDEX: Final[QCIndex] = QCIndex(_QC_TABLE)

QC_LABEL_TO_CURIE: Final[MappingProxyType[str, str]] = MappingProxyType(
    {spec.label: spec.curie for spec in _QC_TABLE}
)
