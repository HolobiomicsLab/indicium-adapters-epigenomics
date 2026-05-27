"""Epigenomics assay name → CURIE resolver.

Pure functions, no I/O. The assay table is kept inline as a Python literal
(mirroring transcriptomics' style for small, hand-curated registries) rather
than as a bundled TSV — the canonical assay set rarely changes, and inline
data keeps the resolver one ``import`` away.

CURIE choice policy in v0.1.0:

* **Verified OBI / EFO CURIEs** (cross-checked against EBI OLS4 on
  2026-05-27) are kept under their canonical prefix.
* **Local ``EpiAssay:`` CURIEs** are used everywhere else — many
  community-popular epigenomics assays (CUT&RUN, CUT&Tag, MethylationEPIC,
  EM-seq, oxBS-seq, single-cell methylation protocols, HiChIP, Micro-C,
  Capture Hi-C) had no verified term in OBI / EFO at the time of writing.
  When upstream coverage improves these promote to the canonical prefix
  and the local namespace is retired (see ROADMAP.md).

Single-cell variants of an assay are kept as separate rows (rather than
inferred from a prefix) so the ``assay_class`` discriminator is explicit;
this is what the three sibling adapters consume to gate which qualifiers
to apply.
"""
from __future__ import annotations

import re
from types import MappingProxyType
from typing import Final, NamedTuple


class _AssaySpec(NamedTuple):
    """One assay table row."""

    pattern: re.Pattern[str]
    curie: str
    label: str
    assay_class: str  # bulk_sequencing | single_cell | array


# Order matters: more-specific assays before more-general ones. The lookup
# walks top-down and returns the first match, so ``scATAC-seq`` must precede
# ``ATAC-seq`` and ``snmCAT-seq`` must precede the underlying ``snmC-seq``.
_ASSAY_TABLE: Final[tuple[_AssaySpec, ...]] = (
    # --- Single-cell variants (must precede their bulk counterparts) -------
    # scATAC-seq: EFO:0010891 verified ("scATAC-seq").
    _AssaySpec(
        re.compile(r"\bsc[- ]?ATAC[- ]?seq\b", re.I),
        "EFO:0010891", "scATAC-seq", "single_cell",
    ),
    # snATAC-seq: OBI:0002762 verified ("single-nucleus ATAC-seq").
    _AssaySpec(
        re.compile(r"\bsn[- ]?ATAC[- ]?seq\b", re.I),
        "OBI:0002762", "snATAC-seq", "single_cell",
    ),
    # Single-cell methylation / co-assays — no verified OBI / EFO; local.
    _AssaySpec(
        re.compile(r"\bsnmCAT[- ]?seq\b", re.I),
        "EpiAssay:snmCAT-seq", "snmCAT-seq", "single_cell",
    ),
    # snmC-seq: EFO:0008939 verified.
    _AssaySpec(
        re.compile(r"\bsnmC[- ]?seq2?\b", re.I),
        "EFO:0008939", "snmC-seq", "single_cell",
    ),
    _AssaySpec(
        re.compile(r"\bsci[- ]?MET\b", re.I),
        "EpiAssay:sci-MET", "sci-MET", "single_cell",
    ),
    _AssaySpec(
        re.compile(r"\bscChIP[- ]?seq\b", re.I),
        "EpiAssay:scChIP-seq", "scChIP-seq", "single_cell",
    ),
    _AssaySpec(
        re.compile(r"\bsingle[- ]cell\s+CUT[ &]?Tag\b", re.I),
        "EpiAssay:single-cell-CUT-Tag", "single-cell CUT&Tag", "single_cell",
    ),
    _AssaySpec(
        re.compile(r"\bscNOMe[- ]?seq\b", re.I),
        "EpiAssay:scNOMe-seq", "scNOMe-seq", "single_cell",
    ),
    # scHi-C: EFO:0009979 verified ("single cell Hi-C").
    _AssaySpec(
        re.compile(r"\bscHi[- ]?C\b", re.I),
        "EFO:0009979", "scHi-C", "single_cell",
    ),
    # --- DNA methylation (bulk sequencing) ---------------------------------
    # WGBS: OBI:0001863 verified ("shotgun bisulfite-seq assay").
    _AssaySpec(
        re.compile(
            r"\b(?:whole[- ]genome\s+bisulfite\s+sequencing|WGBS|"
            r"bs[- ]seq|methylC[- ]seq)\b",
            re.I,
        ),
        "OBI:0001863", "whole-genome bisulfite sequencing", "bulk_sequencing",
    ),
    # RRBS: EFO:0008900 verified ("RRBS").
    _AssaySpec(
        re.compile(
            r"\b(?:reduced\s+representation\s+bisulfite\s+sequencing|RRBS)\b",
            re.I,
        ),
        "EFO:0008900", "reduced representation bisulfite sequencing",
        "bulk_sequencing",
    ),
    _AssaySpec(
        re.compile(r"\bEM[- ]?seq\b", re.I),
        "EpiAssay:EM-seq", "EM-seq", "bulk_sequencing",
    ),
    _AssaySpec(
        re.compile(r"\boxBS[- ]?seq\b", re.I),
        "EpiAssay:oxBS-seq", "oxBS-seq", "bulk_sequencing",
    ),
    # --- Targeted IP (bulk sequencing) -------------------------------------
    # CUT&RUN: EFO:0009973 verified.
    _AssaySpec(
        re.compile(r"\bCUT[ &]?RUN\b", re.I),
        "EFO:0009973", "CUT&RUN", "bulk_sequencing",
    ),
    # CUT&Tag: EFO:0920004 verified.
    _AssaySpec(
        re.compile(r"\bCUT[ &]?Tag\b", re.I),
        "EFO:0920004", "CUT&Tag", "bulk_sequencing",
    ),
    # ChIP-seq: OBI:0000716 verified ("ChIP-seq assay").
    _AssaySpec(
        re.compile(r"\bChIP[- ]?seq\b", re.I),
        "OBI:0000716", "ChIP-seq", "bulk_sequencing",
    ),
    # --- Chromatin accessibility (bulk sequencing) -------------------------
    # ATAC-seq: EFO:0007045 verified ("ATAC-seq").
    _AssaySpec(
        re.compile(r"\bATAC[- ]?seq\b", re.I),
        "EFO:0007045", "ATAC-seq", "bulk_sequencing",
    ),
    # DNase-seq: OBI:0001853 verified ("DNase I hypersensitive sites …").
    _AssaySpec(
        re.compile(r"\bDNase[- ]?seq\b", re.I),
        "OBI:0001853", "DNase-seq", "bulk_sequencing",
    ),
    # MNase-seq: OBI:0001924 verified ("micrococcal nuclease digestion …").
    _AssaySpec(
        re.compile(r"\bMNase[- ]?seq\b", re.I),
        "OBI:0001924", "MNase-seq", "bulk_sequencing",
    ),
    # FAIRE-seq: OBI:0001859 verified ("FAIRE assay").
    _AssaySpec(
        re.compile(r"\bFAIRE[- ]?seq\b", re.I),
        "OBI:0001859", "FAIRE-seq", "bulk_sequencing",
    ),
    # --- 3D / conformation (bulk sequencing) -------------------------------
    _AssaySpec(
        re.compile(r"\bcapture\s+Hi[- ]?C\b", re.I),
        "EpiAssay:Capture-Hi-C", "Capture Hi-C", "bulk_sequencing",
    ),
    # HiChIP: EFO:0010011 verified.
    _AssaySpec(
        re.compile(r"\bHiChIP\b", re.I),
        "EFO:0010011", "HiChIP", "bulk_sequencing",
    ),
    # ChIA-PET: OBI:0001848 verified ("chromatin interaction analysis by
    # paired-end tag sequencing").
    _AssaySpec(
        re.compile(r"\bChIA[- ]?PET\b", re.I),
        "OBI:0001848", "ChIA-PET", "bulk_sequencing",
    ),
    _AssaySpec(
        re.compile(r"\bMicro[- ]?C\b", re.I),
        "EpiAssay:Micro-C", "Micro-C", "bulk_sequencing",
    ),
    # Hi-C: EFO:0007693 verified ("Hi-C").
    _AssaySpec(
        re.compile(r"\bHi[- ]?C\b", re.I),
        "EFO:0007693", "Hi-C", "bulk_sequencing",
    ),
    # --- Methylation arrays ------------------------------------------------
    # No verified EFO/OBI CURIEs as of 2026-05-27 — local namespace.
    _AssaySpec(
        re.compile(
            r"\b(?:HumanMethylation\s*450|Infinium\s*450K?|450K(?:\s+BeadChip)?)\b",
            re.I,
        ),
        "EpiAssay:HumanMethylation450",
        "Illumina HumanMethylation450 BeadChip", "array",
    ),
    _AssaySpec(
        re.compile(
            r"\b(?:EPIC\s*v?2|MethylationEPIC\s*v?2|Infinium\s*EPIC\s*v?2)\b",
            re.I,
        ),
        "EpiAssay:MethylationEPIC-v2",
        "Illumina MethylationEPIC v2.0 BeadChip", "array",
    ),
    # MethylationEPIC (v1): OBI:0002131 verified
    # ("Illumina Infinium MethylationEPIC v1.0 BeadChip").
    _AssaySpec(
        re.compile(
            r"\b(?:HumanMethylationEPIC|MethylationEPIC|EPIC\s*850K?|"
            r"Infinium\s*EPIC)\b",
            re.I,
        ),
        "OBI:0002131",
        "Illumina MethylationEPIC BeadChip", "array",
    ),
    _AssaySpec(
        re.compile(r"\bMouse\s+MethylationEPIC\b", re.I),
        "EpiAssay:Mouse-MethylationEPIC",
        "Illumina Mouse Methylation BeadChip", "array",
    ),
)


class AssayIndex:
    """Case-insensitive epigenomics-assay name → CURIE lookup.

    Pure data; the regex table is compiled at import time and frozen. The
    class only exposes ``lookup()`` so adapters stay decoupled from the
    on-disk format and the ordering rules.
    """

    def __init__(self, table: tuple[_AssaySpec, ...]):
        self._table = table

    def lookup(self, text: str) -> dict[str, str] | None:
        """Return the matching assay row as a frozen dict, or ``None``.

        The first regex in the precedence-ordered table to match wins. Match
        is a ``re.search`` (not anchored): the caller can pass either a bare
        assay name (``"scATAC-seq"``) or a snippet of methods text.
        """
        if not isinstance(text, str) or not text:
            return None
        for spec in self._table:
            if spec.pattern.search(text):
                return {
                    "curie": spec.curie,
                    "label": spec.label,
                    "assay_class": spec.assay_class,
                }
        return None

    def __len__(self) -> int:
        return len(self._table)


ASSAY_INDEX: Final[AssayIndex] = AssayIndex(_ASSAY_TABLE)

# Public view of the (label → curie) mapping for callers that want to iterate
# the supported assays (e.g. tests, extraction_context() rendering). Read-only.
ASSAY_LABEL_TO_CURIE: Final[MappingProxyType[str, str]] = MappingProxyType(
    {spec.label: spec.curie for spec in _ASSAY_TABLE}
)
