"""Shared enrichment back-end for epigenomics adapters.

This module is the single import point for the singleton indices used by the
three sibling adapters (``BulkEpigenomicsAdapter``,
``SingleCellEpigenomicsAdapter``, ``MethylationArrayEpigenomicsAdapter``).
Each index is constructed exactly once — at import time — and frozen via
its dataclass-style wrapper so the runtime side of the package is pure-
function (no I/O at ``enrich_claim()`` call time).

Index ownership map (adapters do not import the underlying modules; they
import the names re-exported here):

- ``HISTONE_INDEX``   — histone PTMs + variants + DNA-modification chemistries
                        (bundled TSV: ``data/histone_marks.tsv``)
- ``TARGET_INDEX``    — TFs + chromatin modifiers, human + mouse
                        (bundled TSV: ``data/chromatin_targets.tsv``)
- ``ASSEMBLY_INDEX``  — reference-genome assemblies (inline)
- ``ASSAY_INDEX``     — epigenomics assay-name → CURIE (inline)
- ``QC_INDEX``        — QC-metric name → CURIE (inline)
"""
from __future__ import annotations

from typing import Final

from indicium_adapters_epigenomics._assays import ASSAY_INDEX, AssayIndex
from indicium_adapters_epigenomics._assemblies import ASSEMBLY_INDEX, AssemblyIndex
from indicium_adapters_epigenomics._histone_marks import HISTONE_INDEX, HistoneIndex
from indicium_adapters_epigenomics._qc_metrics import QC_INDEX, QCIndex
from indicium_adapters_epigenomics._targets import TARGET_INDEX, TargetIndex

__all__ = [
    "ASSAY_INDEX", "AssayIndex",
    "ASSEMBLY_INDEX", "AssemblyIndex",
    "HISTONE_INDEX", "HistoneIndex",
    "QC_INDEX", "QCIndex",
    "TARGET_INDEX", "TargetIndex",
]

# Re-typed at module scope so static analysers see the public surface.
ASSAY_INDEX: Final[AssayIndex] = ASSAY_INDEX
ASSEMBLY_INDEX: Final[AssemblyIndex] = ASSEMBLY_INDEX
HISTONE_INDEX: Final[HistoneIndex] = HISTONE_INDEX
QC_INDEX: Final[QCIndex] = QC_INDEX
TARGET_INDEX: Final[TargetIndex] = TARGET_INDEX
