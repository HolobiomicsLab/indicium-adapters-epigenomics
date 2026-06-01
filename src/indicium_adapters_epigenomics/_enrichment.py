"""Shared enrichment back-end for epigenomics adapters.

This module is the single import point for the singleton indices used by the
three sibling adapters (``BulkEpigenomicsAdapter``,
``SingleCellEpigenomicsAdapter``, ``MethylationArrayEpigenomicsAdapter``).
Each index is constructed exactly once — at import time — and frozen via
its dataclass-style wrapper so the runtime side of the package is pure-
function (no I/O at ``enrich_claim()`` call time).

Index ownership map (adapters do not import the underlying modules; they
import the names re-exported here):

- ``HISTONE_INDEX``      — histone PTMs + variants + DNA-modification chemistries
                           (bundled TSV: ``data/histone_marks.tsv``)
- ``TARGET_INDEX``       — TFs + chromatin modifiers, human + mouse
                           (bundled TSV: ``data/chromatin_targets.tsv``)
- ``ASSEMBLY_INDEX``     — reference-genome assemblies (inline)
- ``ASSAY_INDEX``        — epigenomics assay-name → CURIE (inline)
- ``QC_INDEX``           — QC-metric name → CURIE (inline)
- ``CHROM_STATE_INDEX``  — ChromHMM / Segway chromatin-state labels
                           (bundled TSV: ``data/chromatin_states.tsv``,
                           added v0.2.0 §2)
- ``DNA_MOD_INDEX``      — DNA-base modification chemistries
                           (bundled TSV: ``data/dna_modifications.tsv``,
                           split out of ``HISTONE_INDEX`` in v0.2.0 §3)
- ``IHEC_BIOSAMPLE_INDEX`` — IHEC core biosample taxonomy with CL /
                           UBERON / EFO cross-refs
                           (bundled TSV: ``data/ihec_biosample_index.tsv``,
                           added v0.2.0 §4)
- ``CCRE_CLASS_INDEX``   — ENCODE SCREEN cCRE class labels
                           (bundled TSV: ``data/ccre_classes.tsv``,
                           added v0.2.0 §5)
- ``CONFORMATION_INDEX`` — 3D-genome compartment / sub-compartment labels
                           (bundled TSV: ``data/conformation_labels.tsv``,
                           added v0.2.0 §6)
"""
from __future__ import annotations

from typing import Final

from indicium_adapters_epigenomics._assays import ASSAY_INDEX, AssayIndex
from indicium_adapters_epigenomics._assemblies import ASSEMBLY_INDEX, AssemblyIndex
from indicium_adapters_epigenomics._biosamples import (
    IHEC_BIOSAMPLE_INDEX,
    IHECBiosampleIndex,
)
from indicium_adapters_epigenomics._ccre import CCRE_CLASS_INDEX, CCREClassIndex
from indicium_adapters_epigenomics._chrom_states import (
    CHROM_STATE_INDEX,
    ChromStateIndex,
)
from indicium_adapters_epigenomics._conformation import (
    CONFORMATION_INDEX,
    ConformationIndex,
)
from indicium_adapters_epigenomics._dna_mods import DNA_MOD_INDEX, DNAModIndex
from indicium_adapters_epigenomics._histone_marks import HISTONE_INDEX, HistoneIndex
from indicium_adapters_epigenomics._qc_metrics import QC_INDEX, QCIndex
from indicium_adapters_epigenomics._targets import TARGET_INDEX, TargetIndex

__all__ = [
    "ASSAY_INDEX", "AssayIndex",
    "ASSEMBLY_INDEX", "AssemblyIndex",
    "CCRE_CLASS_INDEX", "CCREClassIndex",
    "CHROM_STATE_INDEX", "ChromStateIndex",
    "CONFORMATION_INDEX", "ConformationIndex",
    "DNA_MOD_INDEX", "DNAModIndex",
    "HISTONE_INDEX", "HistoneIndex",
    "IHEC_BIOSAMPLE_INDEX", "IHECBiosampleIndex",
    "QC_INDEX", "QCIndex",
    "TARGET_INDEX", "TargetIndex",
]

# Re-typed at module scope so static analysers see the public surface.
ASSAY_INDEX: Final[AssayIndex] = ASSAY_INDEX
ASSEMBLY_INDEX: Final[AssemblyIndex] = ASSEMBLY_INDEX
CCRE_CLASS_INDEX: Final[CCREClassIndex] = CCRE_CLASS_INDEX
CHROM_STATE_INDEX: Final[ChromStateIndex] = CHROM_STATE_INDEX
CONFORMATION_INDEX: Final[ConformationIndex] = CONFORMATION_INDEX
DNA_MOD_INDEX: Final[DNAModIndex] = DNA_MOD_INDEX
HISTONE_INDEX: Final[HistoneIndex] = HISTONE_INDEX
IHEC_BIOSAMPLE_INDEX: Final[IHECBiosampleIndex] = IHEC_BIOSAMPLE_INDEX
QC_INDEX: Final[QCIndex] = QC_INDEX
TARGET_INDEX: Final[TargetIndex] = TARGET_INDEX
