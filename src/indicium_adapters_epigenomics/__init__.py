"""indicium-adapters-epigenomics — epigenomics-domain adapters for the indicium ecosystem.

v0.1.0 ships three sibling adapters, all registered under the
``indicium.adapters`` entry-point group:

- ``BulkEpigenomicsAdapter`` — bulk epigenomics assay claims covering
  ChIP-seq / CUT&RUN / CUT&Tag (chromatin binding and histone modification),
  ATAC-seq / MNase-seq / DNase-seq (chromatin accessibility), WGBS / RRBS /
  EM-seq (DNA methylation), and Hi-C / Micro-C / ChIA-PET / HiChIP
  (chromosome conformation). Repository-accession enrichment (ENCODE / EpiRR
  / GEO / SRA / ENA) and chromatin-target resolution via a bundled MOD /
  UniProt index.
- ``SingleCellEpigenomicsAdapter`` — single-cell and single-nucleus
  epigenomics claims covering scATAC-seq / snATAC-seq, snmC-seq / sci-MET /
  snmCAT-seq, scChIP-seq / single-cell CUT&Tag, scNOMe-seq, and scHi-C.
  Cell-type / anatomy (CL + UBERON) enrichment lands in v0.2.0.
- ``MethylationArrayEpigenomicsAdapter`` — Illumina methylation array claims
  for the 450K, EPIC v1, EPIC v2, and MethylationEPIC mouse arrays.

The ``accession_to_curie`` helper from ``_registries`` is re-exported for
bridge layers that want to normalise ENCODE / EpiRR / GEO / SRA / ENA /
ArrayExpress / BioProject / EGA / ChIP-Atlas / Cistrome accessions before
emitting claims.
"""
from __future__ import annotations

from importlib.metadata import version

from indicium_adapters_epigenomics._registries import accession_to_curie
from indicium_adapters_epigenomics.array import MethylationArrayEpigenomicsAdapter
from indicium_adapters_epigenomics.bulk import BulkEpigenomicsAdapter
from indicium_adapters_epigenomics.singlecell import SingleCellEpigenomicsAdapter

__version__ = version("indicium-adapters-epigenomics")
__all__ = [
    "BulkEpigenomicsAdapter",
    "MethylationArrayEpigenomicsAdapter",
    "SingleCellEpigenomicsAdapter",
    "accession_to_curie",
]
