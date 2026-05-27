# indicium-adapters-epigenomics

Epigenomics-domain DomainAdapters for the [indicium](https://github.com/holobiomicslab/indicium) ecosystem.

v0.1.0 ships three sibling adapters:

- **`BulkEpigenomicsAdapter`** — bulk epigenomics assay claims covering
  WGBS / RRBS / EM-seq (DNA methylation), ChIP-seq / CUT&RUN / CUT&Tag
  (chromatin-binding and histone-modification), ATAC-seq / MNase-seq / DNase-seq
  (chromatin-accessibility), and chromosome-conformation capture assays
  (Hi-C / Micro-C / ChIA-PET / HiChIP). Resolves assay names to OBI / EFO
  CURIEs, repository accessions (ENCODE / EpiRR / GEO / SRA / ENA /
  BioProject / ArrayExpress / EGA / ChIP-Atlas / Cistrome), histone marks,
  chromatin-target proteins (107 human + mouse, UniProt), genome assemblies
  (RefSeq), and bulk-sequencing QC metrics.
- **`SingleCellEpigenomicsAdapter`** — single-cell and single-nucleus
  epigenomics claims covering scATAC-seq / snATAC-seq, snmC-seq / sci-MET /
  snmCAT-seq (single-cell methylation), scChIP-seq / single-cell CUT&Tag,
  scNOMe-seq, and scHi-C. Resolves assay names, histone marks, chromatin
  targets, and the same accession / assembly registries as bulk; CL / UBERON
  cell-type and tissue enrichment is on the v0.2.0 roadmap.
- **`MethylationArrayEpigenomicsAdapter`** — Illumina methylation array claims
  covering the 450K, EPIC v1, EPIC v2, and MethylationEPIC mouse arrays.
  Resolves array-platform names to verified OBI CURIE (`OBI:0002131` for
  MethylationEPIC v1) or local `EpiAssay:` CURIEs for the others, and detects
  Illumina ``cg########`` probe IDs to store in `xrefs` (no `ontology_terms`
  entry — by design, no upstream CURIE registry for probes).

All three register under the `indicium.adapters` entry-point group and are
discoverable via `indicium_adapters.discover_adapters()`. Their qualifier sets
are disjoint by construction (enforced by `tests/test_composition.py`).

## Install

```bash
pip install indicium-adapters-epigenomics
```

`indicium` and `indicium-adapters` are not on PyPI yet; if installing from
source, pull them from git first:

```bash
pip install git+https://github.com/HolobiomicsLab/indicium.git@main
pip install git+https://github.com/HolobiomicsLab/indicium-adapters.git@main
pip install -e .[dev]
```

## Usage

```python
from indicium_adapters import compose_adapters, discover_adapters

found = discover_adapters()
bulk   = found["bulk-epigenomics"]
sc     = found["single-cell-epigenomics"]
array  = found["methylation-array-epigenomics"]

# 1) Repository accession → CURIE
enriched = bulk.enrich_claim({
    "subject": "ENCSR000AKB",
    "qualifier": "consolidates_dataset",
    "object": "ENCODE H3K27ac ChIP-seq",
})
# enriched["ontology_terms"]["subject"] == "ENCODE:ENCSR000AKB"

# 2) Histone mark → EpiHist CURIE (v0.1.0 local vocab; will move to MOD/PR
#    when upstream coverage of position-specific histone PTMs improves)
enriched = bulk.enrich_claim({
    "subject": "H3K27ac",
    "qualifier": "marked_by",
    "object": "H3K27ac",
})
# enriched["ontology_terms"]["object"] == "EpiHist:H3K27ac"

# 3) Single-cell anatomy → UBERON CURIE
enriched = sc.enrich_claim({
    "subject": "ATAC-peak-001",
    "qualifier": "accessible_in_cell_type",
    "object": "gastroderm",
})

# 4) Methylation array probe → CpG annotation CURIE
enriched = array.enrich_claim({
    "subject": "cg00000029",
    "qualifier": "methylated_at_locus",
    "object": "TAP2 gene body",
})

# Use all three adapters at once via the composition helper
composite = compose_adapters([bulk, sc, array])
```

## SHACL shapes

All three adapters' `shacl_shapes()` return an `rdflib.Graph` with a NodeShape
that validates CURIE format on ASB `Claim` graphs. Pipeline-level RDF validation
for epigenomics workflow classes (SO / GO / UBERON / PR / MOD shapes) lives
upstream in the consuming graph store, not here.

## Data provenance — what's live-verified and what is local vocabulary

Every CURIE this package emits falls into one of two buckets. **Live-verified
upstream CURIEs** assert that a recognised public database actually holds that
ID with a matching label; each was cross-checked against the relevant REST
API. **Local-namespace CURIEs** are minted by this package for concepts
without a stable upstream term and are clearly prefixed so consumers do not
mistake them for authoritative IDs.

| Bucket | Live source | Coverage in v0.1.0 |
| --- | --- | --- |
| `UniProt:` | UniProt REST | 107 chromatin targets + 5 histone variants |
| `RefSeq:` | NCBI Datasets v2 + Entrez eutils | 7 human + mouse assemblies |
| `CHEBI:` | EBI ChEBI + OLS4 | 3 DNA modifications (5mC, 5hmC, 5fC) |
| `OBI:` | EBI OLS4 | 9 assays (ChIP-seq, DNase-seq, MNase-seq, WGBS, FAIRE-seq, ChIA-PET, snATAC-seq, MethylationEPIC v1) |
| `EFO:` | EBI OLS4 | 8 assays (ATAC-seq, scATAC-seq, Hi-C, CUT&RUN, CUT&Tag, HiChIP, RRBS, snmC-seq, scHi-C) |
| `EDAM:` | EBI OLS4 | 1 QC operation |
| `GEO:` / `SRA:` / `ENCODE:` / `EpiRR:` / `EGA:` / `BioProject:` / `ArrayExpress:` / `ENA:` / `ChIPAtlas:` / `Cistrome:` | identifiers.org schemes; accession regex recognition | Accession-syntax validation only |
| `EpiHist:` (local) | this package | 38 position-specific histone PTMs (H3K27ac etc.) + 2 DNA mods (5caC, 6mA) where no stable upstream term exists |
| `EpiAssay:` (local) | this package | 12 assays without a verified OBI/EFO term (HumanMethylation450, MethylationEPIC v2, Mouse-MethylationEPIC, Capture Hi-C, Micro-C, EM-seq, oxBS-seq, snmCAT-seq, sci-MET, scChIP-seq, scNOMe-seq, single-cell CUT&Tag) |
| `EPIQC:` (local) | this package | ~22 quality-control metric names (FRiP, NSC, RSC, TSS-enrichment, bisulfite-conversion-rate, etc.) — these are metric names, not ontology terms |

Local-namespace CURIEs are documented in `_namespaces.py` and described in
`_histone_marks.py` / `_assays.py` / `_qc_metrics.py`. They are slated for
promotion to the corresponding upstream prefix when upstream coverage of the
concept becomes available — see `ROADMAP.md` for the upgrade path.

Two large bundled snapshots are refreshed by maintainer scripts:

- `data/encode_experiment_index.tsv` — 5,000 epigenomics-only experiments
  from the live ENCODE portal REST API
- `data/epirr_reference_index.tsv` — 325 IHEC reference epigenomes (ENCODE +
  Roadmap) from the ENCODE `ReferenceEpigenome` REST endpoint (the EBI EpiRR
  JSON API was retired in 2024–2025; we federate via ENCODE)

## Development

```bash
git clone https://github.com/holobiomicslab/indicium-adapters-epigenomics
cd indicium-adapters-epigenomics
pip install -e .[dev]
pytest
```

Dev commands:

```bash
ruff check .          # lint
pytest --cov          # test + coverage (must reach 90%)
```

See `CLAUDE.md` for the dev guide and the 4-quadrant filter used to decide
whether a feature belongs in this adapter at all.

## License

MIT — see `LICENSE`.
