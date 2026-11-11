# CLAUDE.md — Dev Guide

## Read this before adding anything

The package's design philosophy and the 4-quadrant filter for new features
live at the top of `ROADMAP.md`. Any new contribution should justify itself
in those terms: does it **constrain** the agent, **resolve** identifiers,
**validate** output, or is it something an agent could do as well freehand?
If the last, skip it.

The adapter is **not** a pipeline orchestrator and does **not** score peaks,
call differentially methylated regions, or run alignment steps. Those belong
upstream in the ASB pipeline before claims are emitted.

## Project shape

Three sibling `DomainAdapter` classes in one wheel:

- `BulkEpigenomicsAdapter` (`bulk.py`) — bulk epigenomics assay claims and
  study-level assay metadata for ChIP-seq / ATAC-seq / WGBS / Hi-C and
  related protocols.
- `SingleCellEpigenomicsAdapter` (`singlecell.py`) — cell-resolution
  epigenomic claims (scATAC / snATAC / snmC-seq / scChIP / scHi-C). Resolves
  assay names, histone marks, chromatin targets, accessions, and assemblies
  in v0.1.0; CL / UBERON cell-type and tissue enrichment is on the v0.2.0
  roadmap.
- `MethylationArrayEpigenomicsAdapter` (`array.py`) — Illumina methylation
  array claims (450K / EPIC v1 / EPIC v2 / MethylationEPIC mouse). Resolves
  array-platform names to OBI / `EpiAssay:` CURIEs and detects ``cg########``
  Illumina probe IDs (stored in `xrefs` only — no upstream CURIE registry).

All three target:

- `api_version = 1` matching `indicium-adapters.CURRENT_API_VERSION`.
- `DomainAdapter` and `SHACLProvider` Protocol conformance.
- Entry-point registration under `indicium.adapters`.
- Hardcoded namespace URI constants in `_namespaces.py`.
- **Disjoint qualifier sets** — enforced by `tests/test_composition.py`.

## Design constraints — what does NOT belong here

- **Pipeline orchestration** — no peak calling, no DMR calling, no alignment.
  These belong upstream in ASB before claims are emitted.
- **Scoring or ranking** — no enrichment scores, no motif scoring. An LLM
  extractor can attach numbers as claim attributes without adapter help.
- **Redundant lookup** — if a resolver already exists in `indicium-adapters`
  core, do not re-implement it here.

## 4-quadrant filter for new code

Before adding any function, method, or module, answer:

1. Does it **constrain** the agent (qualifier set, extraction-context wording)?
2. Does it **resolve** identifiers (histone mark → `EpiHist:` (local) /
   `MOD:` / `PR:` CURIE, ENCODE accession → CURIE, assembly name → RefSeq /
   UCSC identifier)?
3. Does it **validate** output (SHACL NodeShape, CURIE-format check)?
4. Or could an agent do it as well freehand?

If only (4) applies, **skip it**.

## 7-step checklist for adding a qualifier

1. Decide which adapter owns it (or whether it justifies a new sibling).
2. Add it to that adapter's `qualifiers` frozenset.
3. Update `extraction_context()` to mention it.
4. If it needs enrichment behaviour, add a branch in `enrich_claim()`.
5. If it needs validation, add a NodeShape in the adapter's shapes TTL.
6. Add unit tests in `tests/test_<adapter>.py`.
7. Keep qualifier sets disjoint across all three sibling adapters.

## CURIE provenance policy

Two buckets:

1. **Live-verified upstream CURIEs.** Every UniProt / RefSeq / CHEBI / OBI /
   EFO / EDAM CURIE in this package has been cross-checked against the
   relevant REST API (UniProt REST, NCBI Datasets v2, EBI ChEBI, OLS4) with
   a matching label. If you add a new one, run the same verification before
   committing.
2. **Local-namespace CURIEs** (`EpiHist:`, `EpiAssay:`, `EPIQC:`) for
   concepts without a stable upstream term — position-specific histone PTMs
   (PR has only organism + variant proteoforms; MOD doesn't go to
   position-level), a handful of single-cell methylation assays, and most
   QC metric names. These are minted under `https://w3id.org/holobiomicslab/...`
   IRIs declared in `_namespaces.py` and explicitly documented in the
   relevant `_*.py` module docstrings as v0.1.0 local vocabulary scheduled
   for promotion when upstream coverage improves.

**Never assert an upstream CURIE without live verification.** The earlier
draft of this package shipped 39 fabricated `MOD:` CURIEs (e.g.
`MOD:00581 — H3K27 acetylation`, where `MOD:00581` does not exist in PSI-MOD
at all). Use the local namespace instead — it's honest about provenance.

## Refreshing bundled data files

Maintainer-only scripts under `src/indicium_adapters_epigenomics/data/`:

- `build_encode_index.py` — regenerates `encode_experiment_index.tsv` from
  the ENCODE portal REST API, filtered to the epigenomics assay families
  the adapter resolves. Scheduled as a monthly CI cron job from v0.2.0.
- `build_epirr_index.py` — regenerates `epirr_reference_index.tsv` from the
  ENCODE `ReferenceEpigenome` REST endpoint. The EBI EpiRR JSON API was
  retired in 2024–2025; we federate via ENCODE and key on the IHEC
  `IHECRE…` accession found in `dbxrefs`.
- `build_mod_index.py` — pulls the PSI-MOD descendants of MOD:01156 from
  OLS4 to a staging file. NOT a direct overwrite of `histone_marks.tsv`;
  the maintainer compares against the staging file when promoting an
  `EpiHist:` entry to a verified `MOD:` / `PR:` CURIE.

**Do NOT run any of these at install or import time.**

## Test before commit

```bash
ruff check .
pytest --cov=indicium_adapters_epigenomics
```

CI matrix is 3.12 / 3.13.

## Repo relationships

- Core Protocol lives in the `indicium-adapters` sibling repository (do not vendor it).
- The indicium claim/evidence vocabulary lives in the `indicium` sibling repository.
- Pipeline RDF validation for epigenomics workflow classes (SO / GO /
  UBERON / PR / MOD / EFO node shapes) belongs upstream in the consuming
  graph store — we do not duplicate it here.
