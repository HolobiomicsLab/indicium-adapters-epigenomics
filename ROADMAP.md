# Roadmap

## Design philosophy — what belongs here

Before adding anything, ask the 4-quadrant question:

1. **Constrain** the agent (qualifier set, prompt-context wording)?
2. **Resolve** identifiers (histone mark → `EpiHist:` (local) / `MOD:` /
   `PR:` CURIE, ENCODE / EpiRR accession → CURIE, genome assembly name →
   RefSeq / UCSC identifier, chromatin target → UniProt)?
3. **Validate** output (SHACL shape, CURIE-format check)?
4. Or is it something an agent could do as well freehand?

If only the last applies, **skip it**. The adapter is a constraints +
identifier-resolution + validation layer — NOT a pipeline orchestrator,
peak caller, DMR caller, or scorer. Those belong upstream in ASB before
claims are emitted.

The epigenomics adapter spans a wider assay landscape than transcriptomics;
the strategy is to ship *protocol-family-aware* qualifier sets and a
*curated, accession-registry-driven* index rather than an exhaustive
multi-species feature dump. Expand only as benchmark papers anchor the need.

## CURIE provenance policy

Every emitted CURIE is either:

- **Live-verified upstream** (UniProt / RefSeq / CHEBI / OBI / EFO / EDAM /
  MOD / PR) — cross-checked against the relevant REST API with a matching
  label before the entry shipped, OR
- **Locally namespaced** (`EpiHist:` / `EpiAssay:` / `EPIQC:`) — minted by
  this package under a documented `https://w3id.org/holobiomicslab/...`
  IRI because no stable upstream term exists, and explicitly marked as
  v0.1.0 local vocabulary scheduled for promotion.

Never assert an upstream CURIE without live verification. Local namespaces
are honest about provenance; fabricated upstream IDs are not.

## Versions

### v0.1.0 — initial release: three working adapters (current)

- **`BulkEpigenomicsAdapter`** — 19 bulk epigenomics qualifiers covering
  chromatin-binding, chromatin-accessibility, DNA-methylation, and
  chromosome-conformation assay families. `enrich_claim()` resolves
  repository accessions (GEO / SRA / ENCODE / EpiRR / EGA / BioProject /
  ArrayExpress / ENA / ChIP-Atlas / Cistrome), assay names (verified OBI /
  EFO + local `EpiAssay:`), histone marks (local `EpiHist:` + UniProt for
  variants), chromatin targets (107 UniProt IDs), genome assemblies (7
  RefSeq), and QC metrics (1 EDAM + 22 local `EPIQC:`).
- **`SingleCellEpigenomicsAdapter`** — 8 single-cell epigenomics qualifiers.
  `enrich_claim()` resolves the same accession / assay / histone / target /
  assembly indices as bulk; **CL / UBERON cell-type and tissue enrichment
  deferred to v0.2.0** (free-text cell-type labels pass through un-CURIE'd).
- **`MethylationArrayEpigenomicsAdapter`** — 8 array qualifiers. `enrich_claim()`
  resolves array-platform names (MethylationEPIC v1 → `OBI:0002131`; 450K /
  EPIC v2 / Mouse-EPIC → local `EpiAssay:`) and detects Illumina ``cg########``
  probe IDs (stored in `xrefs` only — no upstream CURIE registry for probes).
- **Bundled indices, populated and live-verified:**
  - `chromatin_targets.tsv` — 107 UniProt-verified human + mouse rows
  - `histone_marks.tsv` — 5 UniProt-verified variants, 3 CHEBI-verified DNA
    modifications (5mC / 5hmC / 5fC), 38 position-specific histone PTMs +
    2 DNA mods (5caC, 6mA) under local `EpiHist:` namespace pending upstream
    coverage
  - `assemblies.tsv` — 7 RefSeq-verified assemblies (human + mouse + T2T)
  - `encode_experiment_index.tsv` — 5,000 epigenomics-only experiments
    pulled from the live ENCODE portal REST API
  - `epirr_reference_index.tsv` — 325 IHEC reference epigenomes federated
    via the ENCODE `ReferenceEpigenome` endpoint (EBI EpiRR JSON API was
    retired in 2024–2025)
- **SHACL shapes** — full 385-triple `epigenomics.shapes.ttl` with one
  `sh:NodeShape` per adapter; pyshacl-loadable and CURIE-format constraints
  per qualifier.
- `tests/test_composition.py` enforces disjoint qualifier sets across all
  three sibling adapters.
- Pytest conformance + behaviour suite (165 tests).
- ≥ 90% test coverage enforced by CI (currently at 92.5%).

### v0.2.0 — CL / UBERON cell-type enrichment + monthly index refresh

- `SingleCellEpigenomicsAdapter.enrich_claim()` resolves cell-type labels
  (e.g. "CD8+ T cell") to CL CURIEs and tissue labels to UBERON CURIEs via
  a bundled `cell_type_index.tsv` (built from OLS4 exact-label match —
  same pattern as transcriptomics v0.2.0).
- Monthly GitHub Actions cron job running `data/build_encode_index.py` and
  `data/build_epirr_index.py` against the live ENCODE portal REST API
  (`Experiment` and `ReferenceEpigenome` endpoints respectively). EBI EpiRR
  remains HTML-only; ENCODE federation continues as the EpiRR substrate.
- Promote any `EpiHist:` / `EpiAssay:` entries to verified upstream CURIEs
  if upstream coverage improves (PR for histone proteoforms, EFO for
  newly-onboarded single-cell methylation protocols).

### v0.3.0+ — broader coverage

- Extend `chromatin_targets.tsv` with ChIP-Atlas and Cistrome curated
  antibody–target tables (optional enrichment, gated behind the `build`
  extra).
- Allow alternative ontology prefixes in the epigenomics index (e.g. GO for
  chromatin organisation, SO for sequence features).
- Submit upstream PRs to EFO / OBI for the verified gaps (HumanMethylation450
  is the most obvious — widely cited, no canonical CURIE). When accepted,
  promote the corresponding `EpiAssay:` entry.
- Contribute cnidarian-specific chromatin biology terms upstream to relevant
  ontologies where UBERON / GO gaps exist.
