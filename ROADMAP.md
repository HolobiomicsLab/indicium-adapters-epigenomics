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

### v0.1.1 — indicium v1.7.0 conformance (shipped)

- Documented the snake_case constraint (`^[a-z][a-z0-9_]*$`) on every
  adapter's `qualifiers` frozenset. Required by indicium >= 1.7.0:
  Perspicacité routes adapter qualifiers to the open `asb:domainQualifier`
  slot whose SHACL shape enforces exactly this pattern; non-conforming
  values would be silently rejected by `indicium.validate_graph`.
- Added `test_protocol_conformance::test_every_discovered_adapter_qualifier_is_snake_case`
  — defence-in-depth check that walks every entry-point under
  `indicium.adapters` (not just the hardcoded test table) so a future
  sibling registered in `pyproject.toml` but missing from the tuple is
  still caught.
- Bumped dependency pins: `indicium>=1.7`, `indicium-adapters>=0.5.0`.

### v0.2.0 — qualifier-aware resolver + domain-specific vocabulary expansion (shipped)

This wave applies cross-adapter conventions (qualifier-aware resolver,
bundled-TSV index pattern, local-vs-upstream CURIE provenance) to
epigenomics-specific concepts. Each §addition follows the same shape so
future adapter authors can reuse the structural pattern in a new
domain.

- **§1 — Shared qualifier-aware resolver** (`_resolver.py`). Three
  byte-identical `_resolve(value)` helpers collapsed into one
  `resolve(value, qualifier)` that consults a qualifier → preferred-source
  preference map. Accessions always win first (regex-deterministic). The
  default cascade preserves v0.1.0 behaviour byte-for-byte when no
  qualifier preference exists. Preference is a hint, not a route — a
  preferred source that misses still falls through to the rest of the
  cascade.
- **§2 — ChromHMM / chromatin-state vocabulary** (`EpiState:`, 42 rows).
  ENCODE3 / Roadmap ChromHMM 15-state, 18-state and 25-state-imputed
  segmentation labels plus the commonly interpreted Segway labels, with
  per-row `models` membership metadata. New qualifiers
  `chromatin_state_of` (bulk) and `marker_state_of_cell_type`
  (single-cell). No upstream CURIEs exist for segmentation labels at
  v0.2.0 release.
- **§3 — DNA-modification chemistry depth** (`EpiMod:`, `DNA_MOD_INDEX`).
  5mC / 5hmC / 5fC / 5caC / 6mA split out of `HISTONE_INDEX` into a
  dedicated `_dna_mods.py`. CHEBI-verified entries (5mC `CHEBI:27551`,
  5hmC `CHEBI:76792`, 5fC `CHEBI:76794`) emit under the canonical
  upstream prefix; the rest under `EpiMod:`. Each row carries a
  `verified` column (`chebi` | `local`) recording the provenance
  authority. TET1/2/3 retagged from `chromatin_remodeler` to
  `dna_demethylase`; MBD1-4 to `methyl_cpg_reader` — v0.1.0 had these
  mis-categorised. Methylation-family qualifiers
  (`methylated_in`, `methylated_in_cell_type`,
  `differentially_methylated_in`, `methylated_at_locus`,
  `differentially_methylated_at_probe`) route through `(DNA_MOD,
  TARGET, IHEC_BIOSAMPLE)` so chemistries, enzymes and biosample names
  all resolve cleanly.
- **§4 — IHEC reference biosample matrix** (`IHECBio:`, 33 rows). The
  IHEC core biosample taxonomy is the join key across ENCODE / Roadmap /
  DEEP / BLUEPRINT reference epigenomes. Bundled seed covers 11 ENCODE
  Tier 1+2 cell lines, 9 Roadmap primary / stem cells, 13 reference
  tissues, with CL / UBERON / EFO cross-refs populated only where the
  upstream term is unambiguous (K562 deliberately omits a CL cross-ref —
  K562's CL lineage is contested). Cross-refs travel through
  `enrich_claim` via the `**match` splat into `xrefs[slot]`. New
  qualifier preferences for `accessible_in`, `accessible_in_cell_type`,
  `differentially_accessible_in`, `marker_peak_of_cell_type`,
  `marker_dmr_of_cell_type`, `enriched_in_cluster`,
  `localized_in_spatial_domain`, `co_accessible_with`.
- **§5 — ENCODE SCREEN cCRE classification** (`EpiCCRE:`, 7 classes).
  Full SCREEN class vocabulary (PLS / pELS / dELS / CA-CTCF / CA-H3K4me3
  / CA-TF / CA-only). `EH38E\d{7}` (hg38) and `EM10E\d{7}` (mm10) cCRE
  accessions recognised in `_registries.py` and routed to the canonical
  `ENCODE:` prefix with the `"ENCODE SCREEN"` repository label. New
  `classified_as_ccre` qualifier on bulk.
- **§6 — 3D-genome qualifier split** (`EpiConf:`, 7 compartment labels).
  `interacts_3d_with` collapsed every 3D observation into one slot;
  v0.2.0 splits it by semantics: `contacts_with` (generic Hi-C contact),
  `loops_with` (loop call), `anchored_to` (ChIA-PET / HiChIP anchor),
  `tad_boundary_with` (TAD boundary), `compartment_of` (A / B / A1 / A2
  / B1 / B2 / B3 via `CONFORMATION_INDEX`), `subtad_of` (sub-TAD
  nesting). `interacts_3d_with` retained for one minor version with a
  `DeprecationWarning` emitted from `enrich_claim` — removal target
  v0.3.0.

Qualifier-set union grew from 35 (v0.1.0) to 44 (v0.2.0). Test suite
grew from 165 → 353 with coverage 92.56% → 94.56%. No qualifier name
or namespace was renamed; the only breaking change (`interacts_3d_with`)
is announced via `DeprecationWarning` with the full migration menu in
the message.

### v0.3.0+ — broader coverage and ecosystem integration

- **§7 — Illumina methylation-array probe → genomic context.** Build a
  lazy-loaded probe index from the official 450K / EPIC v1 / EPIC v2 /
  MouseMethylation manifests (450K–935K rows each, gzipped). Extend
  `MethylationArrayEpigenomicsAdapter` so a `methylated_at_locus` claim
  carrying a `cg########` subject attaches an `asb:genomicContext` field
  (chromosome, position, strand, nearest gene symbol, CpG-island
  relation, probe type). Maintainer script
  `data/build_illumina_manifest_index.py` gated behind the `build`
  extra.
- **§8 — Chromatin-organisation GO + SO terms.** OLS4-cascade build
  script scoped to chromatin-relevant subtrees: `GO:0006325 chromatin
  organization`, `GO:0006338 chromatin remodeling`, `GO:0016570 histone
  modification`, `GO:0006342 chromatin silencing`; `SO:0001747
  open_chromatin_region`, `SO:0001911 ATAC-seq_open_chromatin_region`,
  `SO:0001748 CpG_island`, `SO:0001841 binding_site`, `SO:0000167
  promoter`, `SO:0000165 enhancer`. Wire to `enriched_for`,
  `differentially_methylated_in`, `differentially_accessible_in`.
- **§9 — ChIP-Atlas + Cistrome curated antibody–target tables.** Merge
  into `TARGET_INDEX` as `xrefs` attachments without breaking the
  existing 107-row UniProt-verified baseline. Builder scripts gated
  behind the `build` extra.
- **§10 — Opt-in `verify_quote` integration.** With indicium v1.7.0+
  exporting `verify_quote` and `VerificationResult`, accept a
  `verify_quotes=False` constructor flag on each adapter; when on, claims
  carrying `asb:provenance.quote` + `asb:provenance.passage_text` get an
  `asb:anchorStatus ∈ {verified, repaired, unverified, unchecked}` tag,
  with `asb:quoteExact` populated only for verified / repaired matches
  (no laundering). Pattern mirrors Perspicacité PR #10.
- **§11 — Upstream PSI-MOD audit + promotion.** Run `data/audit_psi_mod_coverage.py`
  to query OLS4 for each of the 38 `EpiHist:` position-specific PTMs
  by exact label and report `already-in-PSI-MOD` / `still-missing`.
  Promote in-place for the former; file upstream PSI-MOD submissions for
  the latter via the public ontology-submissions list. Same audit
  pattern applies to `EpiAssay:` (EFO / OBI) and `EpiMod:` (CHEBI).
- **`interacts_3d_with` removal.** The deprecated qualifier emits
  `DeprecationWarning` in v0.2.0 — v0.3.0 removes it from the qualifiers
  frozenset. Migration: any of `contacts_with`, `loops_with`,
  `anchored_to`, `tad_boundary_with`, `compartment_of`, `subtad_of`
  depending on the assay-specific semantics.
- **Monthly index-refresh GitHub Action.** Schedule the existing
  `data/build_encode_index.py` + `data/build_epirr_index.py` (and the
  new §4 / §7 / §8 / §9 build scripts) under
  `.github/workflows/refresh-indices.yml`. Mirrors the transcriptomics
  cron workflow.
- **End-to-end demo.** Pick an ENCODE epigenomics reference experiment
  (e.g. K562 ChIP-seq for CTCF — ENCSR000AKB) and walk extraction →
  enrichment → SHACL validation → indicium emit, end-to-end. Land as
  `examples/encsr000akb_demo.py` + `examples/README.md`. Pattern mirrors
  transcriptomics' Warner2018 demo; content is canonical chromatin
  binding rather than bulk RNA-seq.
- **Contribute cnidarian-specific chromatin biology terms upstream**
  where UBERON / GO gaps exist (carried over from the v0.1.0 roadmap).
