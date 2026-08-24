# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added — classification vocabulary for ASB

- All three adapters now declare `primary_domain = "epigenomics"`, a
  `subdomains` tuple and a `techniques` tuple. ASB's classification layer
  unions these from every installed adapter
  (`agentic_science_builder.classification._adapter_vocab`); until now the
  adapters declared none, so an epigenomics paper classified as "other" and
  carried no technique at all.
- Both vocabularies are the 2026-08-25 cross-domain audit's set, split into
  disjoint slices across the three siblings: terms two sciences both declared
  (`mass-spectrometry`, `sequence-assembly`, `amplicon-sequencing`, ...) were
  dropped everywhere, since a shared term cannot tell the sciences apart.
- `techniques` is served by `techniques_for_class()` from a new
  `_AUDITED_TECHNIQUES` table in `_assays`, no longer derived from
  `_ASSAY_TABLE`: the resolver stays wider on purpose, so assay names whose
  tokens left the classification vocabulary (DNase-seq, MNase-seq, EM-seq,
  EPIC v2, the single-cell methylation protocols) still resolve to a CURIE.
- New public helpers in `_assays`: `techniques_for_class()` and
  `ASSAY_TECHNIQUES_BY_CLASS` (the resolver inventory, kept as its own view).

### Added — SEMAPV mapping-justifications (indicium-adapters A-convention)

- All three adapters' `enrich_claim` now populate an optional
  `ontology_term_justifications` companion-map alongside `ontology_terms`,
  keyed on the same slots, carrying the SEMAPV
  (https://w3id.org/semapv/vocab/) justification for *how* each CURIE was
  obtained: `semapv:CompositeMatching` for accession→CURIE registry
  resolution, `semapv:LexicalMatching` for the curated controlled-vocabulary
  name/label lookups. The mapping is shared in `_resolver.py`
  (`SOURCE_TO_SEMAPV` + `justification_for`), keyed on the resolver's `source`
  label, so the justification stays honest per mechanism. The xrefs-only
  probe-id path (no `ontology_terms`) correctly emits no justification.
  This is the companion-map convention documented by the `indicium-adapters`
  core `DomainAdapter` Protocol; downstream consumers (e.g. Perspicacité)
  promote `ontology_terms` + this map to first-class `sssom:Mapping` nodes.

### Changed — local development environment

- Added `[tool.uv.sources]` pinning `indicium` + `indicium-adapters` to editable
  sibling checkouts (`../indicium`, `../indicium-adapters`), mirroring the
  `indicium-adapters` / `-metabolomics` / `-transcriptomics` repos. Without it,
  `uv sync` resolves `indicium>=1.7` against the unrelated PyPI `indicium`
  package and fails. Ignored by consumers installing from the wheel.

---

## [0.2.0] - 2026-06-01

Domain-specific vocabulary expansion under the cross-adapter conventions
established by the transcriptomics sibling (qualifier-aware resolver,
bundled-TSV index pattern, local-vs-upstream CURIE provenance, snake_case
conformance). Each section below names the reusable pattern and the
epigenomics-specific content it carries; the goal is that a new adapter
author can lift the structural shape into their domain without
inheriting these specific concept choices.

Test suite: 165 → 353 tests. Coverage: 92.56% → 94.56%.

### Added — Shared qualifier-aware resolver (§1)

- `src/indicium_adapters_epigenomics/_resolver.py` — single
  `resolve(value, qualifier)` consulted by the three sibling adapters in
  place of three byte-identical `_resolve(value)` helpers. Walks
  ACCESSION first (regex-deterministic), then qualifier-preferred
  sources from `_QUALIFIER_PREFERENCES`, then the remaining default
  cascade. Accession always wins regardless of qualifier; preference
  is a hint, not a route — a preferred source that misses falls
  through to the default cascade.
- `tests/test_resolver.py` — 15 tests covering backwards-compatibility
  with v0.1.0 cascade (no qualifier), accession-outranks-preference
  invariant, qualifier preferences observed under synthetic ambiguity
  (mocked HISTONE/TARGET cross-hits), fall-through semantics, and
  structural guards on the preference map (snake_case keys, known
  source labels).

### Added — ChromHMM / chromatin-state vocabulary (§2)

- `EpiState:` namespace (`https://w3id.org/holobiomicslab/epistate/`)
  in `_namespaces.py`.
- `data/chromatin_states.tsv` — 42 curated rows covering ChromHMM
  15-state, 18-state and 25-state-imputed segmentation labels
  (TssA, EnhA1, ReprPC, Quies, ZNF/Rpts, …) plus the commonly
  interpreted Segway labels (Active_TSS, Polycomb, Heterochromatin,
  Quiescent). Each row carries a `models` field recording which
  segmentation models include it.
- `_chrom_states.py` with `ChromStateIndex` and broader normaliser
  collapsing `/`, `_`, `'` so `ZNF/Rpts` / `ZNF_Rpts` / `ZnfRpts` and
  `Tx5'` / `Tx 5` converge.
- New qualifiers `chromatin_state_of` (bulk) and
  `marker_state_of_cell_type` (single-cell). Disjoint by construction.
- `tests/test_chrom_states.py` — 29 tests.

### Added — DNA-modification chemistry depth (§3)

- `EpiMod:` namespace (`https://w3id.org/holobiomicslab/epimod/`).
- `data/dna_modifications.tsv` — 5 rows split out of `histone_marks.tsv`
  (5mC / 5hmC / 5fC / 5caC / 6mA). Each row carries an explicit
  `verified` column (`chebi` for upstream-verified, `local` for
  `EpiMod:` fallback). CHEBI-verified entries: 5mC `CHEBI:27551`, 5hmC
  `CHEBI:76792`, 5fC `CHEBI:76794`.
- `_dna_mods.py` with `DNAModIndex` and a paren-stripping normaliser so
  `5-(hydroxymethyl)cytosine` collides cleanly with
  `5-hydroxymethylcytosine`.
- Methylation-family qualifier preferences updated to
  `(DNA_MOD, TARGET, IHEC_BIOSAMPLE)` so chemistries, machinery
  enzymes (DNMTs, TETs, MBDs) and biosample names all resolve in
  the right order.
- `tests/test_dna_mods.py` — 40 tests including a regression guard
  that `HISTONE_INDEX` no longer returns DNA-mod rows.

### Added — IHEC reference biosample matrix (§4)

- `IHECBio:` namespace (`https://w3id.org/holobiomicslab/ihecbio/`).
- `data/ihec_biosample_index.tsv` — 33 curated rows covering 11 ENCODE
  Tier 1+2 cell lines (K562, GM12878, HepG2, HeLa-S3, …), 9 Roadmap
  primary / stem cells (CD4+ T cell, CD8+ T cell, CD14+ monocyte, …)
  and 13 reference tissues. CL / UBERON / EFO cross-refs populated
  only where the upstream term is unambiguous (K562 deliberately omits
  CL — K562's lineage in CL is contested).
- `_biosamples.py` with `IHECBiosampleIndex`. Cross-refs travel
  through `enrich_claim` via the `**match` splat into `xrefs[slot]`.
  Empty cross-refs are *absent* from the returned dict, not
  present-but-empty.
- Eight new accessibility / cell-type qualifier preferences
  (`accessible_in`, `accessible_in_cell_type`,
  `differentially_accessible_in`, `marker_peak_of_cell_type`,
  `marker_dmr_of_cell_type`, `enriched_in_cluster`,
  `localized_in_spatial_domain`, `co_accessible_with`).
- `tests/test_biosamples.py` — 27 tests.

### Added — ENCODE SCREEN cCRE classification (§5)

- `EpiCCRE:` namespace (`https://w3id.org/holobiomicslab/epiccre/`).
- `data/ccre_classes.tsv` — full 7-class SCREEN vocabulary (PLS, pELS,
  dELS, CA-CTCF, CA-H3K4me3, CA-TF, CA-only).
- `_ccre.py` with `CCREClassIndex`.
- `EH38E\d{7}` (hg38) and `EM10E\d{7}` (mm10) cCRE accessions added
  to `_registries._ACCESSION_PATTERNS`, routed to the canonical
  `ENCODE:` prefix with the `"ENCODE SCREEN"` repository label.
- New `classified_as_ccre` qualifier on bulk.
- `tests/test_ccre.py` — 34 tests including positive and negative
  accession-regex cases.

### Added — 3D-genome qualifier split (§6)

- `EpiConf:` namespace (`https://w3id.org/holobiomicslab/epiconf/`).
- `data/conformation_labels.tsv` — 7 compartment labels (A, B,
  A1, A2, B1, B2, B3 per Rao et al. 2014).
- `_conformation.py` with `ConformationIndex`.
- Six new bulk qualifiers split out of `interacts_3d_with`:
  `contacts_with` (generic Hi-C contact), `loops_with` (loop call),
  `anchored_to` (ChIA-PET / HiChIP anchor), `tad_boundary_with` (TAD
  boundary call), `compartment_of` (A/B compartment, resolved via
  CONFORMATION), `subtad_of` (sub-TAD nesting).
- `tests/test_conformation.py` — 38 tests including a `pytest.warns`
  visibility check for the deprecation alias.

### Changed

- `chromatin_targets.tsv` — corrected `category` for the DNA-methylation
  machinery enzymes: TET1/2/3 retagged from `chromatin_remodeler`
  → `dna_demethylase`; MBD1-4 retagged → `methyl_cpg_reader`. v0.1.0
  had these mis-classified.
- `_histone_marks.py` docstring + `histone_marks.tsv` updated to
  reflect that DNA-modification rows moved to `_dna_mods.py`.
- `BulkEpigenomicsAdapter.extraction_context()` updated to mention
  the new qualifiers and the EH38E / EM10E accession formats.
- `SingleCellEpigenomicsAdapter.extraction_context()` updated to
  mention `marker_state_of_cell_type` and ChromHMM state labels.
- `tests/test_singlecell.py::test_cell_type_label_in_object_slot_is_not_curie_d_in_v0_1_0`
  refactored into two parametrised tests — labels NOT in the IHEC
  seed (neuron, hepatocyte, B cell) still pass through un-CURIE'd;
  labels IN the seed (CD8+/CD4+ T cell, K562, GM12878) now resolve
  to `IHECBio:`. Pins the v0.2.0 contract.

### Deprecated

- `BulkEpigenomicsAdapter.qualifiers` still contains
  `interacts_3d_with`, but `enrich_claim()` now emits a
  `DeprecationWarning` when it sees the qualifier. The warning
  message lists every replacement qualifier. The deprecated qualifier
  still enriches under `contacts_with` semantics. Removal target:
  **v0.3.0**.

### Compatibility

Every claim that resolved under v0.1.0 still resolves to the same
source label and CURIE under v0.2.0 — `tests/test_resolver.py` pins
the v0.1.0 fixed-cascade outputs as regression guards. The only
behavioural change to existing qualifiers is `marker_peak_of_cell_type`,
whose preference moved from `(HISTONE, TARGET)` to
`(IHEC_BIOSAMPLE, HISTONE, TARGET)` — a strict superset that gains a
biosample resolution attempt without losing the original ones.

---

## [0.1.1] - 2026-06-01

### Changed

- **Document the snake_case constraint** on the `qualifiers` frozenset of
  `BulkEpigenomicsAdapter`, `SingleCellEpigenomicsAdapter`, and
  `MethylationArrayEpigenomicsAdapter`. Every qualifier MUST match
  `^[a-z][a-z0-9_]*$`; non-conforming values would be silently rejected by
  `indicium.validate_graph` because indicium >= 1.7.0 routes adapter
  qualifiers through the open `asb:domainQualifier` slot, which enforces
  exactly this pattern at the SHACL boundary.
- Bumped dependency pins: `indicium>=1.7`, `indicium-adapters>=0.5.0`.

### Added

- `tests/test_protocol_conformance.py::test_every_discovered_adapter_qualifier_is_snake_case`
  — defence-in-depth conformance test that walks every adapter discovered
  via the `indicium.adapters` entry-point group (not just the hardcoded
  `_ADAPTERS` tuple). Catches the case where a future sibling adapter is
  registered in `pyproject.toml` but accidentally omitted from the
  hardcoded test table.

### Why

Mirrors the snake_case hardening already merged in
[`indicium-adapters-transcriptomics#1`](https://github.com/HolobiomicsLab/indicium-adapters-transcriptomics/issues/1)
(commits `40ff1211` + `1deb1f10`) and the open follow-up on
[`indicium-adapters-metabolomics#9`](https://github.com/HolobiomicsLab/indicium-adapters-metabolomics/issues/9),
both of which are follow-ups to
[`HolobiomicsLab/indicium#1`](https://github.com/HolobiomicsLab/indicium/issues/1)
resolved in indicium 1.7.0. All current epigenomics qualifiers already
comply with the pattern (verified by the regression test); this release
documents the contract and guards against future drift.

### Compatibility

Additive only. No qualifier names, namespaces, or public API surfaces
changed. v0.1.0 consumers continue to work unchanged.

---

## [0.1.0] - 2026-05-27

### Added

Initial release. Three sibling `DomainAdapter` classes registered under the
`indicium.adapters` entry-point group:

- **`BulkEpigenomicsAdapter`** — bulk epigenomics assay claims (WGBS / RRBS /
  EM-seq, ChIP-seq / CUT&RUN / CUT&Tag, ATAC-seq / MNase-seq / DNase-seq,
  Hi-C / Micro-C / ChIA-PET / HiChIP). 19 qualifiers covering chromatin
  binding, marking, accessibility, methylation, 3D interaction, and
  study-level assay metadata.
- **`SingleCellEpigenomicsAdapter`** — single-cell and single-nucleus
  epigenomics claims (scATAC / snATAC / snmC-seq / sci-MET / snmCAT-seq /
  scChIP / scNOMe / scHi-C). 8 cell-resolution qualifiers, disjoint from
  the bulk and array sets.
- **`MethylationArrayEpigenomicsAdapter`** — Illumina methylation array
  claims (450K / EPIC v1 / EPIC v2 / MethylationEPIC mouse). 8
  array-workflow qualifiers; Illumina `cg########` probe-ID recognition
  surfaced via `xrefs`.

Bundled live-verified data:

- 107 chromatin targets (UniProt-verified, human + mouse)
- 5 UniProt-verified histone variants, 3 CHEBI-verified DNA modifications
  (5mC / 5hmC / 5fC), 38 position-specific histone PTMs under local
  `EpiHist:` namespace pending upstream PSI-MOD / PR coverage
- 7 RefSeq-verified reference-genome assemblies
- 5,000 ENCODE experiments (`encode_experiment_index.tsv`) pulled from the
  ENCODE portal REST API
- 325 IHEC reference epigenomes (`epirr_reference_index.tsv`) federated via
  the ENCODE `ReferenceEpigenome` endpoint (EBI EpiRR JSON API retired
  2024–2025)
- 9 OBI-verified + 8 EFO-verified assay terms + 12 local `EpiAssay:`
  entries for assays without a verified upstream term
- 1 EDAM-verified + 22 local `EPIQC:` QC metric names

SHACL shapes (`shapes/epigenomics.shapes.ttl`, 385 triples, one
`sh:NodeShape` per adapter) and disjoint-qualifier-set composition test
(`tests/test_composition.py`). 165-test pytest suite with ≥90% coverage
gate.
