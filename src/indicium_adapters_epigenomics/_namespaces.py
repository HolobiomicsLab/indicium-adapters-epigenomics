"""Hardcoded namespace URI constants for epigenomics ontologies.

We deliberately do not depend on any ontology Python package at runtime; the
URIs are stable and curated here. Update with care.

Cross-adapter consistency: prefixes that overlap with
``indicium-adapters-transcriptomics`` (EFO, OBI, NCBITaxon, CL, UBERON, GO,
SO, MONDO, PATO, STATO, EDAM, RefSeq, UniProt, GEO, SRA, BioProject, ENA,
ArrayExpress, ENCODE) use exactly the same IRIs — claims emitted by sibling
adapters share a single namespace mapping in the consuming graph store.
"""
from __future__ import annotations

from typing import Final

# OBO Foundry — biomedical / experimental / phenotypic
EFO: Final[str] = "http://www.ebi.ac.uk/efo/EFO_"
OBI: Final[str] = "http://purl.obolibrary.org/obo/OBI_"
STATO: Final[str] = "http://purl.obolibrary.org/obo/STATO_"
NCBITAXON: Final[str] = "http://purl.obolibrary.org/obo/NCBITaxon_"
CL: Final[str] = "http://purl.obolibrary.org/obo/CL_"
UBERON: Final[str] = "http://purl.obolibrary.org/obo/UBERON_"
GO: Final[str] = "http://purl.obolibrary.org/obo/GO_"
SO: Final[str] = "http://purl.obolibrary.org/obo/SO_"
MONDO: Final[str] = "http://purl.obolibrary.org/obo/MONDO_"
PATO: Final[str] = "http://purl.obolibrary.org/obo/PATO_"
EDAM: Final[str] = "http://edamontology.org/"

# Chromatin / sequence-feature / protein modifications
PR: Final[str] = "http://purl.obolibrary.org/obo/PR_"
MOD: Final[str] = "http://purl.obolibrary.org/obo/MOD_"      # PSI-MOD
CHEBI: Final[str] = "http://purl.obolibrary.org/obo/CHEBI_"

# Protein authorities
UNIPROT: Final[str] = "https://identifiers.org/uniprot:"

# Genome assemblies / annotation
REFSEQ: Final[str] = "https://identifiers.org/refseq:"
GENCODE: Final[str] = "https://identifiers.org/gencode:"
ENSEMBL: Final[str] = "https://identifiers.org/ensembl:"
UCSC: Final[str] = "https://identifiers.org/ucsc:"

# Repositories / accessions
GEO: Final[str] = "https://identifiers.org/geo:"
SRA: Final[str] = "https://identifiers.org/insdc.sra:"
ARRAYEXPRESS: Final[str] = "https://identifiers.org/arrayexpress:"
ENCODE: Final[str] = "https://identifiers.org/encode:"
IHEC: Final[str] = "https://identifiers.org/epirr:"          # IHEC / EpiRR
EGA: Final[str] = "https://identifiers.org/ega.study:"
BIOPROJECT: Final[str] = "https://identifiers.org/bioproject:"
ENA: Final[str] = "https://identifiers.org/ena.embl:"
CHIP_ATLAS: Final[str] = "https://identifiers.org/chipatlas:"
CISTROME: Final[str] = "https://identifiers.org/cistrome:"

# Local vocabulary; minted in v0.1.0; promote to STATO/EDAM when terms land
# upstream. Used by ``_qc_metrics.py`` for QC metric names that have no
# stable STATO / EDAM operation CURIE at the time of writing (e.g.
# ``EPIQC:bisulfite-conversion-rate``, ``EPIQC:houseman-reference``).
EPIQC: Final[str] = "https://w3id.org/holobiomicslab/epiqc/"

# Local vocabulary for histone-modification SYMBOLS as used in publications
# (e.g. ``EpiHist:H3K27ac``). PSI-MOD (MOD:) coverage of position-specific
# histone PTMs is incomplete and inconsistent across releases; the v0.1.0
# adapter therefore ships local CURIEs and a roadmap entry to promote them
# to PRO (PR:) proteoform identifiers when matched terms are available. The
# v0.1.0 file ALSO uses CHEBI for the small set of DNA modifications and
# UniProt for histone variants, both of which have stable upstream IDs.
EPIHIST: Final[str] = "https://w3id.org/holobiomicslab/epihist/"

# Local vocabulary for assay-type symbols where no stable OBI/EFO CURIE
# exists in v0.1.0 (e.g. ``EpiAssay:CUT&RUN``, ``EpiAssay:HumanMethylation450``).
# The verified-via-OLS4 CURIEs are kept under their canonical OBI/EFO
# prefix in ``_assays.py``; everything else lives here pending upstream
# coverage. Cross-checked against OLS4 on 2026-05-27 — many community-
# popular assays (CUT&RUN, CUT&Tag, MethylationEPIC, single-cell methylation
# protocols, HiChIP, Micro-C, EM-seq, oxBS-seq) had no verified term.
EPIASSAY: Final[str] = "https://w3id.org/holobiomicslab/epiassay/"

# Local vocabulary for chromatin-state labels emitted by genome-segmentation
# methods (ChromHMM 15-state / 18-state / 25-state ENCODE3 imputed models,
# Segway commonly-interpreted labels). ENCODE and Roadmap publish these as
# documentation tables, not queryable APIs, and no upstream ontology has yet
# minted CURIEs for the individual segmentation labels. ``_chrom_states.py``
# ships the curated table; promote a label to a SO / GO sequence-feature
# CURIE if the upstream coverage ever lands.
EPISTATE: Final[str] = "https://w3id.org/holobiomicslab/epistate/"

# Local vocabulary for DNA-modification chemistries without a verified
# upstream CHEBI ID at the time of writing (v0.2.0). Split out of the
# ``EpiHist:`` namespace in v0.2.0 §3 to make the histone-PTM vs
# DNA-chemistry distinction explicit at the CURIE level. Verified
# upstream CHEBI IDs (5mC, 5hmC, 5fC at v0.2.0; 5caC / 6mA / 4mC if
# upstream coverage lands) are emitted under their canonical
# ``CHEBI:`` prefix in ``_dna_mods.py``; only entries without a verified
# upstream live here.
EPIMOD: Final[str] = "https://w3id.org/holobiomicslab/epimod/"

# Local vocabulary for IHEC core biosample terms (v0.2.0 §4). IHEC
# publishes a controlled biosample taxonomy used as the join key across
# the ENCODE, Roadmap Epigenomics, DEEP, and BLUEPRINT reference
# epigenome data sets; the matrix lives in the "IHEC Standard
# Specifications and Required Metadata" document and is not surfaced as
# an OBO-style ontology. Each ``IHECBio:`` row carries cross-references
# to Cell Ontology (CL), UBERON, ENCODE biosample term, and EFO so the
# enriched claim payload is interoperable with consumers that prefer
# those upstream vocabularies. The IRI is reserved here; if IHEC ever
# publishes its taxonomy as a public IRI scheme we promote in place.
IHECBIO: Final[str] = "https://w3id.org/holobiomicslab/ihecbio/"

# Local vocabulary for ENCODE SCREEN candidate cis-Regulatory Element
# (cCRE) class labels (v0.2.0 §5). ENCODE SCREEN publishes a 7-class
# functional classification (PLS, pELS, dELS, CA-CTCF, CA-H3K4me3,
# CA-TF, CA-only) for genome-wide cCRE annotation; the classes are
# documented but no upstream registry mints CURIEs for them. ``EpiCCRE:``
# is the local namespace for the class labels themselves; individual
# cCRE accessions (``EH38E0000001`` for hg38, ``EM10E0000001`` for mm10)
# remain ENCODE accessions and emit under the canonical ``ENCODE:``
# prefix via ``_registries.py``.
EPICCRE: Final[str] = "https://w3id.org/holobiomicslab/epiccre/"

# Local vocabulary for 3D-genome compartment / sub-compartment labels
# (v0.2.0 §6). Hi-C / Micro-C analyses commonly call A / B compartments
# and the Rao et al. (2014) sub-compartment refinement A1 / A2 / B1 /
# B2 / B3. No upstream ontology mints CURIEs for these labels; ``EpiConf:``
# is the local namespace. Used by the ``compartment_of`` qualifier
# introduced in the 3D-genome qualifier split (§6).
EPICONF: Final[str] = "https://w3id.org/holobiomicslab/epiconf/"

EPIGENOMICS_PREFIXES: Final[dict[str, str]] = {
    # Biomedical / experimental / phenotypic
    "EFO": EFO,
    "OBI": OBI,
    "STATO": STATO,
    "NCBITaxon": NCBITAXON,
    "CL": CL,
    "UBERON": UBERON,
    "GO": GO,
    "SO": SO,
    "MONDO": MONDO,
    "PATO": PATO,
    "EDAM": EDAM,
    # Chromatin / sequence-feature / protein modifications
    "PR": PR,
    "MOD": MOD,
    "CHEBI": CHEBI,
    # Protein authorities
    "UniProt": UNIPROT,
    # Genome assemblies / annotation
    "RefSeq": REFSEQ,
    "GENCODE": GENCODE,
    "ENSEMBL": ENSEMBL,
    "UCSC": UCSC,
    # Repositories
    "GEO": GEO,
    "SRA": SRA,
    "ArrayExpress": ARRAYEXPRESS,
    "ENCODE": ENCODE,
    "EpiRR": IHEC,
    "EGA": EGA,
    "BioProject": BIOPROJECT,
    "ENA": ENA,
    "ChIPAtlas": CHIP_ATLAS,
    "Cistrome": CISTROME,
    # Local vocabulary
    "EPIQC": EPIQC,
    "EpiHist": EPIHIST,
    "EpiAssay": EPIASSAY,
    "EpiState": EPISTATE,
    "EpiMod": EPIMOD,
    "IHECBio": IHECBIO,
    "EpiCCRE": EPICCRE,
    "EpiConf": EPICONF,
}
