"""Free-text → CURIE helpers for epigenomics repository accessions.

Pure functions, no I/O. Used by the three epigenomics adapter
``enrich_claim()`` implementations and by upstream bridge layers that want to
normalise study / sample / experiment accessions before emitting claims.

Patterns are ordered from most-specific to least-specific so that e.g. an
ENCODE identifier is never partially matched by the looser SRA pattern, and
the IHEC reference-epigenome pattern wins over the generic ``IHEC*`` family.
"""
from __future__ import annotations

import re
from typing import Final

# Accession patterns. Kept conservative: anchored, well-known prefixes only.
# Each entry maps a compiled regex → (prefix, repository_label). The CURIE
# prefix and the human-readable repository label can diverge (e.g. the IHEC
# reference-epigenome CURIE uses the ``EpiRR`` prefix per identifiers.org but
# the repository is colloquially called "IHEC").
_ACCESSION_PATTERNS: Final[tuple[tuple[re.Pattern[str], str, str], ...]] = (
    # ENCODE — fixed-width type-tag + 6 alphanumerics (digits or letters in any
    # order; ENCODE assigns both forms in the wild). Must precede SRA / GEO.
    (re.compile(r"^ENC[A-Z]{2}[A-Z0-9]{6}$"), "ENCODE", "ENCODE"),
    # ArrayExpress / BioStudies-AE
    (re.compile(r"^E-[A-Z]{4}-\d+$"), "ArrayExpress", "ArrayExpress"),
    # EGA — 11 trailing digits, type letter S/D/F/B/N
    (re.compile(r"^EGA[SDFBN]\d{11}$"), "EGA", "EGA"),
    # BioProject — NA (NCBI), EB (ENA), DB (DDBJ)
    (re.compile(r"^PRJ(?:NA|EB|DB)\d+$"), "BioProject", "BioProject"),
    # IHEC reference epigenomes — the EpiRR-style accession
    # (``IHECRE00000001.1``) must precede the generic ``IHEC*`` family so it
    # always claims the EpiRR CURIE rather than being mis-routed.
    (re.compile(r"^IHECRE\d+\.\d+$"), "EpiRR", "IHEC"),
    (re.compile(r"^IHEC[-A-Z0-9]+\d+$"), "EpiRR", "IHEC"),
    # GEO — GSE = series, GSM = sample, GPL = platform, GDS = dataset
    (re.compile(r"^G(?:SE|SM|PL|DS)\d{3,7}$"), "GEO", "GEO"),
    # SRA / ENA / DDBJ — [SED] prefix + [RXPS] kind + 4+ digits
    (re.compile(r"^SR[RXPS]\d{4,}$"), "SRA", "SRA"),
    (re.compile(r"^ER[RXPS]\d{4,}$"), "ENA", "ENA"),
    (re.compile(r"^DR[RXPS]\d{4,}$"), "SRA", "DDBJ"),
    # ChIP-Atlas portal accessions — ``SRX`` is already matched above; ChIP-
    # Atlas also exposes ``ChIP-Atlas:<id>`` and ``CA<digits>`` external IDs.
    (re.compile(r"^(?:ChIP-Atlas|CA):?(\d+)$"), "ChIPAtlas", "ChIP-Atlas"),
    # Cistrome — numeric IDs alone (e.g. ``45123``) are too loose to
    # autodetect, so we only match the prefixed forms (``Cistrome:45123`` /
    # ``Cistrome_45123`` / ``Cistrome-45123``).
    (re.compile(r"^Cistrome[:_-](\d{3,7})$"), "Cistrome", "Cistrome"),
)


def accession_to_curie(text: str) -> dict[str, str] | None:
    """Return repository CURIE info for a recognised accession string.

    Returns a dict with keys ``curie`` (``"<prefix>:<bare_id>"``),
    ``repository`` (the colloquial repository name) and ``bare_id`` (the
    accession string stripped of any prefix decoration), or ``None`` when
    nothing matches.

    >>> accession_to_curie("GSE168938")
    {'curie': 'GEO:GSE168938', 'repository': 'GEO', 'bare_id': 'GSE168938'}
    >>> accession_to_curie("ENCSR000AKO")
    {'curie': 'ENCODE:ENCSR000AKO', 'repository': 'ENCODE', 'bare_id': 'ENCSR000AKO'}
    >>> accession_to_curie("IHECRE00000001.1")
    {'curie': 'EpiRR:IHECRE00000001.1', 'repository': 'IHEC', 'bare_id': 'IHECRE00000001.1'}
    >>> accession_to_curie("not-a-real-accession") is None
    True
    """
    if not isinstance(text, str):
        return None
    token = text.strip()
    if not token:
        return None
    for pattern, prefix, repository in _ACCESSION_PATTERNS:
        if m := pattern.match(token):
            # Prefixed-form patterns (ChIP-Atlas, Cistrome) capture the bare
            # numeric id in group 1; for the rest the whole match IS the
            # bare id.
            bare_id = m.group(1) if m.groups() else token
            return {
                "curie": f"{prefix}:{bare_id}",
                "repository": repository,
                "bare_id": bare_id,
            }
    return None
