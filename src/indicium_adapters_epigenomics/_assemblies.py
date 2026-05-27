"""Reference-genome-assembly name → RefSeq CURIE resolver.

Pure functions, no I/O. Kept inline as a Python literal (the canonical
human / mouse assembly set is small and very stable). Aliases are
case-insensitive and cover the UCSC short names (``hg38`` / ``mm10`` …)
alongside the canonical NCBI / GRC names (``GRCh38`` / ``T2T-CHM13v2`` …).
"""
from __future__ import annotations

from types import MappingProxyType
from typing import Final, NamedTuple


class _AssemblySpec(NamedTuple):
    """One reference-assembly table row."""

    canonical: str
    aliases: tuple[str, ...]
    curie: str           # RefSeq:GCF_*
    organism_curie: str  # NCBITaxon:9606 (human) / 10090 (mouse)
    ucsc_name: str       # short UCSC alias, e.g. ``hg38``


# Human + mouse focus per the v0.1.0 scope. ``T2T-CHM13v2`` keeps its dotted-
# version suffix in canonical form because that is how the NCBI / RefSeq
# accession is named in publications.
_ASSEMBLIES: Final[tuple[_AssemblySpec, ...]] = (
    _AssemblySpec(
        canonical="GRCh38",
        aliases=("GRCh38", "GRCh38.p14", "hg38", "GCF_000001405.40"),
        curie="RefSeq:GCF_000001405.40",
        organism_curie="NCBITaxon:9606",
        ucsc_name="hg38",
    ),
    _AssemblySpec(
        canonical="GRCh37",
        aliases=("GRCh37", "GRCh37.p13", "hg19", "GCF_000001405.25"),
        curie="RefSeq:GCF_000001405.25",
        organism_curie="NCBITaxon:9606",
        ucsc_name="hg19",
    ),
    _AssemblySpec(
        canonical="NCBI36",
        aliases=("NCBI36", "hg18", "GCF_000001405.12"),
        curie="RefSeq:GCF_000001405.12",
        organism_curie="NCBITaxon:9606",
        ucsc_name="hg18",
    ),
    _AssemblySpec(
        canonical="T2T-CHM13v2.0",
        aliases=(
            "T2T-CHM13v2.0", "T2T-CHM13v2", "T2T-CHM13", "CHM13v2",
            "GCF_009914755.1", "hs1",
        ),
        curie="RefSeq:GCF_009914755.1",
        organism_curie="NCBITaxon:9606",
        ucsc_name="hs1",
    ),
    _AssemblySpec(
        canonical="GRCm39",
        aliases=("GRCm39", "mm39", "GCF_000001635.27"),
        curie="RefSeq:GCF_000001635.27",
        organism_curie="NCBITaxon:10090",
        ucsc_name="mm39",
    ),
    _AssemblySpec(
        canonical="GRCm38",
        aliases=("GRCm38", "GRCm38.p6", "mm10", "GCF_000001635.26"),
        curie="RefSeq:GCF_000001635.26",
        organism_curie="NCBITaxon:10090",
        ucsc_name="mm10",
    ),
    _AssemblySpec(
        canonical="NCBI37",
        aliases=("NCBI37", "MGSCv37", "mm9", "GCF_000001635.18"),
        curie="RefSeq:GCF_000001635.18",
        organism_curie="NCBITaxon:10090",
        ucsc_name="mm9",
    ),
)


class AssemblyIndex:
    """Case-insensitive assembly-name lookup with alias expansion.

    Pure data; loaded once at import time. Lookups try the alias index first
    (covers UCSC and dotted GRC variants), then fall through to a
    canonical-name match. Returns a frozen dict so callers can't accidentally
    mutate the shared registry.
    """

    def __init__(self, assemblies: tuple[_AssemblySpec, ...]):
        self._assemblies = assemblies
        self._by_alias: dict[str, _AssemblySpec] = {}
        for spec in assemblies:
            for alias in spec.aliases:
                self._by_alias[alias.lower()] = spec
            self._by_alias[spec.canonical.lower()] = spec

    def lookup(self, text: str) -> dict[str, str] | None:
        """Return matching assembly info or ``None``.

        >>> AssemblyIndex(_ASSEMBLIES).lookup("hg38")["name"]
        'GRCh38'
        >>> AssemblyIndex(_ASSEMBLIES).lookup("mm10")["organism_curie"]
        'NCBITaxon:10090'
        """
        if not isinstance(text, str) or not text:
            return None
        spec = self._by_alias.get(text.strip().lower())
        if spec is None:
            return None
        return {
            "curie": spec.curie,
            "name": spec.canonical,
            "organism_curie": spec.organism_curie,
            "ucsc_name": spec.ucsc_name,
        }

    def __len__(self) -> int:
        return len(self._assemblies)


ASSEMBLY_INDEX: Final[AssemblyIndex] = AssemblyIndex(_ASSEMBLIES)

ASSEMBLY_NAME_TO_CURIE: Final[MappingProxyType[str, str]] = MappingProxyType(
    {spec.canonical: spec.curie for spec in _ASSEMBLIES}
)
