"""IHEC core biosample → ``IHECBio:`` CURIE resolver with CL / UBERON / EFO
cross-references.

The IHEC core biosample matrix is the lingua franca for joining ENCODE,
Roadmap Epigenomics, DEEP, and BLUEPRINT reference epigenomes. It lives in
the IHEC "Standard Specifications and Required Metadata" document as a
table rather than an OBO-style ontology; no upstream IRI scheme has been
minted for the matrix at v0.2.0 release. ``_biosamples.py`` ships a
curated seed covering the most-cited ENCODE Tier 1 + Tier 2 cell lines,
Roadmap primary hematopoietic and immune cells, and the ENCODE / Roadmap
canonical reference tissues.

Why an IHEC index, not just a Cell Ontology index? Cell-line free-text
labels in the literature (``K562``, ``HeLa-S3``, ``H1-hESC``) routinely
have ambiguous or contested CL mappings — K562's lineage-of-origin
within CL is still debated — while their IHEC biosample-term names are
stable, consistent across the major epigenome consortia, and the actual
key downstream consumers join on. We emit the IHEC term as the primary
CURIE and the CL / UBERON / EFO cross-references as ``xrefs`` payload
so consumers that want the upstream identifier can pull it out.

CURIE policy:
- Primary CURIE: ``IHECBio:<symbol>`` — local namespace, declared under
  ``https://w3id.org/holobiomicslab/ihecbio/`` in ``_namespaces.py``.
- Cross-references: only populated where the upstream term is
  unambiguous and stable. K562's ``cl_curie`` is intentionally empty;
  many tissue-level entries lack a ``cl_curie`` because they describe
  the anatomical site, not a single cell type.

Lookup is case-insensitive and tolerates the major notational variants
that appear in the literature: ``K562`` / ``K-562`` / ``K 562``, ``CD4+
T cell`` / ``CD4-positive T cell``, ``HeLa-S3`` / ``HeLa S3``.
"""
from __future__ import annotations

import csv
import re
from importlib.resources import files
from typing import Final

_IHEC_BIOSAMPLE_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "symbol", "aliases", "curie", "label", "type",
    "cl_curie", "uberon_curie", "efo_curie",
)


def _load_ihec_biosample_index_rows() -> list[dict[str, str]]:
    tsv = files("indicium_adapters_epigenomics.data").joinpath(
        "ihec_biosample_index.tsv"
    )
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _IHEC_BIOSAMPLE_INDEX_COLUMNS:
            raise RuntimeError(
                "ihec_biosample_index.tsv header does not match the declared "
                f"schema. Expected {_IHEC_BIOSAMPLE_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [
            {k: (v or "").strip() for k, v in row.items()} for row in reader
        ]


# Biosample symbols and aliases in the wild mix dashes, plus signs and
# spaces freely (``K-562``, ``CD4+ T cell``, ``HeLa S3``). Collapse all
# of those plus the underscore so canonical-form symbols
# (``CD4_T_cell``) still resolve when the input uses the human-readable
# spelling.
_NORMALISE_RE: Final[re.Pattern[str]] = re.compile(r"[\s\-+_]+")


def _normalise(name: str) -> str:
    return _NORMALISE_RE.sub("", name).lower()


class IHECBiosampleIndex:
    """Case-insensitive IHEC biosample lookup with cross-reference payload.

    Pure data; loaded once at import time. Each row is indexed under its
    canonical ``symbol`` plus every pipe-separated alias, all normalised
    by ``_normalise()`` so spelling variants converge on the same row.
    """

    def __init__(self, rows: list[dict[str, str]]):
        self._rows = rows
        self._by_key: dict[str, dict[str, str]] = {}
        for row in rows:
            if symbol := row.get("symbol"):
                self._by_key.setdefault(_normalise(symbol), row)
            for alias in (row.get("aliases") or "").split("|"):
                alias = alias.strip()
                if alias:
                    self._by_key.setdefault(_normalise(alias), row)

    def lookup(self, text: str) -> dict[str, str] | None:
        """Return the matching IHEC biosample row, or ``None``.

        Returns a dict with the standard adapter contract fields
        (``curie``, ``label``, ``type``) plus the populated cross-reference
        fields (``cl_curie``, ``uberon_curie``, ``efo_curie``) — empty
        cross-refs are omitted from the returned dict so a downstream
        consumer doing ``"cl_curie" in match`` gets the right answer.
        """
        if not isinstance(text, str):
            return None
        key = _normalise(text.strip())
        if not key:
            return None
        if row := self._by_key.get(key):
            result: dict[str, str] = {
                "curie": row["curie"],
                "label": row["label"],
                "type": row["type"],
            }
            for xref in ("cl_curie", "uberon_curie", "efo_curie"):
                if value := row.get(xref):
                    result[xref] = value
            return result
        return None

    def __len__(self) -> int:
        return len(self._rows)


IHEC_BIOSAMPLE_INDEX: Final[IHECBiosampleIndex] = IHECBiosampleIndex(
    _load_ihec_biosample_index_rows()
)
