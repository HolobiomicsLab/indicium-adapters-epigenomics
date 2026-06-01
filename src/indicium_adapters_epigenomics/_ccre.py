"""ENCODE SCREEN cCRE class-label → local ``EpiCCRE:`` CURIE resolver.

ENCODE SCREEN publishes a 7-class functional classification for
candidate cis-Regulatory Elements (cCREs):

- ``PLS`` — Promoter-Like Signature
- ``pELS`` — proximal Enhancer-Like Signature (within 2 kb of a TSS)
- ``dELS`` — distal Enhancer-Like Signature
- ``CA-CTCF`` — Chromatin Accessible + CTCF binding
- ``CA-H3K4me3`` — Chromatin Accessible + H3K4me3 (off-TSS)
- ``CA-TF`` — Chromatin Accessible + transcription-factor binding
- ``CA-only`` — Chromatin Accessible only

The classes are documented in the ENCODE SCREEN release notes; no
upstream registry mints CURIEs for them, so ``_namespaces.EPICCRE``
hosts the local ``EpiCCRE:<symbol>`` namespace. Individual cCRE
accessions (``EH38E0000001`` for hg38, ``EM10E0000001`` for mm10) are
recognised in ``_registries.py`` and emit under the canonical
``ENCODE:`` prefix — they're ENCODE objects, just with a SCREEN-specific
accession family.

Lookup is case-insensitive and tolerates the major notational variants
seen in literature: ``PLS``, ``pls``, ``CA-CTCF``, ``CA_CTCF``,
``CA CTCF``, ``Promoter-Like Signature`` all converge.
"""
from __future__ import annotations

import csv
import re
from importlib.resources import files
from typing import Final

_CCRE_CLASS_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "symbol", "aliases", "curie", "label", "type",
)


def _load_ccre_class_index_rows() -> list[dict[str, str]]:
    tsv = files("indicium_adapters_epigenomics.data").joinpath(
        "ccre_classes.tsv"
    )
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _CCRE_CLASS_INDEX_COLUMNS:
            raise RuntimeError(
                "ccre_classes.tsv header does not match the declared schema. "
                f"Expected {_CCRE_CLASS_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [
            {k: (v or "").strip() for k, v in row.items()} for row in reader
        ]


# Normalisation strips whitespace, dashes, underscores and plus signs so
# ``CA-CTCF``, ``CA_CTCF``, ``CA CTCF`` and ``CA+CTCF`` converge. Mixed
# case folds to lowercase.
_NORMALISE_RE: Final[re.Pattern[str]] = re.compile(r"[\s\-+_]+")


def _normalise(name: str) -> str:
    return _NORMALISE_RE.sub("", name).lower()


class CCREClassIndex:
    """Case-insensitive cCRE class-label lookup.

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
        """Return the matching cCRE class row, or ``None``."""
        if not isinstance(text, str):
            return None
        key = _normalise(text.strip())
        if not key:
            return None
        if row := self._by_key.get(key):
            return {
                "curie": row["curie"],
                "label": row["label"],
                "type": row["type"],
            }
        return None

    def __len__(self) -> int:
        return len(self._rows)


CCRE_CLASS_INDEX: Final[CCREClassIndex] = CCREClassIndex(
    _load_ccre_class_index_rows()
)
