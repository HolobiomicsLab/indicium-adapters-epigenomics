"""Histone-modification name → EpiHist / UniProt / CHEBI CURIE resolver.

Loads ``data/histone_marks.tsv`` once at import time via
``importlib.resources``. The TSV ships ~50 canonical marks plus the major
histone variants (H2A.Z, H2A.X / γH2AX, H3.3, macroH2A, CENP-A) and the
core DNA-modification chemistries (5mC, 5hmC, 5fC, 5caC, 6mA).

Lookup is case-insensitive and tolerant of common notational variants —
``H3K27Ac``, ``H3K27ac``, ``H3 K27 ac`` and ``H3K27-ac`` all resolve to the
same row.

CURIE choice policy in v0.1.0:
- Histone PTM symbols (e.g. H3K27ac) → local ``EpiHist:`` vocabulary. PSI-MOD
  (``MOD:``) coverage of position-specific histone PTMs is incomplete (many
  conventional IDs in the literature point at unrelated PSI-MOD terms or have
  no upstream PSI-MOD entry at all). The local vocabulary is documented and
  promoted via the ROADMAP — see also build_mod_index.py.
- Histone variants → UniProt (variants are protein isoforms, not modifications).
- DNA modifications → CHEBI when verified; otherwise local ``EpiHist:`` (5mC
  is verified as CHEBI:27551; 5hmC/5fC/5caC/6mA stay local pending verified
  ChEBI IDs).
"""
from __future__ import annotations

import csv
import re
from importlib.resources import files
from typing import Final

_HISTONE_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "symbol", "aliases", "curie", "label", "type",
)


def _load_histone_index_rows() -> list[dict[str, str]]:
    tsv = files("indicium_adapters_epigenomics.data").joinpath("histone_marks.tsv")
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _HISTONE_INDEX_COLUMNS:
            raise RuntimeError(
                "histone_marks.tsv header does not match the declared schema. "
                f"Expected {_HISTONE_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


# Normalisation: strip internal whitespace and dashes so ``H3 K27 ac``,
# ``H3K27-ac`` and ``H3K27ac`` all collapse to a single lookup key.
_NORMALISE_RE: Final[re.Pattern[str]] = re.compile(r"[\s\-]+")


def _normalise(name: str) -> str:
    return _NORMALISE_RE.sub("", name).lower()


class HistoneIndex:
    """Case-insensitive histone-mark / variant / DNA-modification lookup.

    Pure data; loaded once at import time. Each row is indexed under its
    canonical ``symbol`` plus every pipe-separated alias, all normalised by
    ``_normalise()`` so spelling variants converge on the same row.
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
        """Return the matching histone-mark row, or ``None``.

        Match is exact-after-normalisation: a full mark/variant name is
        expected, not a free-text snippet. Use ``re.search`` upstream if
        scanning prose.
        """
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


HISTONE_INDEX: Final[HistoneIndex] = HistoneIndex(_load_histone_index_rows())
