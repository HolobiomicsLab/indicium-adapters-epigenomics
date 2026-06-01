"""3D-genome compartment label → local ``EpiConf:`` CURIE resolver.

Hi-C / Micro-C analyses partition the genome into A / B compartments
(Lieberman-Aiden et al. 2009) and the Rao et al. (2014) refinement into
sub-compartments A1 / A2 / B1 / B2 / B3. Other 3D-genome semantics —
contacts, loops, anchors, TAD boundaries, sub-TADs — do not have a
controlled vocabulary of labels; they're typically referenced by
coordinate or anchor ID. Only the (sub-)compartment labels benefit from
a CURIE-resolving index.

CURIE policy:
- Local ``EpiConf:<symbol>`` for every compartment / sub-compartment
  label. No upstream ontology mints CURIEs for these terms at v0.2.0.

Lookup is case-insensitive and tolerates the common spellings: ``A``,
``Compartment A``, ``A_compartment``, ``compartment_A`` all converge.
Sub-compartments accept ``A1``, ``subA1``, ``A1_sub`` and
``A1 subcompartment``.
"""
from __future__ import annotations

import csv
import re
from importlib.resources import files
from typing import Final

_CONFORMATION_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "symbol", "aliases", "curie", "label", "type",
)


def _load_conformation_index_rows() -> list[dict[str, str]]:
    tsv = files("indicium_adapters_epigenomics.data").joinpath(
        "conformation_labels.tsv"
    )
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _CONFORMATION_INDEX_COLUMNS:
            raise RuntimeError(
                "conformation_labels.tsv header does not match the declared "
                f"schema. Expected {_CONFORMATION_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [
            {k: (v or "").strip() for k, v in row.items()} for row in reader
        ]


# Strip whitespace, dashes, underscores so ``A``, ``A_compartment``,
# ``Compartment A`` all converge to ``acompartment`` / ``compartmenta``.
# Note: the bare label ``A`` is the canonical form; longer spellings
# are aliases.
_NORMALISE_RE: Final[re.Pattern[str]] = re.compile(r"[\s\-_]+")


def _normalise(name: str) -> str:
    return _NORMALISE_RE.sub("", name).lower()


class ConformationIndex:
    """Case-insensitive A/B compartment + sub-compartment lookup.

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
        """Return the matching compartment row, or ``None``."""
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


CONFORMATION_INDEX: Final[ConformationIndex] = ConformationIndex(
    _load_conformation_index_rows()
)
