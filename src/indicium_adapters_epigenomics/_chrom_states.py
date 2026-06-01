"""Chromatin-state label → local ``EpiState:`` CURIE resolver.

Covers the canonical ChromHMM 15-state, 18-state and 25-state-imputed
segmentation labels (Roadmap Epigenomics + ENCODE3) plus the commonly
interpreted Segway labels. Both projects publish their state vocabularies
as documentation tables rather than queryable APIs; no upstream OBO
ontology has minted CURIEs for individual segmentation labels at the time
of writing (2026-06).

CURIE policy in v0.2.0:
- Every state resolves to a local ``EpiState:<symbol>`` CURIE, declared
  under the ``https://w3id.org/holobiomicslab/epistate/`` IRI in
  ``_namespaces.py``.
- Each row carries a ``models`` field (pipe-separated) recording which
  segmentation models include the state — preserved verbatim into the
  ``xrefs`` payload so downstream consumers can disambiguate
  shared-label collisions (e.g. ``TssA`` appears in all three ChromHMM
  models with slightly different boundary semantics).
- Promotion path: if SO or GO ever mint sequence-feature CURIEs for
  ChromHMM/Segway states (or if Roadmap publishes its own term IRIs),
  swap the ``curie`` column and add a row in ROADMAP §"CURIE provenance
  policy" — the local IRI then survives as an ``xrefs.legacy_local_iri``
  marker.

Lookup tolerates the major notational variants seen in literature:
``TssA``, ``tssa``, ``Active_TSS``, ``ActiveTSS`` and ``Active TSS`` all
collapse to the same key, as do ``ZNF/Rpts``, ``ZNF_Rpts``, ``ZnfRpts``
and ``ZNFRepeats``.
"""
from __future__ import annotations

import csv
import re
from importlib.resources import files
from typing import Final

_CHROM_STATE_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "symbol", "aliases", "curie", "label", "type", "models",
)


def _load_chrom_state_index_rows() -> list[dict[str, str]]:
    tsv = files("indicium_adapters_epigenomics.data").joinpath(
        "chromatin_states.tsv"
    )
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _CHROM_STATE_INDEX_COLUMNS:
            raise RuntimeError(
                "chromatin_states.tsv header does not match the declared "
                f"schema. Expected {_CHROM_STATE_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [
            {k: (v or "").strip() for k, v in row.items()} for row in reader
        ]


# Normalisation collapses whitespace, dashes, underscores, slashes and
# apostrophes so ``Active TSS``, ``Active_TSS``, ``ActiveTSS``, ``Tx5'``
# and ``Tx 5`` all converge on the same lookup key. This is broader than
# the histone-mark normaliser because segmentation labels in the wild have
# more spelling variance (``ZNF/Rpts`` vs ``ZNF_Rpts`` vs ``ZnfRpts``).
_NORMALISE_RE: Final[re.Pattern[str]] = re.compile(r"[\s\-/_']+")


def _normalise(name: str) -> str:
    return _NORMALISE_RE.sub("", name).lower()


class ChromStateIndex:
    """Case-insensitive chromatin-state-label lookup.

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
        """Return the matching chromatin-state row, or ``None``.

        Match is exact-after-normalisation. Returns a dict with the
        standard adapter contract fields (``curie``, ``label``, ``type``)
        plus a ``models`` field listing which segmentation models the
        state belongs to.
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
                "models": row["models"],
            }
        return None

    def __len__(self) -> int:
        return len(self._rows)


CHROM_STATE_INDEX: Final[ChromStateIndex] = ChromStateIndex(
    _load_chrom_state_index_rows()
)
