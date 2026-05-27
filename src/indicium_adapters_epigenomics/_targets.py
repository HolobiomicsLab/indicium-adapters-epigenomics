"""Chromatin-target (TF + chromatin modifier) name → UniProt CURIE resolver.

Loads ``data/chromatin_targets.tsv`` once at import time via
``importlib.resources``. The TSV covers the most commonly assayed chromatin
targets in ENCODE / ChIP-Atlas across human (NCBITaxon:9606) and mouse
(NCBITaxon:10090) — both orthologs are kept as separate rows so the lookup
can disambiguate by ``organism_hint``.

Lookup is case-insensitive on the canonical symbol and pipe-separated
aliases. When two organisms share a symbol (e.g. ``CTCF`` for human and
``Ctcf`` for mouse), the optional ``organism_hint`` parameter selects the
right row; without a hint, the first row in the TSV (which is the human
ortholog by convention) wins.
"""
from __future__ import annotations

import csv
from importlib.resources import files
from typing import Final

_TARGET_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "symbol", "aliases", "curie", "label",
    "organism_curie", "uniprot_id", "category",
)


def _load_target_index_rows() -> list[dict[str, str]]:
    tsv = files("indicium_adapters_epigenomics.data").joinpath(
        "chromatin_targets.tsv"
    )
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _TARGET_INDEX_COLUMNS:
            raise RuntimeError(
                "chromatin_targets.tsv header does not match the declared "
                f"schema. Expected {_TARGET_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


class TargetIndex:
    """Case-insensitive chromatin-target lookup with organism disambiguation.

    Pure data; loaded once at import time. The class only exposes
    ``lookup()``; adapters stay decoupled from the on-disk layout and the
    disambiguation rules.
    """

    def __init__(self, rows: list[dict[str, str]]):
        self._rows = rows
        # Per-symbol list of rows (multiple rows for orthologs).
        self._by_symbol: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            keys = {row["symbol"].lower()} if row.get("symbol") else set()
            for alias in (row.get("aliases") or "").split("|"):
                if alias.strip():
                    keys.add(alias.strip().lower())
            for k in keys:
                self._by_symbol.setdefault(k, []).append(row)

    def lookup(
        self,
        text: str,
        organism_hint: str | None = None,
    ) -> dict[str, str] | None:
        """Return the matching target row, or ``None``.

        Args:
            text: target name (canonical symbol or alias).
            organism_hint: optional ``NCBITaxon:<id>`` CURIE used to pick
                between orthologs sharing a symbol. If omitted, the first
                row in the TSV order wins.
        """
        if not isinstance(text, str):
            return None
        token = text.strip().lower()
        if not token:
            return None
        candidates = self._by_symbol.get(token)
        if not candidates:
            return None
        chosen: dict[str, str] | None = None
        if organism_hint:
            for row in candidates:
                if row.get("organism_curie") == organism_hint:
                    chosen = row
                    break
        if chosen is None:
            chosen = candidates[0]
        return {
            "curie": chosen["curie"],
            "label": chosen["label"],
            "organism_curie": chosen["organism_curie"],
            "uniprot_id": chosen["uniprot_id"],
            "category": chosen["category"],
        }

    def __len__(self) -> int:
        return len(self._rows)


TARGET_INDEX: Final[TargetIndex] = TargetIndex(_load_target_index_rows())
