"""DNA-modification chemistry → CHEBI / ``EpiMod:`` CURIE resolver.

Split out of ``_histone_marks.py`` in v0.2.0 §3 so the DNA-modification
catalogue is its own concept index — histone PTMs and DNA-base
chemistries are conceptually distinct (different substrate, different
enzymes, different read-out assays), and v0.2.0 introduces methylation-
specific qualifier preferences in ``_resolver.py`` that need a dedicated
source to route to.

CURIE policy:
- Verified upstream CHEBI IDs (5mC ``CHEBI:27551``, 5hmC ``CHEBI:76792``,
  5fC ``CHEBI:76794``) — cross-checked against EBI ChEBI at v0.1.0 release
  and re-validated here. The ``verified`` column records the upstream
  authority (``chebi``).
- Local ``EpiMod:<symbol>`` for entries without a verified upstream term
  (5caC, 6mA at v0.2.0). The IRI ``https://w3id.org/holobiomicslab/epimod/``
  is declared in ``_namespaces.py``. Promotion path: when CHEBI ships a
  matching term, flip the ``curie`` column and the ``verified`` column
  in ``data/dna_modifications.tsv``; downstream consumers that pinned
  against the local IRI keep working until they migrate.

Lookup is case-insensitive and tolerant of the common notational
variants seen in literature: ``5mC``, ``5-mC``, ``5-methylcytosine`` and
``methylcytosine`` all resolve to the same row.
"""
from __future__ import annotations

import csv
import re
from importlib.resources import files
from typing import Final

_DNA_MOD_INDEX_COLUMNS: Final[tuple[str, ...]] = (
    "symbol", "aliases", "curie", "label", "type", "verified",
)


def _load_dna_mod_index_rows() -> list[dict[str, str]]:
    tsv = files("indicium_adapters_epigenomics.data").joinpath(
        "dna_modifications.tsv"
    )
    with tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _DNA_MOD_INDEX_COLUMNS:
            raise RuntimeError(
                "dna_modifications.tsv header does not match the declared "
                f"schema. Expected {_DNA_MOD_INDEX_COLUMNS!r}; got "
                f"{tuple(reader.fieldnames or ())!r}."
            )
        return [
            {k: (v or "").strip() for k, v in row.items()} for row in reader
        ]


# DNA-mod literature uses dashes and parentheses to break up long IUPAC
# names; the normaliser strips whitespace, dashes and parentheses so
# ``5-(hydroxymethyl)cytosine`` and ``5hydroxymethylcytosine`` converge.
_NORMALISE_RE: Final[re.Pattern[str]] = re.compile(r"[\s\-()]+")


def _normalise(name: str) -> str:
    return _NORMALISE_RE.sub("", name).lower()


class DNAModIndex:
    """Case-insensitive DNA-modification chemistry lookup.

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
        """Return the matching DNA-mod row, or ``None``.

        Returns a dict with the standard adapter contract fields
        (``curie``, ``label``, ``type``) plus a ``verified`` field
        recording the provenance authority of the CURIE (``chebi`` for
        verified upstream, ``local`` for ``EpiMod:`` entries).
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
                "verified": row["verified"],
            }
        return None

    def __len__(self) -> int:
        return len(self._rows)


DNA_MOD_INDEX: Final[DNAModIndex] = DNAModIndex(_load_dna_mod_index_rows())
