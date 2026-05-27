"""Maintainer script — regenerate ``histone_marks.tsv``.

Queries the EBI OLS4 API for the PSI-MOD ontology (``mod``) and extracts
the descendants of *histone modified residue* (MOD:01156). Writes one row
per term with columns:

    symbol  aliases  curie  label  type

Notes for the curator:

- OLS4 returns the canonical PSI-MOD label (e.g. *"trimethylated L-lysine"*)
  which is not the symbol used in publications (``H3K4me3``). The script
  prints a draft TSV but ALIASES MUST BE MANUALLY CURATED — the symbol-form
  derivation from the parent histone + position is not encoded in PSI-MOD.
- DNA modifications (5mC, 5hmC, 5fC, 5caC, 6mA) come from CHEBI, not
  PSI-MOD, and must be appended by hand or by a separate CHEBI query.
- Histone variants (H2A.Z, H2A.X, H3.3, macroH2A, CENP-A) come from UniProt
  reviewed entries (per-organism) and must be appended likewise.

Usage (maintainer-only; not run at install or import time)::

    python -m indicium_adapters_epigenomics.data.build_mod_index \\
        [--out src/indicium_adapters_epigenomics/data/histone_marks.tsv]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

OLS_DESCENDANTS_URL = (
    "https://www.ebi.ac.uk/ols4/api/ontologies/mod/terms/"
    "http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FMOD_01156/descendants"
)
OUT_COLUMNS = ("symbol", "aliases", "curie", "label", "type")


def _classify(label: str) -> str:
    """Map a PSI-MOD label to one of the project's ``type`` buckets."""
    lab = label.lower()
    if "lysine" in lab and ("methyl" in lab or "methylated" in lab):
        return "lysine_methylation"
    if "lysine" in lab and ("acetyl" in lab or "acetylated" in lab):
        return "lysine_acetylation"
    if "arginine" in lab and ("methyl" in lab or "methylated" in lab):
        return "arginine_methylation"
    if "serine" in lab and ("phospho" in lab or "phosphorylated" in lab):
        return "serine_phosphorylation"
    return "lysine_methylation"  # conservative fallback


def _fetch_descendants() -> list[dict]:
    # Deferred import: package must import cleanly without the `build` extra.
    import httpx

    print(f"[ols4] fetching {OLS_DESCENDANTS_URL} ...", file=sys.stderr)
    r = httpx.get(OLS_DESCENDANTS_URL, params={"size": 500}, timeout=120.0)
    r.raise_for_status()
    return (r.json().get("_embedded", {}) or {}).get("terms", [])


def build(out_path: Path) -> int:
    terms = _fetch_descendants()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        # Stay header-compatible with the runtime loader; that loader does
        # not tolerate a comment line, so this script writes the schema
        # header directly. The aliases column is left blank — manual
        # curation step REQUIRED before promoting to the bundled TSV.
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t")
        w.writeheader()
        for term in terms:
            curie = (term.get("obo_id") or "").strip()
            label = (term.get("label") or "").strip()
            if not curie.startswith("MOD:") or not label:
                continue
            w.writerow({
                "symbol": label,  # placeholder; curate to H3K?me? symbol
                "aliases": "",
                "curie": curie,
                "label": label,
                "type": _classify(label),
            })
            n += 1
    print(
        f"[out] wrote {out_path}: {n:,} rows. "
        "REMINDER: aliases column needs manual curation; "
        "DNA modifications (CHEBI) and histone variants (UniProt) must be "
        "appended by hand.",
        file=sys.stderr,
    )
    return n


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).parent / "histone_marks.tsv",
        help="Output TSV path (default: bundled location).",
    )
    args = parser.parse_args(argv)
    build(args.out)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
