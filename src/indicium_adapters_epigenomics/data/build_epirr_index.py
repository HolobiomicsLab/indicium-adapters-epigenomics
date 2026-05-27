"""Maintainer script — regenerate ``epirr_reference_index.tsv``.

Originally targeted ``https://www.ebi.ac.uk/vg/epirr/view/all/json``, but
EBI retired the EpiRR JSON API in 2024-2025 and the registry web-front-end
now serves HTML only. We fall back to ENCODE's ``ReferenceEpigenome`` REST
endpoint, which federates IHEC reference epigenomes (Roadmap, ENCODE, and
the IHEC consortia that submitted to ENCODE) and still returns JSON. The
coverage is narrower than the original EpiRR dump (no Blueprint / DEEP /
KNIH / EpiHK / AMED-CREST records that never made it into ENCODE), but the
columns map cleanly and the data is live.

Columns written:

    epirr_id  assay  biosample  organism_curie  consortium

Usage (maintainer-only; not run at install or import time)::

    python -m indicium_adapters_epigenomics.data.build_epirr_index \\
        [--limit 5000] [--out src/indicium_adapters_epigenomics/data/epirr_reference_index.tsv]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ENCODE_REFEPI_URL = "https://www.encodeproject.org/search/"
OUT_COLUMNS = (
    "epirr_id", "assay", "biosample", "organism_curie", "consortium",
)

_TAXON_TO_CURIE: dict[str, str] = {
    "9606": "NCBITaxon:9606",
    "10090": "NCBITaxon:10090",
}


def _fetch_reference_epigenomes(limit: int) -> list[dict]:
    # Deferred import: package must import cleanly without the `build` extra.
    import httpx

    params: list[tuple[str, str]] = [
        ("type", "ReferenceEpigenome"),
        ("format", "json"),
        ("limit", str(limit)),
        ("status", "released"),
        ("field", "accession"),
        ("field", "dbxrefs"),
        ("field", "assay_term_name"),
        ("field", "biosample_ontology.term_name"),
        ("field", "organism.taxon_id"),
        ("field", "award.project"),
        ("field", "related_datasets.assay_title"),
        ("field", "related_datasets.biosample_ontology.term_name"),
    ]
    print(f"[epirr] fetching {ENCODE_REFEPI_URL} (ReferenceEpigenome) ...", file=sys.stderr)
    r = httpx.get(ENCODE_REFEPI_URL, params=params, timeout=180.0)
    r.raise_for_status()
    return r.json().get("@graph", [])


def _epirr_id_from(rec: dict) -> str:
    """Prefer the IHEC EpiRR accession (in ``dbxrefs``); fall back to the
    ENCODE ReferenceEpigenome accession."""
    for x in rec.get("dbxrefs", []) or []:
        if isinstance(x, str) and x.startswith("IHEC:"):
            return x.removeprefix("IHEC:")
    return rec.get("accession") or ""


def _first_term_name(bo) -> str:
    """ReferenceEpigenome's ``biosample_ontology`` slot can be a single dict
    OR a list of dicts (one per related dataset). Handle both."""
    if isinstance(bo, dict):
        return bo.get("term_name", "") or ""
    if isinstance(bo, list):
        for it in bo:
            if isinstance(it, dict):
                tn = it.get("term_name", "") or ""
                if tn:
                    return tn
    return ""


def _row_from_record(rec: dict) -> dict[str, str]:
    related = rec.get("related_datasets") or []
    assays = sorted({
        (ds.get("assay_title") or "")
        for ds in related
        if isinstance(ds, dict) and ds.get("assay_title")
    })
    bios = _first_term_name(rec.get("biosample_ontology"))
    if not bios:
        for ds in related:
            if isinstance(ds, dict):
                bios = _first_term_name(ds.get("biosample_ontology"))
                if bios:
                    break
    organism = rec.get("organism")
    if isinstance(organism, list):
        organism = organism[0] if organism and isinstance(organism[0], dict) else {}
    elif not isinstance(organism, dict):
        organism = {}
    taxon_id = str(organism.get("taxon_id", "") or "")
    consortium = ""
    award = rec.get("award")
    if isinstance(award, dict):
        consortium = award.get("project") or ""
    return {
        "epirr_id": _epirr_id_from(rec),
        "assay": "|".join(assays),
        "biosample": bios,
        "organism_curie": _TAXON_TO_CURIE.get(taxon_id, ""),
        "consortium": consortium,
    }


def build(out_path: Path, limit: int) -> int:
    records = _fetch_reference_epigenomes(limit)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write(
            "# Built by build_epirr_index.py — sourced from the ENCODE "
            "ReferenceEpigenome endpoint after EBI retired the EpiRR JSON "
            "API. Refresh before each release.\n"
        )
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t")
        w.writeheader()
        for rec in records:
            row = _row_from_record(rec)
            if not row["epirr_id"]:
                continue
            w.writerow(row)
            n += 1
    print(f"[out] wrote {out_path}: {n:,} rows", file=sys.stderr)
    return n


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--limit", type=int, default=5000,
        help="Maximum reference epigenomes to fetch (default: 5000).",
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).parent / "epirr_reference_index.tsv",
        help="Output TSV path (default: bundled location).",
    )
    args = parser.parse_args(argv)
    build(args.out, args.limit)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
