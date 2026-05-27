"""Maintainer script — regenerate ``encode_experiment_index.tsv``.

Pulls a snapshot of ENCODE Experiment metadata from the ENCODE portal REST
API (``https://www.encodeproject.org/search/?type=Experiment&format=json``)
and writes one row per experiment with columns:

    accession  assay_type  target  biosample  organism_curie  assembly

Runs offline-only — install with ``pip install
.[build]`` to pull in ``httpx``.

Usage (maintainer-only; not run at install or import time)::

    python -m indicium_adapters_epigenomics.data.build_encode_index \\
        [--limit 5000] [--out src/indicium_adapters_epigenomics/data/encode_experiment_index.tsv]

The bundled ``encode_experiment_index.tsv`` is a small seed (5-10 rows) so
the package wheel stays compact; refresh it before each release.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ENCODE_SEARCH_URL = "https://www.encodeproject.org/search/"
OUT_COLUMNS = (
    "accession", "assay_type", "target", "biosample", "organism_curie",
    "assembly",
)

# Minimal mapping from the NCBITaxon-ID returned by the ENCODE API to the
# CURIE form we write. Add taxa here as the ENCODE corpus broadens.
_TAXON_TO_CURIE: dict[str, str] = {
    "9606": "NCBITaxon:9606",
    "10090": "NCBITaxon:10090",
    "7227": "NCBITaxon:7227",
    "6239": "NCBITaxon:6239",
}


def _fetch_experiments(limit: int) -> list[dict]:
    # Deferred import: the package must import cleanly without the `build`
    # extra. Only the maintainer script needs `httpx`.
    import httpx

    # Restrict the snapshot to assays this adapter actually cares about.
    # Without this filter the ENCODE search returns RNA-seq, miRNA-seq, etc.,
    # which inflates the seed file with rows the adapter has no enrichment for.
    epigenomics_assays = [
        "ChIP-seq", "Histone ChIP-seq", "TF ChIP-seq", "Control ChIP-seq",
        "Mint-ChIP-seq", "Control Mint-ChIP-seq",
        "ATAC-seq", "snATAC-seq",
        "DNase-seq",
        "MNase-seq",
        "WGBS", "RRBS",
        "Hi-C", "intact Hi-C",
        "CUT&RUN", "CUT&Tag",
        "MeDIP-seq", "MRE-seq",
    ]
    params: list[tuple[str, str]] = [
        ("type", "Experiment"),
        ("format", "json"),
        ("limit", str(limit)),
        ("status", "released"),
        ("field", "accession"),
        ("field", "assay_title"),
        ("field", "target.label"),
        ("field", "biosample_ontology.term_name"),
        ("field", "replicates.library.biosample.organism.taxon_id"),
        ("field", "assembly"),
    ]
    params.extend(("assay_title", a) for a in epigenomics_assays)
    print(f"[encode] fetching up to {limit} experiments ...", file=sys.stderr)
    r = httpx.get(ENCODE_SEARCH_URL, params=params, timeout=120.0)
    r.raise_for_status()
    return r.json().get("@graph", [])


def _row_from_experiment(exp: dict) -> dict[str, str]:
    accession = exp.get("accession") or ""
    assay = exp.get("assay_title") or ""
    target_obj = exp.get("target") or {}
    target = target_obj.get("label") or ""
    bio = (exp.get("biosample_ontology") or {}).get("term_name", "")
    assembly_list = exp.get("assembly") or []
    assembly = assembly_list[0] if isinstance(assembly_list, list) and assembly_list else ""
    organism_curie = ""
    for rep in exp.get("replicates") or []:
        taxon_id = (
            ((rep.get("library") or {}).get("biosample") or {})
            .get("organism", {})
            .get("taxon_id")
        )
        if taxon_id and (curie := _TAXON_TO_CURIE.get(str(taxon_id))):
            organism_curie = curie
            break
    return {
        "accession": accession,
        "assay_type": assay,
        "target": target,
        "biosample": bio,
        "organism_curie": organism_curie,
        "assembly": assembly,
    }


def build(out_path: Path, limit: int) -> int:
    experiments = _fetch_experiments(limit)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        # Preserve the seed header comment so consumers can tell the file is
        # auto-generated and never edited by hand.
        fh.write(
            "# Built by build_encode_index.py from https://www.encodeproject.org/ "
            "— refresh before each release.\n"
        )
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t")
        w.writeheader()
        for exp in experiments:
            row = _row_from_experiment(exp)
            if not row["accession"]:
                continue
            w.writerow(row)
            n += 1
    print(f"[out] wrote {out_path}: {n:,} rows", file=sys.stderr)
    return n


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--limit", type=int, default=5000,
        help="Maximum ENCODE experiments to fetch (default: 5000).",
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path(__file__).parent / "encode_experiment_index.tsv",
        help="Output TSV path (default: bundled location).",
    )
    args = parser.parse_args(argv)
    build(args.out, args.limit)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
