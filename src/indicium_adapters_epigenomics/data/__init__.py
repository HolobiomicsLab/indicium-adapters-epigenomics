"""Bundled data + maintainer-only build scripts.

Scripts named ``build_*.py`` regenerate the package resources (ENCODE / EpiRR
accession index TSVs, MOD / UniProt chromatin-target vocabulary subsets) by
hitting external APIs. They must NOT run at package install or import time.
See ROADMAP.md.
"""
