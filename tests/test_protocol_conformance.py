"""Protocol-conformance tests for the three epigenomics adapter classes.

Every adapter shipped in this wheel must satisfy both ``DomainAdapter`` and
``SHACLProvider`` Protocols, expose a stable ``domain_id`` / ``api_version``,
a non-empty lowercase-snake_case qualifier frozenset, an IRI-valued ontology-
prefix map, a domain-naming ``extraction_context()`` string and a non-empty
SHACL graph from ``shacl_shapes()``.
"""
from __future__ import annotations

import re

import pytest
import rdflib

pytest.importorskip("indicium_adapters")

from indicium_adapters import CURRENT_API_VERSION, DomainAdapter, SHACLProvider  # noqa: E402

from indicium_adapters_epigenomics import (  # noqa: E402
    BulkEpigenomicsAdapter,
    MethylationArrayEpigenomicsAdapter,
    SingleCellEpigenomicsAdapter,
)

# Lower-snake_case identifier check (matches the convention used by the
# transcriptomics package's qualifier vocabulary).
_SNAKE_RE = re.compile(r"^[a-z][a-z0-9_]*$")

# Table of (adapter_class, expected_domain_id, expected_domain_phrase).
# ``expected_domain_phrase`` is a substring that ``extraction_context()`` is
# guaranteed to contain (case-insensitive).
_ADAPTERS = [
    (BulkEpigenomicsAdapter, "bulk-epigenomics", "bulk epigenomics"),
    (
        SingleCellEpigenomicsAdapter,
        "single-cell-epigenomics",
        "single-cell epigenomics",
    ),
    (
        MethylationArrayEpigenomicsAdapter,
        "methylation-array-epigenomics",
        "methylation arrays",
    ),
]


@pytest.mark.parametrize(
    "cls",
    [c for c, _, _ in _ADAPTERS],
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_implements_domain_adapter_protocol(cls):
    assert isinstance(cls(), DomainAdapter)


@pytest.mark.parametrize(
    "cls",
    [c for c, _, _ in _ADAPTERS],
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_implements_shacl_provider_protocol(cls):
    assert isinstance(cls(), SHACLProvider)


@pytest.mark.parametrize(
    "cls",
    [c for c, _, _ in _ADAPTERS],
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_api_version_matches_core(cls):
    assert cls.api_version == CURRENT_API_VERSION == 1


@pytest.mark.parametrize(
    "cls, expected_domain_id",
    [(c, d) for c, d, _ in _ADAPTERS],
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_domain_id_is_canonical(cls, expected_domain_id):
    assert cls.domain_id == expected_domain_id


@pytest.mark.parametrize(
    "cls",
    [c for c, _, _ in _ADAPTERS],
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_qualifiers_are_a_nonempty_lowercase_snake_case_frozenset(cls):
    qs = cls.qualifiers
    assert isinstance(qs, frozenset)
    assert len(qs) > 0
    for q in qs:
        assert isinstance(q, str)
        assert _SNAKE_RE.match(q), f"qualifier {q!r} is not lowercase snake_case"


@pytest.mark.parametrize(
    "cls",
    [c for c, _, _ in _ADAPTERS],
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_ontology_prefixes_is_dict_with_iri_values(cls):
    p = cls.ontology_prefixes
    assert isinstance(p, dict)
    # Domain spec requires at least 25 prefixes — the package currently
    # ships 30.
    assert len(p) >= 25
    for prefix, iri in p.items():
        assert isinstance(prefix, str) and prefix
        assert isinstance(iri, str)
        assert iri.startswith(("http://", "https://")), (
            f"prefix {prefix!r} maps to non-HTTP(S) IRI {iri!r}"
        )


@pytest.mark.parametrize(
    "cls, _domain_id, expected_phrase",
    _ADAPTERS,
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_extraction_context_is_non_empty_and_names_the_domain(
    cls, _domain_id, expected_phrase
):
    ctx = cls().extraction_context()
    assert isinstance(ctx, str)
    assert ctx.strip()
    assert expected_phrase.lower() in ctx.lower()


@pytest.mark.parametrize(
    "cls",
    [c for c, _, _ in _ADAPTERS],
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_shacl_shapes_returns_a_nonempty_parsed_graph(cls):
    g = cls().shacl_shapes()
    assert isinstance(g, rdflib.Graph)
    assert len(g) > 0
