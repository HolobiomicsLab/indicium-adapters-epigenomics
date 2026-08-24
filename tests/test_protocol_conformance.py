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

from indicium_adapters import (  # noqa: E402
    CURRENT_API_VERSION,
    DomainAdapter,
    SHACLProvider,
    discover_adapters,
)

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


def test_every_discovered_adapter_qualifier_is_snake_case():
    """Defence-in-depth: the parametrised test above iterates over the
    hardcoded ``_ADAPTERS`` tuple, which a future sibling could bypass if
    its entry-point is added to ``pyproject.toml`` but its class is
    forgotten in the tuple. This variant goes through the
    ``indicium-adapters`` discovery helper instead, so any adapter
    registered under ``indicium.adapters`` is checked, even if the test
    module has not been updated for it.

    Required by indicium >= 1.7.0: ``asb:domainQualifier`` enforces
    ``^[a-z][a-z0-9_]*$`` at the SHACL boundary; non-conforming qualifiers
    would be silently rejected by ``indicium.validate_graph``.
    """
    discovered = discover_adapters()
    offenders: list[tuple[str, str]] = []
    for name, adapter in discovered.items():
        for q in adapter.qualifiers:
            if not _SNAKE_RE.match(q):
                offenders.append((name, q))
    assert not offenders, (
        "non-snake_case qualifiers on discovered adapters: "
        + ", ".join(f"{name}:{q}" for name, q in offenders)
    )


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


# ---------------------------------------------------------------------------
# Classification vocabulary. ASB unions ``primary_domain`` / ``subdomains`` /
# ``techniques`` from every installed adapter into its controlled terms
# (agentic_science_builder.classification._adapter_vocab); an adapter that
# declares none leaves its whole field unclassifiable.
# ---------------------------------------------------------------------------

_SLUG_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")


@pytest.mark.parametrize(
    "cls",
    [c for c, _, _ in _ADAPTERS],
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_primary_domain_is_epigenomics(cls):
    """All three adapters roll up to one selectable ASB primary domain."""
    assert cls.primary_domain == "epigenomics"


@pytest.mark.parametrize(
    "cls",
    [c for c, _, _ in _ADAPTERS],
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_subdomains_are_nonempty_slugs(cls):
    """Each adapter carries a slice of the audited epigenomics subdomains."""
    assert cls.subdomains, f"{cls.__name__} declares no subdomains"
    for subdomain in cls.subdomains:
        assert _SLUG_RE.match(subdomain), f"{subdomain!r} is not a lowercase slug"


def test_subdomains_are_claimed_by_one_adapter_each():
    """Siblings split the subdomain vocabulary rather than overlapping on it."""
    declared = [s for c, _, _ in _ADAPTERS for s in c.subdomains]
    assert len(declared) == len(set(declared)), "a subdomain is claimed twice"


@pytest.mark.parametrize(
    "cls",
    [c for c, _, _ in _ADAPTERS],
    ids=[c.__name__ for c, _, _ in _ADAPTERS],
)
def test_techniques_are_nonempty_slugs(cls):
    """Non-vacuity plus shape: hyphenated lowercase tokens, like ``chip-seq``."""
    assert cls.techniques, f"{cls.__name__} declares no techniques"
    for technique in cls.techniques:
        assert _SLUG_RE.match(technique), f"{technique!r} is not a lowercase slug"


def test_techniques_partition_the_audited_vocabulary():
    """Every audited token reaches exactly one adapter, and none is invented.

    The 2026-08-25 cross-domain audit, not the resolver table, is the source:
    the table still resolves assay names whose tokens no longer separate the
    sciences that shared them.
    """
    from indicium_adapters_epigenomics._assays import _AUDITED_TECHNIQUES

    declared = [t for c, _, _ in _ADAPTERS for t in c.techniques]
    assert len(declared) == len(set(declared)), "a technique is claimed twice"
    audited = {t for tokens in _AUDITED_TECHNIQUES.values() for t in tokens}
    assert set(declared) == audited


def test_spelled_out_assays_use_their_conventional_short_form():
    """WGBS and RRBS enter the vocabulary as acronyms, not as their long labels."""
    bulk = set(BulkEpigenomicsAdapter.techniques)
    assert {"wgbs", "rrbs"} <= bulk
    assert not [t for t in bulk if t.endswith("-bisulfite-sequencing")]
