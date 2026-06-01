"""Behavioural tests for the shared qualifier-aware resolver.

The bundled v0.1.0 indices have no token that genuinely cross-hits two
domains (verified via direct lookup against ``HISTONE_INDEX`` and
``TARGET_INDEX``), so demonstrating the qualifier preference requires
synthetic cross-hits installed via ``unittest.mock``. The infrastructure
is exercised below; the real-data ambiguities it guards against arrive
with the v0.2.0 §2-§6 indices (ChromHMM states, IHEC biosamples, ENCODE
cCRE).

Two regression invariants are also pinned:

1. Accession always wins — even when a qualifier preference would point
   elsewhere.
2. With no qualifier hint, the resolver behaves byte-identically to the
   v0.1.0 ``_resolve(value)`` fixed cascade.
"""
from __future__ import annotations

from unittest.mock import patch

import pytest

from indicium_adapters_epigenomics import _resolver
from indicium_adapters_epigenomics._resolver import resolve

_HISTONE_FAKE = {"curie": "EpiHist:FAKE", "label": "fake histone PTM"}
_TARGET_FAKE = {"curie": "UniProt:FAKE", "label": "fake chromatin target"}


# --- Backwards-compatibility: no qualifier → v0.1.0 fixed-cascade order ---


@pytest.mark.parametrize(
    "value, expected_source",
    [
        ("ENCSR000AKB", "accession"),  # ENCODE accession → accession.
        ("ATAC-seq", "assay"),  # known assay.
        ("GRCh38", "assembly"),  # known assembly.
        ("H3K27ac", "histone"),  # known histone PTM.
        ("EZH2", "target"),  # known chromatin target.
        ("FRiP", "qc"),  # known QC metric.
    ],
)
def test_no_qualifier_reproduces_v0_1_0_cascade(value, expected_source):
    """Without a qualifier hint, each real-data lookup hits the same
    source it did under the v0.1.0 ``_resolve(value)`` cascade. This is
    the regression-guard for the resolver refactor — any drift here
    would silently change every adapter's enrichment output.
    """
    match, source = resolve(value)
    assert match is not None, f"{value!r} should resolve in v0.1.0 cascade"
    assert source == expected_source


def test_unknown_value_returns_none_pair():
    match, source = resolve("definitely-not-a-known-token-1234")
    assert match is None
    assert source == ""


def test_unknown_qualifier_falls_back_to_default_cascade():
    """A qualifier with no entry in ``_QUALIFIER_PREFERENCES`` must not
    change resolution — preference is opt-in, not required."""
    match, source = resolve("EZH2", qualifier="not_a_registered_qualifier")
    assert match is not None
    assert source == "target"


# --- Accession always wins, regardless of qualifier hint ---


def test_accession_outranks_qualifier_preference():
    """ENCSR000AKB is an accession AND ``bound_by`` would prefer
    ``TARGET_INDEX`` first — but accessions outrank the preference."""
    match, source = resolve("ENCSR000AKB", qualifier="bound_by")
    assert match is not None
    assert source == "accession"


# --- Qualifier preference observed when indices cross-hit ---


def test_bound_by_prefers_target_over_histone_on_synthetic_overlap():
    """Inject a synthetic match on both HISTONE and TARGET for the same
    token. Under ``bound_by`` the resolver must return the TARGET hit
    even though HISTONE comes first in the v0.1.0 default order.
    """
    with patch.object(
        _resolver.HISTONE_INDEX, "lookup", return_value=_HISTONE_FAKE
    ), patch.object(
        _resolver.TARGET_INDEX, "lookup", return_value=_TARGET_FAKE
    ):
        match, source = resolve("AMBIGUOUS_TOKEN", qualifier="bound_by")

    assert match == _TARGET_FAKE
    assert source == "target"


def test_marked_by_prefers_histone_over_target_on_synthetic_overlap():
    """The inverse: under ``marked_by`` the resolver must prefer
    HISTONE over TARGET. (HISTONE already wins in the v0.1.0 cascade,
    but this pins the preference contract explicitly.)
    """
    with patch.object(
        _resolver.HISTONE_INDEX, "lookup", return_value=_HISTONE_FAKE
    ), patch.object(
        _resolver.TARGET_INDEX, "lookup", return_value=_TARGET_FAKE
    ):
        match, source = resolve("AMBIGUOUS_TOKEN", qualifier="marked_by")

    assert match == _HISTONE_FAKE
    assert source == "histone"


def test_assay_preference_shadows_target_under_assayed_with():
    """``assayed_with`` lists ASSAY first; if a token cross-hits ASSAY
    and TARGET, ASSAY must win. Cross-cuts an enzyme that's also an
    assay-protocol name (e.g. ``MNase``) — synthetic for now.
    """
    fake_assay = {"curie": "OBI:FAKE", "label": "fake assay"}
    fake_target = {"curie": "UniProt:FAKETGT", "label": "fake target"}
    with patch.object(
        _resolver.ASSAY_INDEX, "lookup", return_value=fake_assay
    ), patch.object(
        _resolver.TARGET_INDEX, "lookup", return_value=fake_target
    ):
        match, source = resolve("MNase", qualifier="assayed_with")

    assert match == fake_assay
    assert source == "assay"


# --- Preference is a hint, not a route ---


def test_preferred_source_misses_then_falls_through_to_default_order():
    """If the preferred source returns ``None`` for the value, the
    resolver must still try the remaining sources in default order so
    backwards-compatibility holds for any token whose preferred index
    just doesn't contain it.
    """
    # TARGET_INDEX returns None for this token (which is what really
    # happens for "H3K27ac" — a histone, not a target), but the
    # qualifier preference for "bound_by" is (TARGET, HISTONE). The
    # resolver should fall through to HISTONE rather than giving up.
    match, source = resolve("H3K27ac", qualifier="bound_by")
    assert match is not None
    assert source == "histone"


# --- Snake_case keyspace contract on the preference map itself ---


def test_every_qualifier_preference_key_is_snake_case():
    """The keyspace of ``_QUALIFIER_PREFERENCES`` is the union of every
    adapter's qualifier vocabulary, so every key must satisfy the same
    ``^[a-z][a-z0-9_]*$`` contract that
    ``test_protocol_conformance.test_every_discovered_adapter_qualifier_is_snake_case``
    pins for the public qualifiers. Otherwise the map would silently
    ignore the intended qualifier at lookup time (key mismatch).
    """
    import re

    snake = re.compile(r"^[a-z][a-z0-9_]*$")
    offenders = [
        q for q in _resolver._QUALIFIER_PREFERENCES if not snake.match(q)
    ]
    assert not offenders, f"non-snake_case preference keys: {offenders}"


def test_every_preference_value_is_a_known_source_label():
    """Each entry in a preference tuple must be one of the canonical
    Source constants — typos would silently disable that preference.
    """
    known = {
        _resolver.ACCESSION,
        _resolver.ASSAY,
        _resolver.ASSEMBLY,
        _resolver.HISTONE,
        _resolver.TARGET,
        _resolver.QC,
        _resolver.CHROM_STATE,
        _resolver.DNA_MOD,
        _resolver.IHEC_BIOSAMPLE,
        _resolver.CCRE,
        _resolver.CONFORMATION,
    }
    offenders: list[tuple[str, str]] = []
    for qualifier, sources in _resolver._QUALIFIER_PREFERENCES.items():
        for s in sources:
            if s not in known:
                offenders.append((qualifier, s))
    assert not offenders, f"unknown source labels in preference map: {offenders}"
