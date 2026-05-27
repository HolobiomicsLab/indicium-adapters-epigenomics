"""SHACL shapes — load + smoke-validation tests.

All three sibling adapters return the same ``shapes/epigenomics.shapes.ttl``
graph from ``shacl_shapes()``. The TTL provides per-qualifier CURIE-prefix
constraints for the ASB claim slots emitted by each adapter.

The negative test deliberately targets ``epi:marked_by`` (severity
``sh:Violation``), not ``epi:localized_in_spatial_domain`` /
``epi:clock_age_estimated_with`` (severity ``sh:Warning``) — those would
report a warning but still set ``conforms=True``.
"""
from __future__ import annotations

import pytest

pyshacl = pytest.importorskip("pyshacl")
import rdflib  # noqa: E402

from indicium_adapters_epigenomics import (  # noqa: E402
    BulkEpigenomicsAdapter,
    MethylationArrayEpigenomicsAdapter,
    SingleCellEpigenomicsAdapter,
)

# Turtle preamble shared by every data-graph fixture.
_PREFIXES = """
@prefix asb: <http://holobiomicslab.org/asb/> .
@prefix epi: <http://holobiomicslab.org/adapters/epigenomics/> .
"""


def _parse(data_ttl: str) -> rdflib.Graph:
    g = rdflib.Graph()
    g.parse(data=_PREFIXES + data_ttl, format="turtle")
    return g


# --- graph-level smoke tests ------------------------------------------------


def test_shacl_graph_is_nonempty_and_parses():
    g = BulkEpigenomicsAdapter().shacl_shapes()
    assert isinstance(g, rdflib.Graph)
    assert len(g) > 0


def test_all_three_adapters_return_the_same_shapes_graph_size():
    """All three sibling adapters share a single TTL; the per-adapter
    ``shacl_shapes()`` calls must therefore return graphs of identical size."""
    sizes = {
        BulkEpigenomicsAdapter().shacl_shapes().__len__(),
        SingleCellEpigenomicsAdapter().shacl_shapes().__len__(),
        MethylationArrayEpigenomicsAdapter().shacl_shapes().__len__(),
    }
    assert len(sizes) == 1, f"divergent shapes-graph sizes: {sizes}"


# --- empty-graph baseline ---------------------------------------------------


def test_validate_empty_data_graph_conforms():
    """An empty data graph has nothing to violate — pyshacl must report
    ``conforms=True`` without raising."""
    shapes = BulkEpigenomicsAdapter().shacl_shapes()
    conforms, _report_graph, _report_text = pyshacl.validate(
        data_graph=rdflib.Graph(),
        shacl_graph=shapes,
    )
    assert conforms is True


# --- positive case: well-formed claim conforms ------------------------------


def test_validate_well_formed_marked_by_claim_conforms():
    """A claim whose ``marked_by`` carries an EpiHist: or MOD: CURIE must
    conform. The v0.1.0 histone vocabulary ships under the local ``EpiHist:``
    namespace pending upstream PR / MOD coverage of position-specific
    histone PTMs."""
    shapes = BulkEpigenomicsAdapter().shacl_shapes()
    data = _parse(
        '<urn:c1> a asb:Claim ; '
        'epi:marked_by "EpiHist:H3K27ac" .'
    )
    conforms, _report_graph, _report_text = pyshacl.validate(
        data_graph=data,
        shacl_graph=shapes,
    )
    assert conforms is True


# --- negative case: malformed CURIE on a Violation-severity qualifier ------


def test_validate_malformed_marked_by_claim_violates():
    """``epi:marked_by`` is a ``sh:Violation`` constraint — a non-MOD / PR /
    CHEBI literal must drive ``conforms=False`` and surface the qualifier
    name in the report text."""
    shapes = BulkEpigenomicsAdapter().shacl_shapes()
    data = _parse(
        '<urn:c2> a asb:Claim ; '
        'epi:marked_by "random_string_not_a_curie" .'
    )
    conforms, _report_graph, report_text = pyshacl.validate(
        data_graph=data,
        shacl_graph=shapes,
    )
    assert conforms is False
    assert "marked_by" in report_text
