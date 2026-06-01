"""Cross-adapter composition invariants for the three epigenomics adapters.

The three sibling adapters in this wheel MUST have pairwise-disjoint qualifier
sets so a composed adapter's qualifier vocabulary is a clean union with no
shadowing. They MUST also expose distinct ``domain_id``s and be discoverable
via the ``indicium.adapters`` entry-point group.
"""
from __future__ import annotations

import pytest

pytest.importorskip("indicium_adapters")

from indicium_adapters import DomainAdapter, discover_adapters  # noqa: E402

from indicium_adapters_epigenomics import (  # noqa: E402
    BulkEpigenomicsAdapter,
    MethylationArrayEpigenomicsAdapter,
    SingleCellEpigenomicsAdapter,
)

# --- discovery / entry-point registration -----------------------------------


def test_all_three_adapters_discoverable_via_entry_points():
    found = discover_adapters()
    for expected in (
        "bulk-epigenomics",
        "single-cell-epigenomics",
        "methylation-array-epigenomics",
    ):
        assert expected in found, (
            f"{expected!r} not registered under the indicium.adapters "
            f"entry-point group; got keys {sorted(found)}"
        )


def test_discovered_adapters_implement_domain_adapter_protocol():
    found = discover_adapters()
    for name in (
        "bulk-epigenomics",
        "single-cell-epigenomics",
        "methylation-array-epigenomics",
    ):
        adapter = found[name]
        assert isinstance(adapter, DomainAdapter), (
            f"discovered adapter {name!r} does not satisfy DomainAdapter"
        )


# --- domain_id distinctness -------------------------------------------------


def test_domain_ids_are_distinct():
    ids = {
        BulkEpigenomicsAdapter.domain_id,
        SingleCellEpigenomicsAdapter.domain_id,
        MethylationArrayEpigenomicsAdapter.domain_id,
    }
    assert len(ids) == 3


# --- pairwise qualifier disjointness ----------------------------------------


@pytest.mark.parametrize(
    "left, right",
    [
        (BulkEpigenomicsAdapter, SingleCellEpigenomicsAdapter),
        (BulkEpigenomicsAdapter, MethylationArrayEpigenomicsAdapter),
        (SingleCellEpigenomicsAdapter, MethylationArrayEpigenomicsAdapter),
    ],
    ids=["bulk_vs_singlecell", "bulk_vs_array", "singlecell_vs_array"],
)
def test_qualifier_sets_are_pairwise_disjoint(left, right):
    overlap = left.qualifiers & right.qualifiers
    assert overlap == frozenset(), (
        f"Qualifiers overlap between {left.__name__} and "
        f"{right.__name__}: {overlap}"
    )


def test_union_of_qualifiers_has_expected_total_count():
    """Sanity check that the disjoint union has the expected total count.

    Running tally:
      v0.1.0 — Bulk (19) + Single-cell (8) + Array (8) = 35.
      v0.2.0 §2 — +chromatin_state_of (bulk) +marker_state_of_cell_type
                   (single-cell) = 37.
      v0.2.0 §5 — +classified_as_ccre (bulk) = 38.
      v0.2.0 §6 — 3D-genome qualifier split: bulk keeps
                   ``interacts_3d_with`` for one minor version (with
                   ``DeprecationWarning``) and adds six replacements
                   (contacts_with, loops_with, anchored_to,
                   tad_boundary_with, compartment_of, subtad_of), net
                   +6 = 44.

    If a qualifier is added or removed from any sibling, update this
    number after re-verifying that the pairwise disjointness invariant
    still holds.
    """
    union = (
        BulkEpigenomicsAdapter.qualifiers
        | SingleCellEpigenomicsAdapter.qualifiers
        | MethylationArrayEpigenomicsAdapter.qualifiers
    )
    assert len(union) == 44
    # The disjoint-union arithmetic must agree — this would only fail if one
    # of the pairwise tests above also failed.
    assert len(union) == (
        len(BulkEpigenomicsAdapter.qualifiers)
        + len(SingleCellEpigenomicsAdapter.qualifiers)
        + len(MethylationArrayEpigenomicsAdapter.qualifiers)
    )
