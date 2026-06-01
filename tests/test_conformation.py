"""3D-genome split + conformation-index behaviour tests.

Covers:
- ``CONFORMATION_INDEX.lookup()`` for the 7 compartment / sub-compartment
  labels across notational variants (``A``, ``Compartment A``,
  ``A_compartment``).
- End-to-end ``compartment_of`` resolution: bare ``A`` / ``B1`` etc.
  routes through CONFORMATION to ``EpiConf:``.
- Six new 3D-genome qualifiers are present on bulk and disjoint from
  the sibling adapters.
- ``interacts_3d_with`` still resolves under the v0.1.0 behaviour
  (treated as a generic contact via the ``_DEPRECATED_QUALIFIERS`` map
  in ``bulk.py``) AND emits a ``DeprecationWarning`` so callers get a
  visible migration nudge.
- The other 5 new 3D qualifiers (``contacts_with``, ``loops_with``,
  ``anchored_to``, ``tad_boundary_with``, ``subtad_of``) accept
  arbitrary subject / object strings without crashing, even when the
  resolver can't CURIE them — they're descriptive qualifiers, not
  controlled-vocabulary ones.
"""
from __future__ import annotations

import re
import warnings

import pytest

from indicium_adapters_epigenomics import (
    BulkEpigenomicsAdapter,
    MethylationArrayEpigenomicsAdapter,
    SingleCellEpigenomicsAdapter,
)
from indicium_adapters_epigenomics._conformation import (
    CONFORMATION_INDEX,
    ConformationIndex,
    _load_conformation_index_rows,
)

# --- Index integrity --------------------------------------------------------


def test_conformation_index_loads_seven_compartment_labels():
    assert isinstance(CONFORMATION_INDEX, ConformationIndex)
    # A / B + A1 / A2 / B1 / B2 / B3 — Rao et al. (2014) refinement.
    assert len(CONFORMATION_INDEX) == 7


def test_every_row_has_epiconf_curie_with_clean_localpart():
    snake_compat = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
    rows = _load_conformation_index_rows()
    bad: list[tuple[str, str]] = []
    for row in rows:
        curie = row["curie"]
        if not curie.startswith("EpiConf:"):
            bad.append((row["symbol"], curie))
            continue
        localpart = curie.split(":", 1)[1]
        if not snake_compat.match(localpart):
            bad.append((row["symbol"], curie))
    assert not bad, f"malformed EpiConf CURIE local parts: {bad}"


def test_every_row_type_is_snake_case():
    snake = re.compile(r"^[a-z][a-z0-9_]*$")
    rows = _load_conformation_index_rows()
    offenders = [row["symbol"] for row in rows if not snake.match(row["type"])]
    assert not offenders, f"non-snake_case conformation types: {offenders}"


# --- Compartment-label lookup variants -------------------------------------


@pytest.mark.parametrize(
    "text, expected_curie",
    [
        # Canonical labels.
        ("A", "EpiConf:A"),
        ("B", "EpiConf:B"),
        ("A1", "EpiConf:A1"),
        ("A2", "EpiConf:A2"),
        ("B1", "EpiConf:B1"),
        ("B2", "EpiConf:B2"),
        ("B3", "EpiConf:B3"),
        # Case-insensitive.
        ("a", "EpiConf:A"),
        ("b3", "EpiConf:B3"),
        # Long-form aliases.
        ("Compartment A", "EpiConf:A"),
        ("A_compartment", "EpiConf:A"),
        ("compartment_A", "EpiConf:A"),
        ("A1 subcompartment", "EpiConf:A1"),
        ("subA1", "EpiConf:A1"),
        ("A1_sub", "EpiConf:A1"),
    ],
)
def test_conformation_lookup_resolves_known_variants(text, expected_curie):
    match = CONFORMATION_INDEX.lookup(text)
    assert match is not None, f"{text!r} should resolve"
    assert match["curie"] == expected_curie


def test_conformation_lookup_returns_none_for_unknown():
    assert CONFORMATION_INDEX.lookup("C") is None  # no C compartment
    assert CONFORMATION_INDEX.lookup("D7") is None


def test_conformation_lookup_returns_none_for_non_string():
    assert CONFORMATION_INDEX.lookup(None) is None  # type: ignore[arg-type]


# --- End-to-end compartment_of routing -------------------------------------


@pytest.mark.parametrize(
    "compartment, expected_curie",
    [
        ("A", "EpiConf:A"),
        ("B", "EpiConf:B"),
        ("A1", "EpiConf:A1"),
        ("Compartment A", "EpiConf:A"),
    ],
)
def test_compartment_of_resolves_to_epiconf(compartment, expected_curie):
    bulk = BulkEpigenomicsAdapter()
    enriched = bulk.enrich_claim(
        {
            "id": "c-conf-1",
            "subject": "region_X",
            "qualifier": "compartment_of",
            "object": compartment,
        }
    )
    assert enriched["ontology_terms"]["object"] == expected_curie
    assert enriched["xrefs"]["object"]["source"] == "conformation"


# --- Six new qualifiers exist on bulk, disjoint elsewhere -------------------


@pytest.mark.parametrize(
    "qualifier",
    [
        "contacts_with",
        "loops_with",
        "anchored_to",
        "tad_boundary_with",
        "compartment_of",
        "subtad_of",
    ],
)
def test_new_3d_qualifier_lives_on_bulk_only(qualifier):
    assert qualifier in BulkEpigenomicsAdapter.qualifiers
    assert qualifier not in SingleCellEpigenomicsAdapter.qualifiers
    assert qualifier not in MethylationArrayEpigenomicsAdapter.qualifiers


@pytest.mark.parametrize(
    "qualifier",
    [
        "contacts_with",
        "loops_with",
        "anchored_to",
        "tad_boundary_with",
        "subtad_of",
    ],
)
def test_descriptive_3d_qualifier_accepts_arbitrary_string_object(qualifier):
    """The 5 non-compartment qualifiers take coordinate / anchor / TAD-ID
    strings that don't resolve to any controlled vocabulary. Verify
    ``enrich_claim`` accepts them without error and leaves the object
    un-CURIE'd.
    """
    bulk = BulkEpigenomicsAdapter()
    enriched = bulk.enrich_claim(
        {
            "id": f"c-{qualifier}",
            "subject": "chr1:100000-200000",
            "qualifier": qualifier,
            "object": "chr1:300000-400000",
        }
    )
    assert isinstance(enriched, dict)
    assert "object" not in enriched["ontology_terms"]


# --- Deprecation alias ------------------------------------------------------


def test_interacts_3d_with_still_resolves_under_contact_semantics():
    """Backwards-compat: a v0.1.0-style claim with the deprecated
    qualifier still enriches as if ``contacts_with`` had been used.
    """
    bulk = BulkEpigenomicsAdapter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        enriched = bulk.enrich_claim(
            {
                "id": "c-legacy",
                "subject": "region_X",
                "qualifier": "interacts_3d_with",
                "object": "region_Y",
            }
        )
    # Legacy qualifier name preserved verbatim in the claim payload.
    assert enriched["qualifier"] == "interacts_3d_with"
    # Object doesn't resolve (no controlled vocab for raw coordinates),
    # but the claim is well-formed and enrich_claim didn't raise.
    assert isinstance(enriched, dict)


def test_interacts_3d_with_emits_deprecation_warning():
    """Visibility: callers must get a runtime nudge when they use the
    deprecated qualifier. Use ``pytest.warns`` rather than
    ``catch_warnings`` so the absence of the warning fails the test.
    """
    bulk = BulkEpigenomicsAdapter()
    with pytest.warns(DeprecationWarning, match="interacts_3d_with"):
        bulk.enrich_claim(
            {
                "id": "c-warn",
                "subject": "region_X",
                "qualifier": "interacts_3d_with",
                "object": "region_Y",
            }
        )


def test_new_3d_qualifiers_do_not_emit_deprecation_warning():
    """Negative control: the replacement qualifiers must NOT emit a
    DeprecationWarning — otherwise users migrating from the legacy
    qualifier get the same warning forever.
    """
    bulk = BulkEpigenomicsAdapter()
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        # Each replacement qualifier should be silent.
        for qualifier in (
            "contacts_with",
            "loops_with",
            "anchored_to",
            "tad_boundary_with",
            "compartment_of",
            "subtad_of",
        ):
            bulk.enrich_claim(
                {
                    "id": f"c-{qualifier}",
                    "subject": "region_X",
                    "qualifier": qualifier,
                    "object": "A" if qualifier == "compartment_of" else "region_Y",
                }
            )
