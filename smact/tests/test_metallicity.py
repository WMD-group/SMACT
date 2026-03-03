"""Tests for the metallicity module."""

from __future__ import annotations

import math

import pytest
from pymatgen.core import Composition

import smact
from smact.metallicity import (
    get_d_block_element_fraction,
    get_distinct_metal_count,
    get_element_fraction,
    get_metal_fraction,
    get_pauling_test_mismatch,
    metallicity_score,
)


@pytest.fixture
def metallic_compounds():
    return [
        "Fe3Al",  # Classic intermetallic
        "Ni3Ti",  # Superalloy component
        "Cu3Au",  # Ordered alloy
        "Fe2Nb",  # Laves phase
    ]


@pytest.fixture
def non_metallic_compounds():
    return [
        "NaCl",  # Ionic
        "SiO2",  # Covalent
        "Fe2O3",  # Metal oxide
        "CuSO4",  # Complex ionic
    ]


def test_get_element_fraction():
    # Test with metals set
    assert get_element_fraction(Composition("Fe3Al"), smact.metals) == pytest.approx(1.0, abs=1e-6)

    # Test with d-block set
    assert get_element_fraction(Composition("Fe2Nb"), smact.d_block) == pytest.approx(1.0, abs=1e-6)

    # Test with empty set
    assert get_element_fraction(Composition("Fe3Al"), set()) == pytest.approx(0.0, abs=1e-6)


def test_get_metal_fraction():
    # Should be 1.0 for pure metallic compounds
    assert get_metal_fraction(Composition("Fe3Al")) == pytest.approx(1.0, abs=1e-6)

    # Should be 0.0 for compounds with no metals
    assert get_metal_fraction(Composition("SiO2")) == pytest.approx(0.0, abs=1e-6)

    # Should be fractional for mixed compounds
    fe2o3 = get_metal_fraction(Composition("Fe2O3"))
    assert 0 < fe2o3 < 1


def test_get_d_block_element_fraction():
    # Should be 1.0 for pure transition metal compounds
    assert get_d_block_element_fraction(Composition("Fe2Nb")) == pytest.approx(1.0, abs=1e-6)

    # Should be 0.0 for compounds with no d-block elements
    assert get_d_block_element_fraction(Composition("NaCl")) == pytest.approx(0.0, abs=1e-6)


def test_get_distinct_metal_count():
    assert get_distinct_metal_count(Composition("Fe3Al")) == 2
    assert get_distinct_metal_count(Composition("NaCl")) == 1
    assert get_distinct_metal_count(Composition("SiO2")) == 0


def test_get_pauling_test_mismatch():
    # Ionic compounds should have high mismatch
    nacl_mismatch = get_pauling_test_mismatch(Composition("NaCl"))

    # Metallic compounds should have lower mismatch
    fe3al_mismatch = get_pauling_test_mismatch(Composition("Fe3Al"))

    assert fe3al_mismatch < nacl_mismatch


@pytest.mark.parametrize(
    "formula",
    [
        "Fe3Al",  # Classic intermetallic
        "Ni3Ti",  # Superalloy component
        "Cu3Au",  # Ordered alloy
        "Fe2Nb",  # Laves phase
    ],
)
def test_metallicity_score_metallic(formula):
    score = metallicity_score(formula)
    assert score > 0.7


@pytest.mark.parametrize(
    "formula",
    [
        "NaCl",  # Ionic
        "SiO2",  # Covalent
        "Fe2O3",  # Metal oxide
        "CuSO4",  # Complex ionic
    ],
)
def test_metallicity_score_non_metallic(formula):
    score = metallicity_score(formula)
    assert score < 0.5


def test_pauling_mismatch_none_electronegativity():
    # He has no Pauling electronegativity
    result = get_pauling_test_mismatch(Composition("HeNe"))
    assert math.isnan(result)


def test_metallicity_score_vec_value_error():
    # Lr has no valence data, triggering ValueError in VEC -> vec_factor=0.5 fallback
    score = metallicity_score("LrFe")
    assert 0.0 <= score <= 1.0


def test_metallicity_score_pauling_nan_path():
    # He has no Pauling eneg -> pauling_mismatch is NaN -> pauling_term=0.5
    score = metallicity_score("HeFe")
    assert 0.0 <= score <= 1.0


def test_edge_cases():
    # Single element
    score = metallicity_score("Fe")
    assert 0.0 <= score <= 1.0

    # Empty composition -> expect ValueError("Empty composition")
    with pytest.raises(ValueError, match="Empty composition"):
        metallicity_score("")

    # Invalid formula -> e.g. "NotAnElement"
    with pytest.raises(ValueError, match="Invalid formula"):
        metallicity_score("NotAnElement")
