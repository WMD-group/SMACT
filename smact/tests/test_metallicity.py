"""Tests for the metallicity module."""

from __future__ import annotations

import math
import unittest

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


class TestMetallicity(unittest.TestCase):
    """Test the metallicity module functionality."""

    def setUp(self):
        """Set up test cases."""
        # Known metallic/alloy compounds
        self.metallic_compounds = [
            "Fe3Al",  # Classic intermetallic
            "Ni3Ti",  # Superalloy component
            "Cu3Au",  # Ordered alloy
            "Fe2Nb",  # Laves phase
        ]

        # Known non-metallic compounds
        self.non_metallic_compounds = [
            "NaCl",  # Ionic
            "SiO2",  # Covalent
            "Fe2O3",  # Metal oxide
            "CuSO4",  # Complex ionic
        ]

    def test_get_element_fraction(self):
        """Test the helper function for element fraction calculations."""
        # Test with metals set
        self.assertAlmostEqual(
            get_element_fraction(Composition("Fe3Al"), smact.metals),
            1.0,
            places=6,
            msg="Expected all elements in Fe3Al to be metals",
        )

        # Test with d-block set
        self.assertAlmostEqual(
            get_element_fraction(Composition("Fe2Nb"), smact.d_block),
            1.0,
            places=6,
            msg="Expected all elements in Fe2Nb to be d-block",
        )

        # Test with empty set
        self.assertAlmostEqual(
            get_element_fraction(Composition("Fe3Al"), set()),
            0.0,
            places=6,
            msg="Expected zero fraction for empty element set",
        )

    def test_get_metal_fraction(self):
        """Test metal fraction calculation."""
        # Should be 1.0 for pure metallic compounds
        self.assertAlmostEqual(
            get_metal_fraction(Composition("Fe3Al")),
            1.0,
            places=6,
            msg="Expected pure metallic composition for Fe3Al",
        )

        # Should be 0.0 for compounds with no metals
        self.assertAlmostEqual(
            get_metal_fraction(Composition("SiO2")),
            0.0,
            places=6,
            msg="Expected no metallic elements in SiO2",
        )

        # Should be fractional for mixed compounds
        fe2o3 = get_metal_fraction(Composition("Fe2O3"))
        self.assertTrue(
            0 < fe2o3 < 1,
            msg=f"Expected fractional metal content for Fe2O3, got {fe2o3:.2f}",
        )

    def test_get_d_block_element_fraction(self):
        """Test d-block element fraction calculation."""
        # Should be 1.0 for pure transition metal compounds
        self.assertAlmostEqual(
            get_d_block_element_fraction(Composition("Fe2Nb")),
            1.0,
            places=6,
            msg="Expected all d-block elements in Fe2Nb",
        )

        # Should be 0.0 for compounds with no d-block elements
        self.assertAlmostEqual(
            get_d_block_element_fraction(Composition("NaCl")),
            0.0,
            places=6,
            msg="Expected no d-block elements in NaCl",
        )

    def test_get_distinct_metal_count(self):
        """Test counting of distinct metals."""
        self.assertEqual(
            get_distinct_metal_count(Composition("Fe3Al")),
            2,
            msg="Expected 2 distinct metals in Fe3Al",
        )
        self.assertEqual(
            get_distinct_metal_count(Composition("NaCl")),
            1,
            msg="Expected 1 distinct metal in NaCl",
        )
        self.assertEqual(
            get_distinct_metal_count(Composition("SiO2")),
            0,
            msg="Expected no metals in SiO2",
        )

    def test_get_pauling_test_mismatch(self):
        """Test Pauling electronegativity mismatch calculation."""
        # Ionic compounds should have high mismatch
        nacl_mismatch = get_pauling_test_mismatch(Composition("NaCl"))

        # Metallic compounds should have lower mismatch
        fe3al_mismatch = get_pauling_test_mismatch(Composition("Fe3Al"))

        self.assertTrue(
            fe3al_mismatch < nacl_mismatch,
            msg=f"Expected lower Pauling mismatch for Fe3Al ({fe3al_mismatch:.2f}) than NaCl ({nacl_mismatch:.2f})",
        )

    def test_metallicity_score(self):
        """Test the overall metallicity scoring function."""
        # Known metallic compounds should score high
        for formula in self.metallic_compounds:
            score = metallicity_score(formula)
            self.assertTrue(
                score > 0.7,
                msg=f"Expected high metallicity score (>0.7) for {formula}, but got {score:.2f}",
            )

        # Non-metallic compounds should score low
        for formula in self.non_metallic_compounds:
            score = metallicity_score(formula)
            self.assertTrue(
                score < 0.5,
                msg=f"Expected low metallicity score (<0.5) for {formula}, but got {score:.2f}",
            )

    def test_pauling_mismatch_none_electronegativity(self):
        """get_pauling_test_mismatch returns NaN when element lacks Pauling eneg (line 108)."""
        # He has no Pauling electronegativity
        result = get_pauling_test_mismatch(Composition("HeNe"))
        self.assertTrue(math.isnan(result))

    def test_metallicity_score_vec_value_error(self):
        """metallicity_score handles ValueError from VEC (line 143-144)."""
        # Lr has no valence data, triggering ValueError in VEC → vec_factor=0.5 fallback
        score = metallicity_score("LrFe")
        self.assertTrue(0.0 <= score <= 1.0)

    def test_metallicity_score_pauling_nan_path(self):
        """metallicity_score handles NaN Pauling mismatch (line 149)."""
        # He has no Pauling eneg → pauling_mismatch is NaN → pauling_term=0.5
        score = metallicity_score("HeFe")
        self.assertTrue(0.0 <= score <= 1.0)

    def test_edge_cases(self):
        """Test edge cases and error handling."""
        # Single element
        score = metallicity_score("Fe")
        self.assertTrue(0.0 <= score <= 1.0, msg=f"Expected score between 0 and 1 for Fe, got {score:.2f}")

        # Empty composition -> expect ValueError("Empty composition")
        with pytest.raises(ValueError, match="Empty composition"):
            metallicity_score("")

        # Invalid formula -> e.g. "NotAnElement"
        with pytest.raises(ValueError, match="Invalid formula"):
            metallicity_score("NotAnElement")
