"""Tests for the intermetallics module."""

from __future__ import annotations

import unittest

import pytest
from pymatgen.core import Composition

import smact
from smact.intermetallics import (
    get_d_electron_fraction,
    get_distinct_metal_count,
    get_element_fraction,
    get_metal_fraction,
    get_pauling_test_mismatch,
    intermetallic_score,
)


class TestIntermetallics(unittest.TestCase):
    """Test the intermetallics module functionality."""

    def setUp(self):
        """Set up test cases."""
        # Known intermetallics
        self.intermetallics = [
            "Fe3Al",  # Classic intermetallic
            "Ni3Ti",  # Superalloy component
            "Cu3Au",  # Ordered alloy
            "Fe2Nb",  # Laves phase
        ]

        # Known non-intermetallics
        self.non_intermetallics = [
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
        # Should be 1.0 for pure intermetallics
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

    def test_get_d_electron_fraction(self):
        """Test d-electron fraction calculation."""
        # Should be 1.0 for pure transition metal compounds
        self.assertAlmostEqual(
            get_d_electron_fraction(Composition("Fe2Nb")),
            1.0,
            places=6,
            msg="Expected all d-block elements in Fe2Nb",
        )

        # Should be 0.0 for compounds with no d-block elements
        self.assertAlmostEqual(
            get_d_electron_fraction(Composition("NaCl")),
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

        # Intermetallics should have lower mismatch
        fe3al_mismatch = get_pauling_test_mismatch(Composition("Fe3Al"))

        self.assertTrue(
            fe3al_mismatch < nacl_mismatch,
            msg=f"Expected lower Pauling mismatch for Fe3Al ({fe3al_mismatch:.2f}) than NaCl ({nacl_mismatch:.2f})",
        )

    def test_intermetallic_score(self):
        """Test the overall intermetallic scoring function."""
        # Known intermetallics should score high
        for formula in self.intermetallics:
            score = intermetallic_score(formula)
            self.assertTrue(
                score > 0.7,
                msg=f"Expected high intermetallic score (>0.7) for {formula}, but got {score:.2f}",
            )

        # Non-intermetallics should score low
        for formula in self.non_intermetallics:
            score = intermetallic_score(formula)
            self.assertTrue(
                score < 0.5,
                msg=f"Expected low intermetallic score (<0.5) for {formula}, but got {score:.2f}",
            )

    def test_edge_cases(self):
        """Test edge cases and error handling."""
        # Single element
        score = intermetallic_score("Fe")
        self.assertTrue(
            0.0 <= score <= 1.0,
            msg=f"Expected score between 0 and 1 for Fe, got {score:.2f}",
        )

        # Empty composition should not crash
        with pytest.raises(ValueError, match="Empty composition"):
            intermetallic_score("")

        # Invalid formula should raise error
        with pytest.raises(ValueError, match="Invalid formula"):
            intermetallic_score("NotAnElement")
