"""Tests for the intermetallics module."""

from __future__ import annotations

import unittest

import pytest
from pymatgen.core import Composition

from smact.intermetallics import (
    get_d_electron_fraction,
    get_distinct_metal_count,
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

    def test_get_metal_fraction(self):
        """Test metal fraction calculation."""
        # Should be 1.0 for pure intermetallics
        self.assertEqual(get_metal_fraction(Composition("Fe3Al")), 1.0)

        # Should be 0.0 for compounds with no metals
        self.assertEqual(get_metal_fraction(Composition("SiO2")), 0.0)

        # Should be fractional for mixed compounds
        fe2o3 = get_metal_fraction(Composition("Fe2O3"))
        self.assertTrue(0 < fe2o3 < 1)

    def test_get_d_electron_fraction(self):
        """Test d-electron fraction calculation."""
        # Should be 1.0 for pure transition metal compounds
        self.assertEqual(get_d_electron_fraction(Composition("Fe2Nb")), 1.0)

        # Should be 0.0 for compounds with no d-block elements
        self.assertEqual(get_d_electron_fraction(Composition("NaCl")), 0.0)

    def test_get_distinct_metal_count(self):
        """Test counting of distinct metals."""
        self.assertEqual(get_distinct_metal_count(Composition("Fe3Al")), 2)
        self.assertEqual(get_distinct_metal_count(Composition("NaCl")), 1)
        self.assertEqual(get_distinct_metal_count(Composition("SiO2")), 0)

    def test_get_pauling_test_mismatch(self):
        """Test Pauling electronegativity mismatch calculation."""
        # Ionic compounds should have high mismatch
        nacl_mismatch = get_pauling_test_mismatch(Composition("NaCl"))

        # Intermetallics should have lower mismatch
        fe3al_mismatch = get_pauling_test_mismatch(Composition("Fe3Al"))

        self.assertTrue(fe3al_mismatch < nacl_mismatch)

    def test_intermetallic_score(self):
        """Test the overall intermetallic scoring function."""
        # Known intermetallics should score high
        for formula in self.intermetallics:
            score = intermetallic_score(formula)
            self.assertTrue(score > 0.7, f"Expected high score (>0.7) for {formula}, got {score}")

        # Non-intermetallics should score low
        for formula in self.non_intermetallics:
            score = intermetallic_score(formula)
            self.assertTrue(score < 0.5, f"Expected low score (<0.5) for {formula}, got {score}")

    def test_edge_cases(self):
        """Test edge cases and error handling."""
        # Single element
        self.assertTrue(0.0 <= intermetallic_score("Fe") <= 1.0)

        # Empty composition should not crash
        with pytest.raises(ValueError):
            intermetallic_score("")

        # Invalid formula should raise error
        with pytest.raises(ValueError):
            intermetallic_score("NotAnElement")


if __name__ == "__main__":
    unittest.main()
