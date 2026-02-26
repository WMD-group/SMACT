"""Tests for the probability_models module (structure_prediction)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from smact.structure_prediction.probability_models import RadiusModel


class TestRadiusModel:
    """Test the RadiusModel substitution probability model."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.model = RadiusModel()

    def test_init_loads_shannon_data(self):
        """Shannon data is loaded as a DataFrame with expected columns."""
        assert isinstance(self.model.shannon_data, pd.DataFrame)
        assert "ionic_radius" in self.model.shannon_data.columns
        assert "charge" in self.model.shannon_data.columns

    def test_spring_constant_positive(self):
        """Spring constant k must be positive."""
        assert self.model.k > 0

    def test_sub_prob_same_species_is_one(self):
        """Substitution probability of a species with itself should be 1.0."""
        prob = self.model.sub_prob("Na1+", "Na1+")
        assert prob == pytest.approx(1.0)

    def test_sub_prob_between_zero_and_one(self):
        """Substitution probability should be in [0, 1]."""
        prob = self.model.sub_prob("Na1+", "K1+")
        assert 0.0 <= prob <= 1.0

    def test_sub_prob_symmetric(self):
        """Substitution probability should be symmetric: p(A,B) == p(B,A)."""
        p1 = self.model.sub_prob("Na1+", "K1+")
        p2 = self.model.sub_prob("K1+", "Na1+")
        assert p1 == pytest.approx(p2)

    def test_sub_prob_unknown_element_raises(self):
        """KeyError raised for species not in Shannon radius data."""
        with pytest.raises(KeyError, match="Element not in Shannon radius data"):
            self.model.sub_prob("Xx9+", "Na1+")

    def test_gen_lambda_shape(self):
        """gen_lambda returns a square DataFrame for given species."""
        species = ["Na1+", "K1+", "Li1+"]
        table = self.model.gen_lambda(species)
        assert isinstance(table, pd.DataFrame)
        assert table.shape == (3, 3)

    def test_gen_lambda_symmetric(self):
        """gen_lambda table should be symmetric."""
        species = ["Na1+", "K1+"]
        table = self.model.gen_lambda(species)
        assert table.loc["Na1+", "K1+"] == pytest.approx(table.loc["K1+", "Na1+"])

    def test_gen_lambda_diagonal_is_one(self):
        """Diagonal entries (self-substitution) should be 1.0."""
        species = ["Na1+", "K1+"]
        table = self.model.gen_lambda(species)
        for sp in species:
            assert table.loc[sp, sp] == pytest.approx(1.0)

    def test_gen_lambda_single_species(self):
        """gen_lambda with one species returns 1x1 table with value 1.0."""
        table = self.model.gen_lambda(["Na1+"])
        assert table.shape == (1, 1)
        assert table.iloc[0, 0] == pytest.approx(1.0)

    def test_sub_prob_different_radii_lower(self):
        """Species with different radii should have lower probability than self-substitution."""
        prob_diff = self.model.sub_prob("Ba2+", "Be2+")
        prob_same = self.model.sub_prob("Ba2+", "Ba2+")
        assert prob_diff < prob_same
