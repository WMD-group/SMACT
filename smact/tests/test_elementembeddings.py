"""Unit tests for smact.io.elementembeddings module."""

from __future__ import annotations

import unittest
from unittest.mock import MagicMock, patch

import pytest

from smact.io.elementembeddings import (
    AllowedElementEmbeddings,
    AllowedSpeciesEmbeddings,
    PoolingStats,
    composition_featuriser,
    species_composition_featuriser,
)


class TestEnums(unittest.TestCase):
    """Test the enum classes defined in the module."""

    def test_allowed_element_embeddings_members(self):
        expected = {
            "magpie": "magpie",
            "mat2vec": "mat2vec",
            "skipatom": "skipatom",
            "cgnf": "cgnf",
            "xenonpy": "xenonpy",
            "random": "random",
            "oliynyk": "oliynyk",
            "matscholar": "matscholar",
            "crystallm": "crystallm",
            "megnet": "megnet",
        }
        self.assertEqual(len(AllowedElementEmbeddings), 10)
        for name, value in expected.items():
            self.assertEqual(AllowedElementEmbeddings[name].value, value)

    def test_allowed_species_embeddings_members(self):
        self.assertEqual(len(AllowedSpeciesEmbeddings), 1)
        self.assertEqual(AllowedSpeciesEmbeddings.skipspecies.value, "skipspecies")

    def test_pooling_stats_members(self):
        expected = {
            "mean": "mean",
            "variance": "variance",
            "minpool": "minpool",
            "maxpool": "maxpool",
            "range": "range",
            "sum": "sum",
            "geometric_mean": "geometric_mean",
            "harmonic_mean": "harmonic_mean",
        }
        self.assertEqual(len(PoolingStats), 8)
        for name, value in expected.items():
            self.assertEqual(PoolingStats[name].value, value)


class TestHasElementEmbeddingsFlag(unittest.TestCase):
    """Test the HAS_ELEMENTEMBEDDINGS flag behaviour."""

    def test_raises_import_error_when_flag_false(self):
        with (
            patch("smact.io.elementembeddings.HAS_ELEMENTEMBEDDINGS", new=False),
            pytest.raises(ImportError, match="ElementEmbeddings"),
        ):
            composition_featuriser(composition_data=["NaCl"])

    def test_flag_true_does_not_raise(self):
        mock_featuriser = MagicMock(return_value="result")
        with (
            patch("smact.io.elementembeddings.HAS_ELEMENTEMBEDDINGS", new=True),
            patch("smact.io.elementembeddings.ee_composition_featuriser", mock_featuriser),
        ):
            composition_featuriser(composition_data=["NaCl"])


class TestCompositionFeaturiser(unittest.TestCase):
    """Test composition_featuriser wrapper."""

    def test_delegates_with_correct_args(self):
        mock_featuriser = MagicMock(return_value="featurised_result")
        with (
            patch("smact.io.elementembeddings.HAS_ELEMENTEMBEDDINGS", new=True),
            patch("smact.io.elementembeddings.ee_composition_featuriser", mock_featuriser),
        ):
            result = composition_featuriser(
                composition_data=["NaCl"],
                formula_column="formula",
                embedding=AllowedElementEmbeddings.mat2vec,
                stats=PoolingStats.variance,
                inplace=True,
            )
        mock_featuriser.assert_called_once_with(
            data=["NaCl"],
            formula_column="formula",
            embedding=AllowedElementEmbeddings.mat2vec,
            stats=PoolingStats.variance,
            inplace=True,
        )
        self.assertEqual(result, "featurised_result")

    def test_passes_default_args(self):
        mock_featuriser = MagicMock(return_value="featurised_result")
        with (
            patch("smact.io.elementembeddings.HAS_ELEMENTEMBEDDINGS", new=True),
            patch("smact.io.elementembeddings.ee_composition_featuriser", mock_featuriser),
        ):
            composition_featuriser(composition_data=["NaCl"])
        mock_featuriser.assert_called_once_with(
            data=["NaCl"],
            formula_column="formula",
            embedding=AllowedElementEmbeddings.magpie,
            stats=PoolingStats.mean,
            inplace=False,
        )

    def test_raises_import_error_when_not_installed(self):
        with (
            patch("smact.io.elementembeddings.HAS_ELEMENTEMBEDDINGS", new=False),
            pytest.raises(ImportError, match="ElementEmbeddings"),
        ):
            composition_featuriser(composition_data=["NaCl"])


class TestSpeciesCompositionFeaturiser(unittest.TestCase):
    """Test species_composition_featuriser wrapper."""

    def test_delegates_with_correct_args(self):
        mock_featuriser = MagicMock(return_value="species_result")
        with (
            patch("smact.io.elementembeddings.HAS_ELEMENTEMBEDDINGS", new=True),
            patch("smact.io.elementembeddings.ee_species_composition_featuriser", mock_featuriser),
        ):
            result = species_composition_featuriser(
                composition_data=[{"Na+": 1, "Cl-": 1}],
                embedding="skipspecies",
                stats=[PoolingStats.mean, PoolingStats.variance],
                to_dataframe=True,
            )
        mock_featuriser.assert_called_once_with(
            data=[{"Na+": 1, "Cl-": 1}],
            embedding="skipspecies",
            stats=[PoolingStats.mean, PoolingStats.variance],
            to_dataframe=True,
        )
        self.assertEqual(result, "species_result")

    def test_passes_default_args(self):
        mock_featuriser = MagicMock(return_value="species_result")
        with (
            patch("smact.io.elementembeddings.HAS_ELEMENTEMBEDDINGS", new=True),
            patch("smact.io.elementembeddings.ee_species_composition_featuriser", mock_featuriser),
        ):
            species_composition_featuriser(
                composition_data=[{"Fe2+": 1, "O2-": 1}],
            )
        mock_featuriser.assert_called_once_with(
            data=[{"Fe2+": 1, "O2-": 1}],
            embedding=AllowedSpeciesEmbeddings.skipspecies,
            stats=PoolingStats.mean,
            to_dataframe=False,
        )

    def test_raises_import_error_when_not_installed(self):
        with (
            patch("smact.io.elementembeddings.HAS_ELEMENTEMBEDDINGS", new=False),
            pytest.raises(ImportError, match="ElementEmbeddings"),
        ):
            species_composition_featuriser(
                composition_data=[{"Fe2+": 1, "O2-": 1}],
            )
