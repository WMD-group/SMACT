"""Unit tests for smact.io.elementembeddings module."""

from __future__ import annotations

import unittest
from unittest.mock import MagicMock, patch

from smact.io.elementembeddings import (
    AllowedElementEmbeddings,
    AllowedSpeciesEmbeddings,
    PoolingStats,
    _check_elementembeddings,
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


class TestCheckElementEmbeddings(unittest.TestCase):
    """Test _check_elementembeddings helper."""

    def test_raises_import_error_when_not_installed(self):
        with patch.dict("sys.modules", {"elementembeddings": None}):
            with self.assertRaises(ImportError) as ctx:
                _check_elementembeddings()
            self.assertIn("ElementEmbeddings", str(ctx.exception))
            self.assertIn("pip install", str(ctx.exception))

    def test_passes_when_installed(self):
        mock_ee = MagicMock()
        with patch.dict("sys.modules", {"elementembeddings": mock_ee}):
            # Should not raise
            _check_elementembeddings()


class TestCompositionFeaturiser(unittest.TestCase):
    """Test composition_featuriser wrapper."""

    def _setup_mocks(self):
        """Create mock modules for elementembeddings."""
        mock_ee = MagicMock()
        mock_comp_mod = MagicMock()
        mock_featuriser = MagicMock(return_value="featurised_result")
        mock_comp_mod.composition_featuriser = mock_featuriser
        return mock_ee, mock_comp_mod, mock_featuriser

    def test_delegates_with_correct_args(self):
        mock_ee, mock_comp_mod, mock_featuriser = self._setup_mocks()
        with patch.dict(
            "sys.modules",
            {
                "elementembeddings": mock_ee,
                "elementembeddings.composition": mock_comp_mod,
            },
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
        mock_ee, mock_comp_mod, mock_featuriser = self._setup_mocks()
        with patch.dict(
            "sys.modules",
            {
                "elementembeddings": mock_ee,
                "elementembeddings.composition": mock_comp_mod,
            },
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
        with patch.dict("sys.modules", {"elementembeddings": None}), self.assertRaises(ImportError):
            composition_featuriser(composition_data=["NaCl"])


class TestSpeciesCompositionFeaturiser(unittest.TestCase):
    """Test species_composition_featuriser wrapper."""

    def _setup_mocks(self):
        """Create mock modules for elementembeddings."""
        mock_ee = MagicMock()
        mock_comp_mod = MagicMock()
        mock_featuriser = MagicMock(return_value="species_result")
        mock_comp_mod.species_composition_featuriser = mock_featuriser
        return mock_ee, mock_comp_mod, mock_featuriser

    def test_delegates_with_correct_args(self):
        mock_ee, mock_comp_mod, mock_featuriser = self._setup_mocks()
        with patch.dict(
            "sys.modules",
            {
                "elementembeddings": mock_ee,
                "elementembeddings.composition": mock_comp_mod,
            },
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
        mock_ee, mock_comp_mod, mock_featuriser = self._setup_mocks()
        with patch.dict(
            "sys.modules",
            {
                "elementembeddings": mock_ee,
                "elementembeddings.composition": mock_comp_mod,
            },
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
        with patch.dict("sys.modules", {"elementembeddings": None}), self.assertRaises(ImportError):
            species_composition_featuriser(
                composition_data=[{"Fe2+": 1, "O2-": 1}],
            )
