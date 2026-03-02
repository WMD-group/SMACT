"""Unit tests for smact.io.elementembeddings module."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from smact.io.elementembeddings import (
    AllowedElementEmbeddings,
    AllowedSpeciesEmbeddings,
    PoolingStats,
    composition_featuriser,
    species_composition_featuriser,
)


def test_allowed_element_embeddings_members():
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
    assert len(AllowedElementEmbeddings) == 10
    for name, value in expected.items():
        assert AllowedElementEmbeddings[name].value == value


def test_allowed_species_embeddings_members():
    assert len(AllowedSpeciesEmbeddings) == 1
    assert AllowedSpeciesEmbeddings.skipspecies.value == "skipspecies"


def test_pooling_stats_members():
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
    assert len(PoolingStats) == 8
    for name, value in expected.items():
        assert PoolingStats[name].value == value


def test_raises_import_error_when_not_available():
    with (
        patch("smact.io.elementembeddings._ee_composition_featuriser", new=None),
        pytest.raises(ImportError, match="ElementEmbeddings"),
    ):
        composition_featuriser(composition_data=["NaCl"])


def test_does_not_raise_when_available():
    mock_featuriser = MagicMock(return_value="result")
    with patch("smact.io.elementembeddings._ee_composition_featuriser", mock_featuriser):
        composition_featuriser(composition_data=["NaCl"])


def test_composition_featuriser_delegates_with_correct_args():
    mock_featuriser = MagicMock(return_value="featurised_result")
    with patch("smact.io.elementembeddings._ee_composition_featuriser", mock_featuriser):
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
    assert result == "featurised_result"


def test_composition_featuriser_passes_default_args():
    mock_featuriser = MagicMock(return_value="featurised_result")
    with patch("smact.io.elementembeddings._ee_composition_featuriser", mock_featuriser):
        composition_featuriser(composition_data=["NaCl"])
    mock_featuriser.assert_called_once_with(
        data=["NaCl"],
        formula_column="formula",
        embedding=AllowedElementEmbeddings.magpie,
        stats=PoolingStats.mean,
        inplace=False,
    )


def test_composition_featuriser_raises_import_error_when_not_installed():
    with (
        patch("smact.io.elementembeddings._ee_composition_featuriser", new=None),
        pytest.raises(ImportError, match="ElementEmbeddings"),
    ):
        composition_featuriser(composition_data=["NaCl"])


def test_species_composition_featuriser_delegates_with_correct_args():
    mock_featuriser = MagicMock(return_value="species_result")
    with patch("smact.io.elementembeddings._ee_species_composition_featuriser", mock_featuriser):
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
    assert result == "species_result"


def test_species_composition_featuriser_passes_default_args():
    mock_featuriser = MagicMock(return_value="species_result")
    with patch("smact.io.elementembeddings._ee_species_composition_featuriser", mock_featuriser):
        species_composition_featuriser(
            composition_data=[{"Fe2+": 1, "O2-": 1}],
        )
    mock_featuriser.assert_called_once_with(
        data=[{"Fe2+": 1, "O2-": 1}],
        embedding=AllowedSpeciesEmbeddings.skipspecies,
        stats=PoolingStats.mean,
        to_dataframe=False,
    )


def test_species_composition_featuriser_raises_import_error_when_not_installed():
    with (
        patch("smact.io.elementembeddings._ee_species_composition_featuriser", new=None),
        pytest.raises(ImportError, match="ElementEmbeddings"),
    ):
        species_composition_featuriser(
            composition_data=[{"Fe2+": 1, "O2-": 1}],
        )
