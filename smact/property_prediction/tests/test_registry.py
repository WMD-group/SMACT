"""Tests for the model registry module."""

from __future__ import annotations

import pytest

from smact.property_prediction.config import DEFAULT_MODELS, PROPERTY_METADATA
from smact.property_prediction.registry import (
    get_default_model,
    get_property_description,
    get_property_fidelities,
    get_property_unit,
    get_supported_properties,
    parse_model_name,
)


class TestGetSupportedProperties:
    """Tests for get_supported_properties function."""

    def test_returns_list(self):
        """Should return a list of property names."""
        properties = get_supported_properties()
        assert isinstance(properties, list)
        assert len(properties) > 0

    def test_contains_band_gap(self):
        """Should contain band_gap."""
        properties = get_supported_properties()
        assert "band_gap" in properties


class TestGetDefaultModel:
    """Tests for get_default_model function."""

    def test_band_gap_default(self):
        """Should return default model for band_gap."""
        model = get_default_model("band_gap")
        assert "band_gap" in model.lower()
        assert "Roost" in model

    def test_invalid_property_raises(self):
        """Should raise ValueError for unknown property."""
        with pytest.raises(ValueError, match="Unknown property"):
            get_default_model("nonexistent_property")


class TestGetPropertyFidelities:
    """Tests for get_property_fidelities function."""

    def test_band_gap_fidelities(self):
        """Band gap should return fidelities or None."""
        fidelities = get_property_fidelities("band_gap")
        # Currently no fidelity-specific models, so should be None
        assert fidelities is None or isinstance(fidelities, list)

    def test_invalid_property_raises(self):
        """Should raise ValueError for unknown property."""
        with pytest.raises(ValueError, match="Unknown property"):
            get_property_fidelities("nonexistent_property")


class TestGetPropertyUnit:
    """Tests for get_property_unit function."""

    def test_band_gap_unit(self):
        """Band gap should be in eV."""
        unit = get_property_unit("band_gap")
        assert unit == "eV"

    def test_unknown_property_returns_empty(self):
        """Unknown property should return empty string."""
        unit = get_property_unit("nonexistent_property")
        assert unit == ""


class TestGetPropertyDescription:
    """Tests for get_property_description function."""

    def test_band_gap_description(self):
        """Band gap should have a description."""
        desc = get_property_description("band_gap")
        assert len(desc) > 0
        assert "band" in desc.lower()

    def test_unknown_property_returns_empty(self):
        """Unknown property should return empty string."""
        desc = get_property_description("nonexistent_property")
        assert desc == ""


class TestParseModelName:
    """Tests for parse_model_name function."""

    def test_standard_model_name(self):
        """Should parse standard model name correctly."""
        result = parse_model_name("Roost-MP-2024.12.0-band_gap")
        assert result["model_type"] == "Roost"
        assert result["dataset"] == "MP"
        assert result["version"] == "2024.12.0"

    def test_short_model_name(self):
        """Should handle short/invalid model names gracefully."""
        result = parse_model_name("invalid")
        assert result["model_type"] is None


class TestConfigConsistency:
    """Tests for consistency between config dicts."""

    def test_default_models_have_metadata(self):
        """All properties in DEFAULT_MODELS should have metadata."""
        for prop in DEFAULT_MODELS:
            assert prop in PROPERTY_METADATA, f"Missing metadata for {prop}"

    def test_metadata_has_required_fields(self):
        """All metadata entries should have required fields."""
        required_fields = ["unit", "description"]
        for prop, meta in PROPERTY_METADATA.items():
            for field in required_fields:
                assert field in meta, f"Missing {field} for {prop}"
