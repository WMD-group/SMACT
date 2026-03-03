"""Tests for the base predictor module."""

from __future__ import annotations

import numpy as np
import pytest

from smact.property_prediction.base_predictor import (
    BasePropertyPredictor,
    PredictionResult,
)


class TestPredictionResult:
    """Tests for the PredictionResult dataclass."""

    def test_basic_creation(self):
        """Should create a basic prediction result."""
        preds = np.array([1.0, 2.0, 3.0])
        result = PredictionResult(predictions=preds)

        assert len(result) == 3
        assert result[0] == 1.0
        assert result.uncertainties is None
        assert result.unit == ""

    def test_with_uncertainties(self):
        """Should handle uncertainty estimates."""
        preds = np.array([1.0, 2.0])
        uncs = np.array([0.1, 0.2])
        result = PredictionResult(predictions=preds, uncertainties=uncs, unit="eV")

        assert result.uncertainties is not None
        assert len(result.uncertainties) == 2
        assert result.unit == "eV"

    def test_with_compositions(self):
        """Should store composition strings."""
        preds = np.array([1.0, 2.0])
        comps = ["NaCl", "TiO2"]
        result = PredictionResult(predictions=preds, compositions=comps)

        assert result.compositions == comps

    def test_repr(self):
        """Should have informative repr."""
        preds = np.array([1.0, 2.0])
        result = PredictionResult(predictions=preds, unit="eV")

        repr_str = repr(result)
        assert "n=2" in repr_str
        assert "has_uncertainty=False" in repr_str

    def test_to_dict(self):
        """Should convert to dictionary."""
        preds = np.array([1.0, 2.0])
        uncs = np.array([0.1, 0.2])
        result = PredictionResult(
            predictions=preds,
            uncertainties=uncs,
            unit="eV",
            compositions=["NaCl", "TiO2"],
        )

        d = result.to_dict()
        assert "predictions" in d
        assert "uncertainties" in d
        assert "unit" in d
        assert "compositions" in d
        assert d["predictions"] == [1.0, 2.0]


class ConcretePredictor(BasePropertyPredictor):
    """Concrete implementation for testing."""

    @property
    def supported_properties(self) -> list[str]:
        """Return supported properties."""
        return ["test_property"]

    def predict(
        self,
        compositions: str | list[str],
        return_uncertainty: bool = False,
    ) -> np.ndarray | PredictionResult:
        """Make predictions."""
        comps = self._validate_compositions(compositions, validate_smact=False)
        preds = np.ones(len(comps))
        if return_uncertainty:
            return PredictionResult(predictions=preds, compositions=comps)
        return preds


class TestBasePropertyPredictor:
    """Tests for the BasePropertyPredictor class."""

    def test_init_valid_property(self):
        """Should initialise with valid property."""
        predictor = ConcretePredictor(property_name="test_property")
        assert predictor.property_name == "test_property"
        assert predictor.device == "cpu"

    def test_init_invalid_property_raises(self):
        """Should raise error for invalid property."""
        with pytest.raises(ValueError, match="not supported"):
            ConcretePredictor(property_name="invalid_property")

    def test_init_with_fidelity(self):
        """Should store fidelity parameter."""
        predictor = ConcretePredictor(property_name="test_property", fidelity="pbe")
        assert predictor.fidelity == "pbe"

    def test_init_with_model_name(self):
        """Should store model name."""
        predictor = ConcretePredictor(property_name="test_property", model_name="test-model")
        assert predictor.model_name == "test-model"


class TestValidateCompositions:
    """Tests for the _validate_compositions method."""

    def test_single_string(self):
        """Should convert single string to list."""
        predictor = ConcretePredictor(property_name="test_property")
        result = predictor._validate_compositions("NaCl", validate_smact=False)
        assert result == ["NaCl"]

    def test_list_of_strings(self):
        """Should keep list as is."""
        predictor = ConcretePredictor(property_name="test_property")
        result = predictor._validate_compositions(["NaCl", "TiO2"], validate_smact=False)
        assert result == ["NaCl", "TiO2"]

    def test_empty_list_raises(self):
        """Should raise error for empty list."""
        predictor = ConcretePredictor(property_name="test_property")
        with pytest.raises(ValueError, match="(?i)at least one"):
            predictor._validate_compositions([], validate_smact=False)

    def test_smact_validation(self):
        """Should validate compositions using SMACT when enabled."""
        # Test that validation runs without crashing on valid composition
        # Skip if SMACT has compatibility issues (Python version, etc.)
        try:
            from smact.screening import smact_validity

            # First check if smact_validity works at all
            smact_validity("NaCl")
        except (ImportError, ModuleNotFoundError) as e:
            pytest.skip(f"SMACT validation not available in this environment: {e}")

        predictor = ConcretePredictor(property_name="test_property")
        result = predictor._validate_compositions("NaCl", validate_smact=True)
        assert result == ["NaCl"]

    def test_invalid_smact_raises(self):
        """Should raise error for SMACT-invalid compositions."""
        # Skip if SMACT has compatibility issues
        try:
            from smact.screening import smact_validity

            # First check if smact_validity works at all
            smact_validity("NaCl")
        except (ImportError, ModuleNotFoundError) as e:
            pytest.skip(f"SMACT validation not available in this environment: {e}")

        predictor = ConcretePredictor(property_name="test_property")
        # This should fail SMACT validation (nonsensical composition)
        with pytest.raises(ValueError, match="Invalid compositions"):
            predictor._validate_compositions("XxYyZz999", validate_smact=True)


class TestMetadata:
    """Tests for metadata property."""

    def test_metadata_default_empty(self):
        """Metadata should be empty dict by default."""
        predictor = ConcretePredictor(property_name="test_property")
        assert predictor.metadata == {}

    def test_metadata_readonly(self):
        """Metadata property should return the stored dict."""
        predictor = ConcretePredictor(property_name="test_property")
        predictor._metadata = {"key": "value"}
        assert predictor.metadata == {"key": "value"}
