"""Base class for property predictors in SMACT.

This module provides the abstract base class that all property predictors
must inherit from, ensuring a consistent interface across different
prediction models.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class PredictionResult:
    """Container for prediction results with optional uncertainty estimates.

    Attributes:
        predictions: Array of predicted property values.
        uncertainties: Aleatoric (per-sample) uncertainties from robust models.
        epistemic_std: Epistemic uncertainty from ensemble predictions.
        unit: Unit string for the predicted property (e.g., "eV", "GPa").
        compositions: List of composition strings that were predicted.
    """

    predictions: np.ndarray
    uncertainties: np.ndarray | None = None
    epistemic_std: np.ndarray | None = None
    unit: str = ""
    compositions: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        """Return the number of predictions."""
        return len(self.predictions)

    def __getitem__(self, idx: int) -> float:
        """Get prediction at index."""
        return float(self.predictions[idx])

    def __repr__(self) -> str:
        """Return string representation."""
        n = len(self.predictions)
        has_unc = self.uncertainties is not None
        return f"PredictionResult(n={n}, has_uncertainty={has_unc}, unit='{self.unit}')"

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for serialisation."""
        result = {
            "predictions": self.predictions.tolist(),
            "unit": self.unit,
            "compositions": self.compositions,
        }
        if self.uncertainties is not None:
            result["uncertainties"] = self.uncertainties.tolist()
        if self.epistemic_std is not None:
            result["epistemic_std"] = self.epistemic_std.tolist()
        return result


class BasePropertyPredictor(ABC):
    """Abstract base class for property predictors.

    This class defines the interface for all property predictors in SMACT.
    Subclasses must implement the `supported_properties` property and
    the `predict` method for specific models.

    Attributes:
        property_name: Name of the property being predicted.
        fidelity: Fidelity level for the prediction (e.g., "pbe", "hse06").
        model_name: Name of the specific model version.
        model_path: Path to a local model checkpoint.
        device: Device to run the model on ("cpu" or "cuda").
        model: The loaded model instance.
    """

    def __init__(
        self,
        property_name: str,
        fidelity: str | None = None,
        model_name: str | None = None,
        model_path: str | Path | None = None,
        device: str = "cpu",
        **kwargs: Any,
    ):
        """Initialise the property predictor.

        Args:
            property_name: Name of the property to predict
                (e.g., "band_gap", "bulk_modulus").
            fidelity: Fidelity level for properties with multiple models
                (e.g., "pbe", "hse06" for band gap). If None, uses default.
            model_name: Specific model version to use. If None, uses default
                for the property/fidelity combination.
            model_path: Path to a local model checkpoint directory.
                Overrides model_name if provided.
            device: Device to run the model on ("cpu" or "cuda").
            **kwargs: Additional keyword arguments for specific implementations.
        """
        self.property_name = property_name
        self.fidelity = fidelity
        self.model_name = model_name
        self.model_path = Path(model_path) if model_path else None
        self.device = device
        self.model: Any = None
        self._metadata: dict[str, Any] = {}

        # Validate property name
        if property_name not in self.supported_properties:
            raise ValueError(
                f"Property '{property_name}' not supported. "
                f"Supported properties: {', '.join(self.supported_properties)}"
            )

    @property
    @abstractmethod
    def supported_properties(self) -> list[str]:
        """List of properties supported by this predictor."""

    @abstractmethod
    def predict(
        self,
        compositions: str | list[str],
        return_uncertainty: bool = False,
    ) -> np.ndarray | PredictionResult:
        """Predict property values for given compositions.

        Args:
            compositions: Single composition string or list of
                composition strings (e.g., "NaCl", ["TiO2", "GaN"]).
            return_uncertainty: If True, return PredictionResult with
                uncertainty estimates. If False, return numpy array.

        Returns:
            Property predictions as numpy array, or PredictionResult
            object if return_uncertainty is True.
        """

    def _validate_compositions(
        self,
        compositions: str | list[str],
        validate_smact: bool = True,
    ) -> list[str]:
        """Validate and normalise composition inputs.

        Converts input to a list of composition strings and optionally
        validates each composition using SMACT validity checks.

        Args:
            compositions: Single composition string or list of strings.
            validate_smact: Whether to apply SMACT validity checks.
                If True, invalid compositions will raise an error.

        Returns:
            List of validated composition strings.

        Raises:
            ValueError: If validate_smact is True and any composition
                fails SMACT validity checks.
        """
        # Convert to list
        if isinstance(compositions, str):
            compositions = [compositions]

        # Ensure all are strings
        compositions = [str(c) for c in compositions]

        if not compositions:
            raise ValueError("At least one composition must be provided.")

        if validate_smact:
            from smact.screening import smact_validity

            invalid = []
            for comp in compositions:
                try:
                    if not smact_validity(comp):
                        invalid.append(comp)
                except (ValueError, TypeError):
                    # If smact_validity raises an exception, the composition
                    # is likely malformed
                    invalid.append(comp)

            if invalid:
                raise ValueError(
                    f"Invalid compositions according to SMACT: {invalid}. Set validate_smact=False to skip validation."
                )

        return compositions

    @property
    def metadata(self) -> dict[str, Any]:
        """Model metadata including training information."""
        return self._metadata
