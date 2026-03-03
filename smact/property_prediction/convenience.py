"""Convenience functions for property prediction.

This module provides simple wrapper functions for common property
predictions, allowing users to predict properties with minimal setup.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from smact.property_prediction.base_predictor import PredictionResult


def predict_band_gap(
    compositions: str | list[str],
    model_path: str | Path | None = None,
    return_uncertainty: bool = False,
    device: str = "cpu",
) -> np.ndarray | PredictionResult:
    """Predict band gap for given compositions using ROOST.

    Uses a pre-trained ROOST model (0.28 eV MAE on Materials Project
    103K DFT band gaps) to predict band gaps from composition alone.

    Args:
        compositions: Single composition string or list of composition
            strings (e.g., "NaCl", ["TiO2", "GaN"]).
        model_path: Path to a local model checkpoint. If None, uses the
            default pre-trained model.
        return_uncertainty: If True, return PredictionResult with
            uncertainty estimates.
        device: Device to run the model on ("cpu" or "cuda").

    Returns:
        Band gap predictions in eV as numpy array, or PredictionResult
        if return_uncertainty is True.

    Example:
        >>> predict_band_gap("GaN")
        array([1.59])

        >>> result = predict_band_gap(["NaCl", "TiO2"], return_uncertainty=True)
        >>> print(result.predictions, result.uncertainties)
    """
    from smact.property_prediction.roost import RoostPropertyPredictor

    predictor = RoostPropertyPredictor(
        property_name="band_gap",
        model_path=model_path,
        device=device,
    )
    return predictor.predict(compositions, return_uncertainty=return_uncertainty)
