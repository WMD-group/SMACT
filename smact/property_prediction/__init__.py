"""Property prediction module for materials using pre-trained models.

This module provides utilities for predicting material properties
from chemical composition using pre-trained ROOST models.

Example:
    >>> from smact.property_prediction import predict_band_gap
    >>> predictions = predict_band_gap(["NaCl", "TiO2", "GaN"])

    >>> from smact.property_prediction import RoostPropertyPredictor
    >>> predictor = RoostPropertyPredictor(property_name="band_gap")
    >>> result = predictor.predict(["GaN"], return_uncertainty=True)
"""

from __future__ import annotations

import logging

from smact.property_prediction.base_predictor import (
    BasePropertyPredictor,
    PredictionResult,
)
from smact.property_prediction.convenience import predict_band_gap
from smact.property_prediction.io import clear_cache, list_cached_models
from smact.property_prediction.registry import (
    get_available_models,
    get_property_fidelities,
    get_property_unit,
    get_supported_properties,
)
from smact.property_prediction.roost import RoostPropertyPredictor

__all__ = [
    "BasePropertyPredictor",
    "PredictionResult",
    "RoostPropertyPredictor",
    "clear_cache",
    "get_available_models",
    "get_property_fidelities",
    "get_property_unit",
    "get_supported_properties",
    "list_cached_models",
    "predict_band_gap",
]

__authors__ = (
    "WMD Group",
    "Imperial College London",
    "Ahmed Ismail",
    "Matthew Walker",
    "Jamie Swaine",
    "Hyunsoo Park",
    "Ry Nduma",
)
__credits__ = "Rhys Goodall"
__maintainer__ = "SMACT Development Team"
__status__ = "Development"

logger = logging.getLogger(__name__)
