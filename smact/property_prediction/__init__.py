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
from typing import TYPE_CHECKING

from smact.property_prediction.base_predictor import (
    BasePropertyPredictor,
    PredictionResult,
)
from smact.property_prediction.io import clear_cache, list_cached_models
from smact.property_prediction.registry import (
    get_available_models,
    get_property_fidelities,
    get_property_unit,
    get_supported_properties,
)

if TYPE_CHECKING:
    from smact.property_prediction.convenience import predict_band_gap
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

_LAZY_IMPORTS = {
    "RoostPropertyPredictor": "smact.property_prediction.roost",
    "predict_band_gap": "smact.property_prediction.convenience",
}


def __getattr__(name: str):
    """Lazy import for torch-dependent symbols."""
    if name in _LAZY_IMPORTS:
        import importlib

        module = importlib.import_module(_LAZY_IMPORTS[name])
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
