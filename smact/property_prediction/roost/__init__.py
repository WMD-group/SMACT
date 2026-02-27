"""ROOST-based property prediction."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from smact.property_prediction.roost.predictor import RoostPropertyPredictor


def __getattr__(name: str):
    """Lazy import to avoid requiring torch at package import time."""
    if name in ("RoostPropertyPredictor", "Roost"):
        from smact.property_prediction.roost.predictor import RoostPropertyPredictor

        if name == "Roost":
            return RoostPropertyPredictor
        return RoostPropertyPredictor

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["RoostPropertyPredictor"]
