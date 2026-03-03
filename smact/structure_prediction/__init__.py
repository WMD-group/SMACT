"""Minimalist ionic compound prediction tools for materials design."""

from __future__ import annotations

import logging

__all__ = ["database", "mutation", "prediction", "probability_models", "structure", "utilities"]

__author__ = "Alexander Moriarty"
__credits__ = {
    "WMD Group",
    "Imperial College London",
    "Andrew Jackson",
    "Dan Davies",
    "Keith Butler",
    "Aron Walsh",
    "Alexander Moriarty",
}
__status__ = "Development"

logger = logging.getLogger(__name__)
