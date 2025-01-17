"""Utility functions for handling intermetallic compounds in SMACT."""

from __future__ import annotations

import numpy as np
from pymatgen.core import Composition

import smact
from smact import Element
from smact.properties import valence_electron_count


def _ensure_composition(composition: str | Composition) -> Composition:
    """Convert input to pymatgen Composition if it isn't already.

    Args:
        composition: Chemical formula as string or pymatgen Composition

    Returns:
        Composition: A pymatgen Composition object

    Raises:
        ValueError: If the composition string is empty
    """
    if isinstance(composition, str):
        if not composition.strip():
            raise ValueError("Empty composition")
        return Composition(composition)
    return composition


def get_element_fraction(composition: str | Composition, element_set: set[str]) -> float:
    """Calculate the fraction of elements from a given set in a composition.
    This helper function is used to avoid code duplication in functions that
    calculate fractions of specific element types (e.g., metals, d-block elements).

    Args:
        composition: Chemical formula as string or pymatgen Composition
        element_set: Set of element symbols to check for

    Returns:
        float: Fraction of the composition that consists of elements from the set (0-1)
    """
    comp = _ensure_composition(composition)
    total_amt = sum(comp.values())
    target_amt = sum(amt for el, amt in comp.items() if el.symbol in element_set)
    return target_amt / total_amt


def get_metal_fraction(composition: str | Composition) -> float:
    """Calculate the fraction of metallic elements in a composition.

    Implemented using get_element_fraction helper with smact.metals set.
    """
    return get_element_fraction(composition, smact.metals)


def get_d_electron_fraction(composition: str | Composition) -> float:
    """Calculate the fraction of d-block elements in a composition.

    Implemented using get_element_fraction helper with smact.d_block set.
    """
    return get_element_fraction(composition, smact.d_block)


def get_distinct_metal_count(composition: str | Composition) -> int:
    """Count the number of distinct metallic elements in a composition."""
    comp = _ensure_composition(composition)
    return sum(1 for el in comp.elements if el.symbol in smact.metals)


def get_pauling_test_mismatch(composition: str | Composition) -> float:
    """Calculate a score for how much the composition deviates from ideal
    Pauling electronegativity ordering.

    Higher mismatch => more difference (ionic, e.g. NaCl).
    Lower mismatch => metal-metal bonds (e.g. Fe-Al).

    Returns:
        float: Mismatch score (0=perfect match, higher=more deviation, NaN=missing data)
    """
    comp = _ensure_composition(composition)
    elements = [Element(el.symbol) for el in comp.elements]
    electronegativities = [el.pauling_eneg for el in elements]

    # If any element lacks a known electronegativity, return NaN
    if None in electronegativities:
        return float("nan")
    else:
        mismatches = []
        for i, (_el1, eneg1) in enumerate(zip(elements, electronegativities, strict=False)):
            for _el2, eneg2 in zip(elements[i + 1 :], electronegativities[i + 1 :], strict=False):
                # Always use absolute difference
                mismatch = abs(eneg1 - eneg2)
                mismatches.append(mismatch)

        # Return average mismatch across all unique pairs
        return np.mean(mismatches) if mismatches else 0.0


def intermetallic_score(composition: str | Composition) -> float:
    """Calculate a score (0-1) indicating how intermetallic a composition is.

    1. Fraction of metallic elements
    2. Number of distinct metals
    3. d-electron fraction
    4. Electronegativity difference (Pauling mismatch)
    5. Valence electron count proximity to 8

    Args:
        composition: Chemical formula or pymatgen Composition
    """
    comp = _ensure_composition(composition)

    # Basic metrics
    metal_fraction = get_metal_fraction(comp)
    d_electron_fraction = get_d_electron_fraction(comp)
    n_metals = get_distinct_metal_count(comp)

    # Valence electron count factor
    try:
        vec = valence_electron_count(comp.reduced_formula)
        vec_factor = 1.0 - abs(vec - 8.0) / 8.0
    except ValueError:
        vec_factor = 0.5

    # Pauling mismatch => large => penalize
    pauling_mismatch = get_pauling_test_mismatch(comp)
    if np.isnan(pauling_mismatch):
        pauling_term = 0.5
    else:
        scale = 3.0
        penalty = min(pauling_mismatch / scale, 1.0)
        pauling_term = 1.0 - penalty

    # Weighted sum
    weights = {
        "metal_fraction": 0.3,
        "d_electron": 0.2,
        "n_metals": 0.2,
        "vec": 0.15,
        "pauling": 0.15,
    }
    score = (
        weights["metal_fraction"] * metal_fraction
        + weights["d_electron"] * d_electron_fraction
        + weights["n_metals"] * min(n_metals / 3.0, 1.0)
        + weights["vec"] * vec_factor
        + weights["pauling"] * pauling_term
    )
    return max(0.0, min(1.0, score))
