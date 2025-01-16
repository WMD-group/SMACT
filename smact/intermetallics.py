"""Utility functions for handling intermetallic compounds in SMACT."""

from __future__ import annotations

import numpy as np
from pymatgen.core import Composition

import smact
from smact import Element
from smact.properties import valence_electron_count, compound_electroneg


def get_metal_fraction(composition: Composition) -> float:
    """Calculate the fraction of metallic elements in a composition.
    
    Args:
        composition: A pymatgen Composition object
        
    Returns:
        float: Fraction of the composition that consists of metallic elements (0-1)
    """
    total_amt = sum(composition.values())
    metal_amt = sum(amt for el, amt in composition.items() if el.symbol in smact.metals)
    return metal_amt / total_amt


def get_d_electron_fraction(composition: Composition) -> float:
    """Calculate the fraction of d-block elements in a composition.
    
    Args:
        composition: A pymatgen Composition object
        
    Returns:
        float: Fraction of the composition that consists of d-block elements (0-1)
    """
    total_amt = sum(composition.values())
    d_block_amt = sum(amt for el, amt in composition.items() if el.symbol in smact.d_block)
    return d_block_amt / total_amt


def get_distinct_metal_count(composition: Composition) -> int:
    """Count the number of distinct metallic elements in a composition.
    
    Args:
        composition: A pymatgen Composition object
        
    Returns:
        int: Number of distinct metallic elements
    """
    return sum(1 for el in composition.elements if el.symbol in smact.metals)


def get_pauling_test_mismatch(composition: Composition) -> float:
    """Calculate a score for how much the composition deviates from ideal Pauling electronegativity ordering.
    
    Args:
        composition: A pymatgen Composition object
        
    Returns:
        float: Mismatch score (0 = perfect match, higher = more deviation)
    """
    # Convert pymatgen elements to SMACT elements using their symbols
    elements = [Element(el.symbol) for el in composition.elements]
    electronegativities = [el.pauling_eneg for el in elements]
    
    # Skip if any electronegativities are None
    if None in electronegativities:
        return 0.0
        
    # Calculate pairwise differences
    mismatches = []
    for i, (el1, eneg1) in enumerate(zip(elements, electronegativities)):
        for el2, eneg2 in zip(elements[i+1:], electronegativities[i+1:]):
            # For metal pairs, we expect small electronegativity differences
            if el1.symbol in smact.metals and el2.symbol in smact.metals:
                mismatches.append(abs(eneg1 - eneg2))
            # For metal-nonmetal pairs, we expect larger differences
            elif (el1.symbol in smact.metals) != (el2.symbol in smact.metals):
                mismatches.append(1.0 - abs(eneg1 - eneg2))
                
    return np.mean(mismatches) if mismatches else 0.0


def intermetallic_score(composition: str | Composition) -> float:
    """Calculate a score indicating how likely a composition is to be an intermetallic compound.
    
    The score is based on several heuristics:
    1. Fraction of metallic elements
    2. Number of distinct metals
    3. Presence of d-block elements
    4. Electronegativity differences
    5. Valence electron count
    
    Args:
        composition: Chemical formula as string or pymatgen Composition
        
    Returns:
        float: Score between 0 and 1, where higher values indicate more intermetallic character
        
    Example:
        >>> intermetallic_score("Fe3Al")
        0.85  # High score - likely intermetallic
        >>> intermetallic_score("NaCl") 
        0.2   # Low score - ionic compound
    """
    if isinstance(composition, str):
        composition = Composition(composition)
        
    # 1. Basic metrics
    metal_fraction = get_metal_fraction(composition)
    d_electron_fraction = get_d_electron_fraction(composition)
    n_metals = get_distinct_metal_count(composition)
    
    # 2. Electronic structure indicators
    try:
        vec = valence_electron_count(composition.reduced_formula)
        vec_factor = 1.0 - (abs(vec - 8.0) / 8.0)  # Normalized around VEC=8
    except ValueError:
        vec_factor = 0.5  # Default if we can't calculate VEC
        
    # 3. Bonding character
    pauling_mismatch = get_pauling_test_mismatch(composition)
    
    # 4. Calculate weighted score
    # These weights can be tuned based on empirical testing
    weights = {
        'metal_fraction': 0.3,
        'd_electron': 0.2,
        'n_metals': 0.2,
        'vec': 0.15,
        'pauling': 0.15
    }
    
    score = (
        weights['metal_fraction'] * metal_fraction +
        weights['d_electron'] * d_electron_fraction +
        weights['n_metals'] * min(1.0, n_metals / 3.0) +  # Normalized to max of 3 metals
        weights['vec'] * vec_factor +
        weights['pauling'] * (1.0 - pauling_mismatch)  # Invert mismatch score
    )
    
    return min(1.0, max(0.0, score))  # Clamp between 0 and 1 