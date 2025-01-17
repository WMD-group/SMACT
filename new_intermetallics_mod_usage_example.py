"""Example usage of SMACT's intermetallic functionality."""

from __future__ import annotations

from smact.intermetallics import intermetallic_score
from smact.screening import smact_validity

# Test some known compounds
compounds = [
    "Fe3Al",  # Classic intermetallic
    "Ni3Ti",  # Superalloy
    "NaCl",  # Ionic compound
    "Fe2O3",  # Metal oxide
    "Mg2Si",  # Semiconductor
    "Cr3Si",  # Intermetallic
    "CoAl",  # Intermetallic
    "Fe3C",  # Interstitial compound
    "NbTiAlCr",  # High-entropy alloy
    "Cu2MgSn",  # Heusler alloy
    "AcFe3",  # Contains element with missing electronegativity
    "Cu",  # Single element (metal)
    "Si",  # Single element (non-metal)
]

print("Testing SMACT validity with different validation methods:")
print("-" * 50)

for compound in compounds:
    # Test with different validation methods
    is_valid_standard = smact_validity(compound, use_pauling_test=True, include_alloys=False, check_intermetallic=False)

    is_valid_alloy = smact_validity(compound, use_pauling_test=True, include_alloys=True, check_intermetallic=False)

    is_valid_intermetallic = smact_validity(
        compound, use_pauling_test=True, include_alloys=False, check_intermetallic=True, intermetallic_threshold=0.7
    )

    # Get detailed intermetallic score
    score = intermetallic_score(compound)

    print(f"\nCompound: {compound}")
    print(f"Standard validity (no alloys/intermetallics): {is_valid_standard}")
    print(f"With alloy detection: {is_valid_alloy}")
    print(f"With intermetallic detection: {is_valid_intermetallic}")
    print(f"Intermetallic score: {score:.2f}")

# Test some edge cases
print("\nTesting edge cases:")
print("-" * 50)

edge_cases = [
    "Fe2O3",  # Mixed metal-nonmetal (should be valid standard)
    "CuZn",  # Simple binary alloy (should be valid with alloys)
    "Fe3Al",  # Intermetallic (should be valid with intermetallics)
    "NaCl",  # Ionic (should be valid standard)
    "SiO2",  # Covalent (should be valid standard)
]

for compound in edge_cases:
    # Test all validation methods together
    is_valid = smact_validity(
        compound, use_pauling_test=True, include_alloys=True, check_intermetallic=True, intermetallic_threshold=0.7
    )

    score = intermetallic_score(compound)

    print(f"\nCompound: {compound}")
    print(f"Valid by any method: {is_valid}")
    print(f"Intermetallic score: {score:.2f}")
