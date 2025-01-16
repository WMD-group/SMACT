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
    "Mg2Si",
    "Cr3Si",
    "CoAl",
    "Fe3C",
    "NbTiAlCr",
    "Cu2MgSn",
]

print("Testing SMACT validity with intermetallic detection:")
print("-" * 50)

for compound in compounds:
    # Test with standard validity check
    is_valid_standard = smact_validity(compound, include_alloys=False)

    # Test with intermetallic detection
    is_valid_intermetallic = smact_validity(compound, include_alloys=True)

    # Get detailed intermetallic score
    score = intermetallic_score(compound)

    print(f"\nCompound: {compound}")
    print(f"Standard validity: {is_valid_standard}")
    print(f"With intermetallic detection: {is_valid_intermetallic}")
    print(f"Intermetallic score: {score:.2f}")
