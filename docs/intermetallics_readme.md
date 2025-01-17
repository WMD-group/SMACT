# SMACT Intermetallics Enhancement

This document describes the enhancements made to SMACT to better handle intermetallic compounds and alloys. The changes introduce a scoring system for identifying and validating intermetallic compounds, moving beyond the simple binary classification previously used.

## Key Changes

### 1. New Module: `smact.intermetallics`

A dedicated module for handling intermetallic compounds with several specialized functions:

#### `get_metal_fraction(composition)`

- Calculates the fraction of metallic elements in a composition
- Input: Chemical formula (str) or pymatgen Composition object
- Output: Float between 0-1
- Example:

```python
from smact.intermetallics import get_metal_fraction

# Pure intermetallic - returns 1.0
print(get_metal_fraction("Fe3Al"))  # Works with string formula (& Composition object)

# Mixed compound - returns fraction of how much of the compound is metallic
print(get_metal_fraction("Fe2O3"))  # 0.4
```

#### `get_d_electron_fraction(composition)`

- Calculates the fraction of d-block elements
- Important for transition metal intermetallics
- Example:

```python
from smact.intermetallics import get_d_electron_fraction

# Pure transition metal compound - returns 1.0
print(get_d_electron_fraction("Fe2Nb"))

# Mixed compound - returns fraction
print(get_d_electron_fraction("Fe3Al"))  # 0.75 (Fe is d-block, Al is not)
```

#### `get_distinct_metal_count(composition)`

- Counts unique metallic elements
- Useful for identifying complex intermetallics
- Example:

```python
from smact.intermetallics import get_distinct_metal_count

# Binary intermetallic
print(get_distinct_metal_count("Fe3Al"))  # Returns 2

# Complex HEA-like composition
print(get_distinct_metal_count("NbTiAlCr"))  # Returns 4

# Non-metallic compound
print(get_distinct_metal_count("SiO2"))  # Returns 0
```

#### `get_pauling_test_mismatch(composition)`

- Calculates deviation from ideal electronegativity ordering
- Handles both metal-metal and metal-nonmetal pairs
- Lower scores indicate more intermetallic-like bonding
- Example:

```python
from smact.intermetallics import get_pauling_test_mismatch

# Intermetallic - low mismatch
print(get_pauling_test_mismatch("Fe3Al"))  # 0.22

# Ionic compound - high mismatch
print(get_pauling_test_mismatch("NaCl"))  # -1.23
```

#### `intermetallic_score(composition)`

- Main scoring function combining multiple metrics
- Returns a score between 0-1
- Higher scores indicate more intermetallic character depending on the set threshold and weighting of the metrics
- Example:

```python
from smact.intermetallics import intermetallic_score

# Classic intermetallics - high scores
print(intermetallic_score("Fe3Al"))  # ~0.83
print(intermetallic_score("Ni3Ti"))  # ~0.87

# Non-intermetallics - low scores
print(intermetallic_score("NaCl"))  # 0.63
print(intermetallic_score("Fe2O3"))  # 0.64
print(intermetallic_score("SiO2"))  # 0.25
```

### 2. Enhanced `smact_validity`

The existing `smact_validity` function in `smact.screening` has been enhanced with three validation paths:

1. Standard validation (charge neutrality and electronegativity)
2. Simple metal alloy validation
3. Intermetallic scoring validation

Example usage showing all validation paths:

```python
from smact.screening import smact_validity

# Test with different validation methods
compound = "Fe3Al"

# 1. Standard validation (no alloys/intermetallics)
is_valid_standard = smact_validity(
    compound, use_pauling_test=True, include_alloys=False, check_intermetallic=False
)

# 2. With alloy detection
is_valid_alloy = smact_validity(
    compound, use_pauling_test=True, include_alloys=True, check_intermetallic=False
)

# 3. With intermetallic detection
is_valid_intermetallic = smact_validity(
    compound,
    use_pauling_test=True,
    include_alloys=False,
    check_intermetallic=True,
    intermetallic_threshold=0.7,
)

# Or combine methods
is_valid = smact_validity(
    compound,
    use_pauling_test=True,
    include_alloys=True,
    check_intermetallic=True,
    intermetallic_threshold=0.7,
)
```

### 3. Comprehensive Analysis Example

Here's how to perform a detailed analysis of a compound:

```python
from smact.intermetallics import *
from smact.screening import smact_validity


def analyze_compound(formula):
    """Perform comprehensive analysis of a compound."""
    # Basic metrics
    metal_frac = get_metal_fraction(formula)
    d_frac = get_d_electron_fraction(formula)
    n_metals = get_distinct_metal_count(formula)
    pauling = get_pauling_test_mismatch(formula)
    score = intermetallic_score(formula)

    # Validity checks
    valid_standard = smact_validity(
        formula, use_pauling_test=True, include_alloys=False, check_intermetallic=False
    )
    valid_alloy = smact_validity(
        formula, use_pauling_test=True, include_alloys=True, check_intermetallic=False
    )
    valid_intermetallic = smact_validity(
        formula, use_pauling_test=True, include_alloys=False, check_intermetallic=True
    )

    print(f"Analysis of {formula}:")
    print(f"Metal fraction: {metal_frac:.2f}")
    print(f"d-electron fraction: {d_frac:.2f}")
    print(f"Distinct metals: {n_metals}")
    print(f"Pauling mismatch: {'nan' if math.isnan(pauling) else f'{pauling:.2f}'}")
    print(f"Intermetallic score: {score:.2f}")
    print(f"Valid (standard): {valid_standard}")
    print(f"Valid (alloy): {valid_alloy}")
    print(f"Valid (intermetallic): {valid_intermetallic}")


# Example usage
compounds = [
    "Fe3Al",  # Classic intermetallic
    "NaCl",  # Ionic
    "Fe2O3",  # Metal oxide
    "Cu2MgSn",  # Heusler alloy
    "NbTiAlCr",  # High-entropy alloy
]

for compound in compounds:
    analyze_compound(compound)
    print("-" * 50)
```

## Differences from Previous Version

### Before

- Simple binary check for all-metal compositions
- No distinction between intermetallics and other metallic phases
- Limited handling of mixed bonding character
- Binary valid/invalid classification

### After

- Scoring system using multiple chemical descriptors
- Better handling of partial metallic character
- Consideration of d-electron contributions
- Continuous scoring (0-1) for more nuanced classification
- Adjustable threshold for different applications and weightings for said combined rule\*\*
  -The intermetallic_threshold in smact_validity() can be raised or lowered. Literature shows that in real materials, the line between "intermetallic" and "ionic/metallic" can be fuzzy. Having a tunable threshold aligns with different research needs (e.g., searching for strongly metallic Heuslers vs. half-metallic systems)

## Usage Examples

### Basic Screening

```python
from smact.screening import smact_validity
from smact.intermetallics import intermetallic_score

compounds = [
    "Fe3Al",  # Classic intermetallic
    "Ni3Ti",  # Superalloy
    "NaCl",  # Ionic
    "Fe2O3",  # Metal oxide
]

for compound in compounds:
    score = intermetallic_score(compound)
    is_valid = smact_validity(compound, include_alloys=True)
    print(f"{compound}: score={score:.2f}, valid={is_valid}")
```

### Advanced Usage

```python
from pymatgen.core import Composition
from smact.intermetallics import *

# Detailed analysis of a compound
comp = Composition("Fe3Al")
metrics = {
    "metal_fraction": get_metal_fraction(comp),
    "d_electron_fraction": get_d_electron_fraction(comp),
    "distinct_metals": get_distinct_metal_count(comp),
    "pauling_mismatch": get_pauling_test_mismatch(comp),
    "overall_score": intermetallic_score(comp),
}
```

## Known Limitations and Pitfalls

1. **Electronegativity Data**

   - Some elements may lack Pauling electronegativity data
   - Falls back to default behavior in these cases

2. **VEC Calculation**

   - Assumes simple electron counting rules
   - May not capture complex electronic structures

3. **Threshold Selection**

   - Default threshold (0.7) may need adjustment for specific applications
   - Consider domain-specific validation

4. **Complex Compositions**
   - High-entropy alloys may need different weighting schemes
   - Current weights optimized for binary/ternary systems

## Future Development Directions

1. **Additional Features**

   - Incorporate atomic size factors
   - Add structure prediction capabilities
   - Include formation energy estimates - linking back to Miedema's model [DOI Link](https://doi.org/10.1016/j.cpc.2016.08.013)

2. **Validation and Refinement**

   - Benchmark against experimental databases
   - Refine scoring weights with more data
   - Add support for mixed-valence compounds

3. **Extended Functionality**
   - Add support for partial occupancies
   - Include temperature-dependent properties
   - Integrate with phase diagram tools

## Contributing

Contributions to improve the intermetallics functionality are welcome! Areas particularly in need of development:

1. Additional test cases and validation
2. Refinement of scoring weights
3. Integration with external databases
4. Performance optimization for large-scale screening

## References

1. Original SMACT paper: [SMACT: Semiconducting Materials by Analogy and Chemical Theory](https://joss.theoj.org/papers/10.21105/joss.01361)
2. Intermetallics theory and classification literature sources:

   - D.G. Pettifor introduced the concept of a single "chemical scale" or "structure map" coordinate (Pettifor number) to systematically separate compound classes [1, p. 31]. The new intermetallicscore is a step in that direction but customized to SMACT's internal data structures.

     - Reference: D.G. Pettifor, "A chemical scale for crystal-structure maps," Solid State Communications. 51 (1984) 31–34. [DOI Link](<https://doi.org/10.1016/0038-1098(84)90765-8>)

   - The role of charge transfer and atomic size mismatch is pivotal in stabilizing intermetallic phases. Miedema's framework quantifies these effects, making it useful for predicting alloying behaviors and crystal structure. The parameters coded here, while conceptually similar, have not implemented Miedema's model directly.
   - Reference: A.R. Miedema, Cohesion in alloys - fundamentals of a semi-empirical model.[DOI Link](<https://doi.org/10.1016/0378-4363(80)90054-6>)

3. Electronegativity scales (pauling electronegativity)
