# SMACT Intermetallics Enhancement

This document describes the enhancements made to SMACT to better handle intermetallic compounds and alloys. The changes introduce a scoring system for identifying and validating intermetallic compounds, moving beyond the simple binary classification previously used.

## Key Changes

### 1. New Module: `smact.intermetallics`

A dedicated module for handling intermetallic compounds with several specialized functions:

#### `get_metal_fraction(composition)`

- Calculates the fraction of metallic elements in a composition
- Input: pymatgen Composition object
- Output: Float between 0-1
- Example:

```python
from pymatgen.core import Composition
from smact.intermetallics import get_metal_fraction

# Pure intermetallic - returns 1.0
print(get_metal_fraction(Composition("Fe3Al")))

# Mixed compound - returns fraction
print(get_metal_fraction(Composition("Fe2O3")))
```

#### `get_d_electron_fraction(composition)`

- Calculates the fraction of d-block elements
- Important for transition metal intermetallics
- Example:

```python
from smact.intermetallics import get_d_electron_fraction

# Pure transition metal compound - returns 1.0
print(get_d_electron_fraction(Composition("Fe2Nb")))

# Main group compound - returns 0.0
print(get_d_electron_fraction(Composition("Mg2Si")))
```

#### `get_distinct_metal_count(composition)`

- Counts unique metallic elements
- Useful for identifying complex intermetallics
- Example:

```python
from smact.intermetallics import get_distinct_metal_count

# Binary intermetallic
print(get_distinct_metal_count(Composition("Fe3Al")))  # Returns 2

# Complex HEA-like composition
print(get_distinct_metal_count(Composition("NbTiAlCr")))  # Returns 4
```

#### `get_pauling_test_mismatch(composition)`

- Calculates deviation from ideal electronegativity ordering
- Handles both metal-metal and metal-nonmetal pairs
- Lower scores indicate more intermetallic-like bonding
- Example:

```python
from smact.intermetallics import get_pauling_test_mismatch

# Intermetallic - low mismatch
print(get_pauling_test_mismatch(Composition("Fe3Al")))

# Ionic compound - high mismatch
print(get_pauling_test_mismatch(Composition("NaCl")))
```

#### `intermetallic_score(composition)`

- Main scoring function combining multiple metrics
- Returns a score between 0-1
- Higher scores indicate more intermetallic character
- Example:

```python
from smact.intermetallics import intermetallic_score

# Classic intermetallics - high scores
print(intermetallic_score("Fe3Al"))  # ~0.85
print(intermetallic_score("Ni3Ti"))  # ~0.82

# Non-intermetallics - low scores
print(intermetallic_score("NaCl"))  # ~0.20
print(intermetallic_score("Fe2O3"))  # ~0.45
```

### 2. Enhanced `smact_validity`

The existing `smact_validity` function in `smact.screening` has been enhanced:

- New parameter `intermetallic_threshold` (default: 0.7) # TODO: change to 0.5 as the classification task suggests that this is the best parameter for the task
- Uses the scoring system when `include_alloys=True`
- More nuanced than previous binary metal check
- Example:

```python
from smact.screening import smact_validity

# Check with intermetallic detection
print(smact_validity("Fe3Al", include_alloys=True))  # True
print(smact_validity("NaCl", include_alloys=True))  # False

# Adjust threshold for stricter filtering
print(smact_validity("Fe3Al", include_alloys=True, intermetallic_threshold=0.8))
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
  -The intermetallic_threshold in smact_validity() can be raised or lowered. Literature shows that in real materials, the line between “intermetallic” and “ionic/metallic” can be fuzzy. Having a tunable threshold aligns with different research needs (e.g., searching for strongly metallic Heuslers vs. half-metallic systems)

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
   - Include formation energy estimates

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
2. Intermetallics theory and classification: Various literature sources
3. Electronegativity scales and their application to intermetallics
