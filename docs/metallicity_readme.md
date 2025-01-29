# SMACT Metallicity Module

This document describes the enhancements made to SMACT to better handle metallic compounds and their classification. The changes introduce a scoring system for identifying and validating metallic compounds, moving beyond simple binary classification to provide a more nuanced understanding of a compound's metallic character.

## Key Changes

### 1. Module: `smact.metallicity`

A dedicated module for handling metallic compounds with several specialized functions:

#### `get_metal_fraction(composition)`

- Calculates the fraction of metallic elements in a composition
- Input: Chemical formula (str) or pymatgen Composition object
- Output: Float between 0-1
- Example:

```python
from smact.metallicity import get_metal_fraction

# Pure metallic compound - returns 1.0
print(get_metal_fraction("Fe3Al"))  # Works with string formula (& Composition object)

# Mixed compound - returns fraction of how much of the compound is metallic
print(get_metal_fraction("Fe2O3"))  # 0.4
```

#### `get_d_electron_fraction(composition)`

- Calculates the fraction of d-block elements
- Important for transition metal compounds
- Example:

```python
from smact.metallicity import get_d_electron_fraction

# Pure transition metal compound - returns 1.0
print(get_d_electron_fraction("Fe2Nb"))

# Mixed compound - returns fraction
print(get_d_electron_fraction("Fe3Al"))  # 0.75 (Fe is d-block, Al is not)
```

#### `get_distinct_metal_count(composition)`

- Counts unique metallic elements
- Useful for identifying complex metallic systems
- Example:

```python
from smact.metallicity import get_distinct_metal_count

# Binary metallic compound
print(get_distinct_metal_count("Fe3Al"))  # Returns 2

# Complex multi-component composition
print(get_distinct_metal_count("NbTiAlCr"))  # Returns 4

# Non-metallic compound
print(get_distinct_metal_count("SiO2"))  # Returns 0
```

#### `get_pauling_test_mismatch(composition)`

- Calculates deviation from ideal electronegativity ordering
- Handles both metal-metal and metal-nonmetal pairs
- Lower scores indicate more metallic-like bonding
- Example:

```python
from smact.metallicity import get_pauling_test_mismatch

# Metallic compound - low mismatch
print(get_pauling_test_mismatch("Fe3Al"))  # 0.22

# Ionic compound - high mismatch
print(get_pauling_test_mismatch("NaCl"))  # -1.23
```

#### `metallicity_score(composition)`

- Main scoring function combining multiple metrics
- Returns a score between 0-1
- Higher scores indicate more metallic character
- Example:

```python
from smact.metallicity import metallicity_score

# Metallic compounds - high scores
print(metallicity_score("Fe3Al"))  # ~0.83
print(metallicity_score("Ni3Ti"))  # ~0.87

# Non-metallic compounds - low scores
print(metallicity_score("NaCl"))  # 0.63
print(metallicity_score("Fe2O3"))  # 0.64
print(metallicity_score("SiO2"))  # 0.25
```

### 2. Enhanced `smact_validity`

The existing `smact_validity` function in `smact.screening` has been enhanced with three validation paths:

1. Standard validation (charge neutrality and electronegativity)
2. Simple metal validation
3. Metallicity scoring validation

Example usage showing all validation paths:

```python
from smact.screening import smact_validity

# Test with different validation methods
compound = "Fe3Al"

# 1. Standard validation (no alloy/metallicity check)
is_valid_standard = smact_validity(
    compound, use_pauling_test=True, include_alloys=False, check_metallicity=False
)

# 2. With alloy detection
is_valid_alloy = smact_validity(
    compound, use_pauling_test=True, include_alloys=True, check_metallicity=False
)

# 3. With metallicity detection
is_valid_metallic = smact_validity(
    compound,
    use_pauling_test=True,
    include_alloys=False,
    check_metallicity=True,
    metallicity_threshold=0.7,
)

# Or combine methods
is_valid = smact_validity(
    compound,
    use_pauling_test=True,
    include_alloys=True,
    check_metallicity=True,
    metallicity_threshold=0.7,
)
```

### 3. Comprehensive Analysis Example

Here's how to perform a detailed analysis of a compound:

```python
from smact.metallicity import *
from smact.screening import smact_validity


def analyze_compound(formula):
    """Perform comprehensive analysis of a compound."""
    # Basic metrics
    metal_frac = get_metal_fraction(formula)
    d_frac = get_d_electron_fraction(formula)
    n_metals = get_distinct_metal_count(formula)
    pauling = get_pauling_test_mismatch(formula)
    score = metallicity_score(formula)

    # Validity checks
    valid_standard = smact_validity(
        formula, use_pauling_test=True, include_alloys=False, check_metallicity=False
    )
    valid_alloy = smact_validity(
        formula, use_pauling_test=True, include_alloys=True, check_metallicity=False
    )
    valid_metallic = smact_validity(
        formula, use_pauling_test=True, include_alloys=False, check_metallicity=True
    )

    print(f"Analysis of {formula}:")
    print(f"Metal fraction: {metal_frac:.2f}")
    print(f"d-electron fraction: {d_frac:.2f}")
    print(f"Distinct metals: {n_metals}")
    print(f"Pauling mismatch: {'nan' if math.isnan(pauling) else f'{pauling:.2f}'}")
    print(f"Metallicity score: {score:.2f}")
    print(f"Valid (standard): {valid_standard}")
    print(f"Valid (alloy): {valid_alloy}")
    print(f"Valid (metallic): {valid_metallic}")


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

## Machine Learning Classification

The metallicity module has been used in conjunction with machine learning approaches to predict the metallic nature of compounds. The `metallicity_classification.ipynb` notebook demonstrates how these chemical descriptors can be used to train models for predicting whether a compound will exhibit metallic behavior.

The classification task uses features derived from the metallicity module including:

- Metal fraction
- d-electron fraction
- Number of distinct metals
- Pauling electronegativity mismatch
- Overall metallicity score

These features are used to train models that can predict the likelihood of metallic behavior in compounds, which is more general than specifically predicting intermetallic behavior. This approach allows for a broader understanding of metallic character in materials, encompassing various types of metallic systems including:

- Traditional metals and alloys
- Intermetallic compounds
- Solid solutions
- Metallic glasses
- Mixed-character compounds

## Future Development Directions

1. **Specialized Submodules**

   - Development of specific modules for different types of metallic systems:
     - Intermetallic compound prediction
     - Solid solution formation
     - Amorphous metal formation
   - Integration with structure prediction

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

Contributions to improve the metallicity functionality are welcome! Areas particularly in need of development:

1. Additional test cases and validation
2. Refinement of scoring weights
3. Integration with external databases
4. Development of specialized submodules for specific metallic systems
5. Performance optimization for large-scale screening

## References

1. Original SMACT paper: [SMACT: Semiconducting Materials by Analogy and Chemical Theory](https://joss.theoj.org/papers/10.21105/joss.01361)
2. Intermetallics theory and classification literature sources:

   - D.G. Pettifor introduced the concept of a single "chemical scale" or "structure map" coordinate (Pettifor number) to systematically separate compound classes [1, p. 31]. The new intermetallicscore is a step in that direction but customized to SMACT's internal data structures.

     - Reference: D.G. Pettifor, "A chemical scale for crystal-structure maps," Solid State Communications. 51 (1984) 31–34. [DOI Link](<https://doi.org/10.1016/0038-1098(84)90765-8>)

   - The role of charge transfer and atomic size mismatch is pivotal in stabilizing intermetallic phases. Miedema's framework quantifies these effects, making it useful for predicting alloying behaviors and crystal structure. The parameters coded here, while conceptually similar, have not implemented Miedema's model directly.
   - Reference: A.R. Miedema, Cohesion in alloys - fundamentals of a semi-empirical model.[DOI Link](<https://doi.org/10.1016/0378-4363(80)90054-6>)

3. Electronegativity scales and their role in predicting bonding character (Pauling electronegativity)
