# SMACT Property Prediction Module

## Module Structure

```text
smact/property_prediction/
├── __init__.py              # Main module interface
├── base_predictor.py        # Abstract base class
├── convenience.py           # Convenience functions
├── roost/
│   ├── __init__.py         # ROOST module
│   └── predictor.py        # ROOST implementation
└── example_usage.py        # Example usage demonstration
```

## API Implementation

Your exact usage pattern now works:

```python
from smact.property_prediction import RoostPropertyPredictor

# Check available properties
print(RoostPropertyPredictor.available_properties)  # ['band_gap', 'bulk_modulus']

# Create model with just property name and device
model = RoostPropertyPredictor(property_name="band_gap", device="cpu")

# Predict single or multiple compositions
model.predict("NaCl")  # Single composition
model.predict(["NaCl", "TiO2", "GaN"])  # Multiple compositions
```

## Key Features

1. **Clean API**: No model_path required - handled internally with defaults
2. **Class-level properties**: `RoostPropertyPredictor.available_properties`
3. **SMACT validation**: All compositions validated using `smact_validity()`
4. **Extensible design**: Easy to add new model types (Wren, CGCNN, etc.)
5. **Type hints**: Full typing support throughout
6. **Error handling**: Clear error messages for invalid compositions/properties

## Extensibility for Multiple Models

The architecture is designed for multiple models:

```python
# Future implementations will work the same way:
from smact.property_prediction import WrenPropertyPredictor, CGCNNPropertyPredictor

wren_model = WrenPropertyPredictor(property_name="band_gap")
cgcnn_model = CGCNNPropertyPredictor(property_name="bulk_modulus")
```

## Abstract Base Class Architecture

In `/home/ryan/smact4properties/SMACT/smact/property_prediction/base_predictor.py`:

```python
from abc import ABC, abstractmethod


class BasePropertyPredictor(ABC):  # Abstract base class
    """Abstract base class for property predictors."""

    @property
    @abstractmethod
    def supported_properties(self) -> list[str]:  # Must be implemented
        """List of properties supported by this predictor."""
        pass

    @abstractmethod
    def predict(
        self, compositions: str | list[str]
    ) -> np.ndarray:  # Must be implemented
        """Predict property values for given compositions."""
        pass
```

## Concrete Implementation

The `RoostPropertyPredictor` **inherits from and implements** the abstract base class:

```python
class RoostPropertyPredictor(BasePropertyPredictor):  # Concrete implementation

    @property
    def supported_properties(self) -> list[str]:  # Implements abstract method
        return self.available_properties

    def predict(
        self, compositions: str | list[str]
    ) -> np.ndarray:  # Implements abstract method
        compositions = self._validate_compositions(compositions)
        # Implementation here...
        return np.zeros(len(compositions))  # Placeholder for now
```

## Extensible Pattern

This allows for multiple concrete implementations:

```python
# Future models will also inherit from BasePropertyPredictor
class WrenPropertyPredictor(BasePropertyPredictor):
    def supported_properties(self) -> list[str]:
        return ["band_gap", "formation_energy"]

    def predict(self, compositions: str | list[str]) -> np.ndarray:
        # Wren-specific implementation
        pass


class CGCNNPropertyPredictor(BasePropertyPredictor):
    def supported_properties(self) -> list[str]:
        return ["bulk_modulus", "shear_modulus"]

    def predict(self, compositions: str | list[str]) -> np.ndarray:
        # CGCNN-specific implementation
        pass
```

## Current Status

- Clean API matching your specification
- SMACT composition validation
- Extensible architecture for multiple models
- Proper error handling and type hints
- Abstract base class enforcing interface contract
- Model loading logic ready for actual ROOST integration
- Placeholder predictions (returns zeros) until real models added

The abstract base class enforces the interface contract while allowing flexibility in
implementation details for each model type. This enables the clean, extensible API and
ensures all future predictors follow the same pattern.

The module is ready for integration with actual ROOST model checkpoints.
