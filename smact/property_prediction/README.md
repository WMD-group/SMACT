# SMACT Property Prediction Module

## Module Structure

```text
smact/property_prediction/
├── __init__.py              # Main module interface
├── base_predictor.py        # Abstract base class
├── config.py                # Property metadata and defaults
├── convenience.py           # Convenience functions (predict_band_gap)
├── io.py                    # Model loading, caching, and downloading
├── registry.py              # Model registry and discovery
├── roost/
│   ├── __init__.py         # ROOST module
│   ├── predictor.py        # ROOST implementation
│   └── train.py            # Training pipeline
└── scripts/
    └── convert_checkpoint.py  # Checkpoint conversion utility
```

## Quick Start

```python
from smact.property_prediction import predict_band_gap

# Single composition
result = predict_band_gap("GaN")
# array([1.59])

# Multiple compositions
results = predict_band_gap(["NaCl", "TiO2", "GaN", "Si"])
# array([5.67, 2.13, 1.59, 0.72])
```

## Predictor API

```python
from smact.property_prediction import RoostPropertyPredictor

# Check available properties
print(RoostPropertyPredictor.available_properties)  # ['band_gap']

# Create predictor
predictor = RoostPropertyPredictor(property_name="band_gap", device="cpu")

# Predict
predictor.predict("NaCl")  # Single composition
predictor.predict(["NaCl", "TiO2", "GaN"])  # Multiple compositions

# With uncertainty quantification (requires robust model)
result = predictor.predict(["NaCl", "TiO2"], return_uncertainty=True)
print(result.predictions)  # [5.67, 2.13]
print(result.uncertainties)  # [0.42, 0.31]
print(result.unit)  # "eV"
```

## Key Features

1. **Clean API**: No model_path required — handled internally with defaults
2. **Class-level properties**: `RoostPropertyPredictor.available_properties`
3. **SMACT validation**: All compositions validated using `smact_validity()`
4. **Uncertainty quantification**: Aleatoric uncertainty via heteroscedastic (robust) loss
5. **Model registry**: Query available models, properties, and units
6. **Extensible design**: ABC-based architecture for adding new backends (Wren, CGCNN, etc.)
7. **Lazy imports**: torch/aviary only imported when actually used

## Architecture

The module uses an Abstract Base Class pattern (`BasePropertyPredictor` → `RoostPropertyPredictor`)
to support future model backends without API changes. All predictors share:

- `predict(compositions, return_uncertainty=False)` → `np.ndarray | PredictionResult`
- `PredictionResult` dataclass with `.predictions`, `.uncertainties`, `.unit`, `.to_dict()`

## Model Resolution

When loading a model by name, SMACT checks in order:

1. **Local directory path** — user-supplied explicit path
2. **`pretrained_models/`** — bundled in the repo (for source installs)
3. **`~/.cache/smact/models/`** — previously downloaded models
4. **Remote download** — fetches from GitHub Releases
