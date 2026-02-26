# SMACT Pretrained Models

This directory contains pretrained models distributed with the SMACT repository, following the [MatGL](https://github.com/materialsvirtuallab/matgl) convention of per-model directories with README documentation.

## Available Models

| Model                                                       | Property      | Architecture | Dataset | Test MAE |
| ----------------------------------------------------------- | ------------- | ------------ | ------- | -------- |
| [Roost-MP-2024.12.0-band_gap](Roost-MP-2024.12.0-band_gap/) | Band gap (eV) | ROOST        | MP 103K | 0.283 eV |

## Model Resolution Order

When loading a model by name, SMACT checks the following locations in order:

1. **Local directory path** -- if you pass an explicit path
2. **`pretrained_models/`** -- this directory (available when running from the repo)
3. **Cache `~/.cache/smact/models/`** -- previously downloaded models
4. **Remote download** -- fetched from GitHub releases

For pip installs (where `pretrained_models/` is not bundled), models are downloaded on first use and cached locally.

## Adding a New Model

1. Create a subdirectory named `<Architecture>-<Dataset>-<Version>-<Property>` (e.g., `Roost-MP-2025.1.0-bulk_modulus`)
2. Place the three checkpoint files: `model.json`, `model.pt`, `state.pt`
3. Add a `README.md` following the existing per-model template (Aim / Training Dataset / Model Details / Performance / References)
4. Update the table above
5. Register the model in `smact/property_prediction/config.py` under `DEFAULT_MODELS`
