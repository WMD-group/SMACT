"""Configuration constants for the property prediction module."""

from __future__ import annotations

import os
from pathlib import Path

# Pretrained models bundled in the repository (available when running from source)
PRETRAINED_MODELS_DIR = Path(__file__).resolve().parent.parent.parent / "pretrained_models"

# Cache locations (following MatGL pattern)
SMACT_CACHE = Path(os.path.expanduser("~")) / ".cache" / "smact"
MODELS_CACHE = SMACT_CACHE / "models"

# Remote URL base for pretrained models
# Models will be hosted as GitHub releases
PRETRAINED_MODELS_BASE_URL = "https://github.com/WMD-group/SMACT/releases/download/models-v1/"

# Model naming convention: Roost-<dataset>-<version>-<property>[-<fidelity>]
# Examples:
#   Roost-MP-2024.1.0-band_gap-pbe
#   Roost-MP-2024.1.0-band_gap-hse06
#   Roost-MP-2024.1.0-bulk_modulus

# Default model versions for each property+fidelity combination
DEFAULT_MODELS: dict[str, dict[str, str]] = {
    "band_gap": {
        # DFT band gap model trained on 103K MP samples (MAE: 0.28 eV)
        "default": "Roost-MP-2024.12.0-band_gap",
    },
}

# Property metadata including units and descriptions
PROPERTY_METADATA: dict[str, dict[str, str | list[str] | None]] = {
    "band_gap": {
        "unit": "eV",
        "description": "Electronic band gap",
        "fidelities": None,
    },
}

# Element embeddings available in aviary
AVAILABLE_EMBEDDINGS = [
    "matscholar200",  # Default, 200-dim from MatScholar
    "cgcnn92",  # 92-dim from CGCNN
    "megnet16",  # 16-dim from MEGNet
    "onehot112",  # 112-dim one-hot encoding
]

# Default model hyperparameters for training
DEFAULT_MODEL_PARAMS = {
    "elem_embedding": "matscholar200",
    "elem_fea_len": 64,
    "n_graph": 3,
    "elem_heads": 3,
    "elem_gate": [256],
    "elem_msg": [256],
    "cry_heads": 3,
    "cry_gate": [256],
    "cry_msg": [256],
    "trunk_hidden": [1024, 512],
    "out_hidden": [256, 128, 64],
}

# Default training parameters
DEFAULT_TRAINING_PARAMS = {
    "epochs": 100,
    "batch_size": 128,
    "learning_rate": 3e-4,
    "weight_decay": 1e-6,
    "patience": None,  # Early stopping patience (None = no early stopping)
    "robust": True,  # Use heteroscedastic loss for uncertainty
}
