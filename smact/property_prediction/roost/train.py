"""Training utilities for ROOST models.

This module provides functions for training ROOST models from scratch
using composition-property data.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import pandas as pd
import torch
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, Subset

from smact.property_prediction.config import (
    DEFAULT_MODEL_PARAMS,
    DEFAULT_TRAINING_PARAMS,
    PROPERTY_METADATA,
)
from smact.property_prediction.io import save_checkpoint

logger = logging.getLogger(__name__)


def train_roost_model(
    data_path: str | Path,
    target_column: str,
    property_name: str,
    output_dir: str | Path,
    fidelity: str | None = None,
    # Model hyperparameters
    elem_embedding: str | None = None,
    elem_fea_len: int | None = None,
    n_graph: int | None = None,
    robust: bool | None = None,
    # Training parameters
    epochs: int | None = None,
    batch_size: int | None = None,
    learning_rate: float | None = None,
    weight_decay: float | None = None,
    patience: int | None = None,
    # Data splits
    test_size: float = 0.1,
    val_size: float = 0.1,
    data_seed: int = 42,
    # Misc
    device: str | None = None,
    ensemble_size: int = 1,
) -> Path:
    """Train a ROOST model for property prediction.

    Args:
        data_path: Path to CSV with composition and target columns.
        target_column: Column name for target property in the CSV.
        property_name: Property name for model naming and metadata.
        output_dir: Directory to save trained model.
        fidelity: Optional fidelity level for model naming.
        elem_embedding: Element embedding type. Defaults to "matscholar200".
        elem_fea_len: Element feature length. Defaults to 64.
        n_graph: Number of message passing layers. Defaults to 3.
        robust: Use heteroscedastic loss for uncertainty. Defaults to True.
        epochs: Number of training epochs. Defaults to 100.
        batch_size: Training batch size. Defaults to 128.
        learning_rate: Initial learning rate. Defaults to 3e-4.
        weight_decay: Weight decay for optimiser. Defaults to 1e-6.
        patience: Early stopping patience (None = no early stopping).
        test_size: Fraction of data for testing.
        val_size: Fraction of data for validation.
        data_seed: Random seed for data splitting.
        device: Device to train on (auto-detect if None).
        ensemble_size: Number of models to train in ensemble.

    Returns:
        Path to saved model directory.

    Example:
        >>> train_roost_model(
        ...     data_path="band_gaps.csv",
        ...     target_column="band_gap",
        ...     property_name="band_gap",
        ...     output_dir="./models",
        ...     fidelity="pbe",
        ...     epochs=100,
        ...     robust=True,
        ... )
    """
    from aviary.data import Normalizer
    from aviary.roost.data import CompositionData, collate_batch
    from aviary.roost.model import Roost

    # Apply defaults
    elem_embedding = elem_embedding or DEFAULT_MODEL_PARAMS["elem_embedding"]
    elem_fea_len = elem_fea_len or DEFAULT_MODEL_PARAMS["elem_fea_len"]
    n_graph = n_graph or DEFAULT_MODEL_PARAMS["n_graph"]
    robust = robust if robust is not None else DEFAULT_TRAINING_PARAMS["robust"]
    epochs = epochs or DEFAULT_TRAINING_PARAMS["epochs"]
    batch_size = batch_size or DEFAULT_TRAINING_PARAMS["batch_size"]
    learning_rate = learning_rate or DEFAULT_TRAINING_PARAMS["learning_rate"]
    weight_decay = weight_decay or DEFAULT_TRAINING_PARAMS["weight_decay"]
    patience = patience if patience is not None else DEFAULT_TRAINING_PARAMS["patience"]

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    logger.info(f"Training ROOST model on {device}")

    # Load data
    train_data = pd.read_csv(data_path, keep_default_na=False, na_values=[])

    if target_column not in train_data.columns:
        raise ValueError(f"Target column '{target_column}' not found in data")

    if "composition" not in train_data.columns:
        raise ValueError("Data must have 'composition' column")

    # Add material_id if not present
    if "material_id" not in train_data.columns:
        train_data["material_id"] = train_data.index.astype(str)

    logger.info(f"Loaded {len(train_data)} samples from {data_path}")

    # Setup task
    task_dict = {target_column: "regression"}

    # Create dataset
    dataset = CompositionData(
        df=train_data,
        task_dict=task_dict,
        inputs="composition",
        identifiers=["material_id", "composition"],
    )
    n_targets = dataset.n_targets

    # Split data
    indices = list(range(len(dataset)))
    train_idx, test_idx = train_test_split(indices, test_size=test_size, random_state=data_seed)
    train_idx, val_idx = train_test_split(train_idx, test_size=val_size / (1 - test_size), random_state=data_seed)

    train_set = Subset(dataset, train_idx)
    val_set = Subset(dataset, val_idx)

    logger.info(f"Data split: train={len(train_set)}, val={len(val_set)}, test={len(test_idx)}")

    # Create dataloaders
    data_params = {
        "batch_size": batch_size,
        "num_workers": 0,
        "pin_memory": False,
        "collate_fn": collate_batch,
    }

    train_loader = DataLoader(train_set, shuffle=True, **data_params)
    val_loader = DataLoader(val_set, shuffle=False, **data_params)

    # Model parameters
    model_params = {
        "task_dict": task_dict,
        "robust": robust,
        "n_targets": n_targets,
        "elem_embedding": elem_embedding,
        "elem_fea_len": elem_fea_len,
        "n_graph": n_graph,
        "elem_heads": DEFAULT_MODEL_PARAMS["elem_heads"],
        "elem_gate": DEFAULT_MODEL_PARAMS["elem_gate"],
        "elem_msg": DEFAULT_MODEL_PARAMS["elem_msg"],
        "cry_heads": DEFAULT_MODEL_PARAMS["cry_heads"],
        "cry_gate": DEFAULT_MODEL_PARAMS["cry_gate"],
        "cry_msg": DEFAULT_MODEL_PARAMS["cry_msg"],
        "trunk_hidden": DEFAULT_MODEL_PARAMS["trunk_hidden"],
        "out_hidden": DEFAULT_MODEL_PARAMS["out_hidden"],
    }

    # Initialise model
    model = Roost(**model_params)
    model.to(device)

    # Setup optimiser
    optimiser = torch.optim.AdamW(
        model.parameters(),
        lr=learning_rate,
        weight_decay=weight_decay,
    )

    # Setup scheduler
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimiser,
        mode="min",
        factor=0.5,
        patience=10,
    )

    # Setup loss and normaliser
    if robust:
        from aviary.losses import robust_l1_loss

        criterion = robust_l1_loss
    else:
        criterion = torch.nn.L1Loss()

    # Calculate normaliser from training data
    train_targets = train_data.iloc[train_idx][target_column].to_numpy()
    normaliser = Normalizer()
    normaliser.fit(torch.tensor(train_targets, dtype=torch.float32))
    normaliser_dict = {target_column: normaliser.state_dict()}

    logger.info(f"Starting training for {epochs} epochs...")

    # Training loop
    best_val_loss = float("inf")
    best_state_dict = None
    epochs_without_improvement = 0

    for epoch in range(epochs):
        # Training
        model.train()
        train_loss = 0.0

        for batch in train_loader:
            inputs = batch[0]
            targets = batch[1]

            # Move to device
            inputs = tuple(t.to(device) if isinstance(t, torch.Tensor) else t for t in inputs)

            # Normalise target
            target = targets[0].to(device)
            target_norm = normaliser.norm(target)

            # Forward pass
            optimiser.zero_grad()
            outputs = model(*inputs)
            output = outputs[0]

            # Compute loss
            if robust:
                preds, log_std = output.unbind(dim=1)
                loss = criterion(preds, log_std, target_norm)
            else:
                loss = criterion(output.squeeze(1), target_norm)

            # Backward pass
            loss.backward()
            optimiser.step()

            train_loss += loss.item()

        train_loss /= len(train_loader)

        # Validation
        model.eval()
        val_loss = 0.0

        with torch.no_grad():
            for batch in val_loader:
                inputs = batch[0]
                targets = batch[1]

                inputs = tuple(t.to(device) if isinstance(t, torch.Tensor) else t for t in inputs)
                target = targets[0].to(device)
                target_norm = normaliser.norm(target)

                outputs = model(*inputs)
                output = outputs[0]

                if robust:
                    preds, log_std = output.unbind(dim=1)
                    loss = criterion(preds, log_std, target_norm)
                else:
                    loss = criterion(output.squeeze(1), target_norm)

                val_loss += loss.item()

        val_loss /= len(val_loader)

        # Update scheduler
        scheduler.step(val_loss)

        # Check for improvement
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state_dict = model.state_dict().copy()
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1

        # Early stopping
        if patience is not None and epochs_without_improvement >= patience:
            logger.info(f"Early stopping at epoch {epoch + 1}")
            break

        if (epoch + 1) % 10 == 0:
            logger.info(f"Epoch {epoch + 1}/{epochs}: train_loss={train_loss:.4f}, val_loss={val_loss:.4f}")

    # Load best model
    if best_state_dict is not None:
        model.load_state_dict(best_state_dict)

    # Generate model name and save
    model_name = f"Roost-Custom-{property_name}"
    if fidelity:
        model_name += f"-{fidelity}"

    output_path = Path(output_dir) / model_name

    # Build metadata
    metadata = {
        "property": property_name,
        "fidelity": fidelity,
        "model_type": "Roost",
        "training": {
            "data_path": str(data_path),
            "train_size": len(train_set),
            "val_size": len(val_set),
            "test_size": len(test_idx),
            "epochs": epochs,
            "batch_size": batch_size,
            "learning_rate": learning_rate,
            "best_val_loss": best_val_loss,
        },
    }

    # Save checkpoint
    save_checkpoint(
        model_params=model_params,
        state_dict=model.state_dict(),
        normalizer_dict=normaliser_dict,
        path=output_path,
        metadata=metadata,
    )

    # Write README
    _write_readme(output_path, property_name, fidelity, metadata)

    logger.info(f"Model saved to {output_path}")

    return output_path


def _write_readme(
    output_path: Path,
    property_name: str,
    fidelity: str | None,
    metadata: dict[str, Any],
) -> None:
    """Write a README file for the trained model."""
    unit = PROPERTY_METADATA.get(property_name, {}).get("unit", "")
    training = metadata.get("training", {})

    readme = f"""# {output_path.name}

## Property

- **Property**: {property_name}
- **Fidelity**: {fidelity or "N/A"}
- **Unit**: {unit}

## Training Details

- **Training samples**: {training.get("train_size", "N/A")}
- **Validation samples**: {training.get("val_size", "N/A")}
- **Test samples**: {training.get("test_size", "N/A")}
- **Epochs**: {training.get("epochs", "N/A")}
- **Best validation loss**: {training.get("best_val_loss", "N/A"):.4f}

## Usage

```python
from smact.property_prediction import RoostPropertyPredictor

predictor = RoostPropertyPredictor(
    property_name="{property_name}", model_path="{output_path}"
)
predictions = predictor.predict(["NaCl", "TiO2"])
```
"""
    with open(output_path / "README.md", "w") as f:
        f.write(readme)


def main() -> None:
    """Command-line interface for training ROOST models."""
    parser = argparse.ArgumentParser(description="Train a ROOST model for SMACT property prediction")

    parser.add_argument(
        "--data-path",
        required=True,
        help="Path to training CSV with composition and target columns",
    )
    parser.add_argument(
        "--target",
        required=True,
        help="Target column name in the CSV",
    )
    parser.add_argument(
        "--property",
        required=True,
        help="Property name (e.g., band_gap, bulk_modulus)",
    )
    parser.add_argument(
        "--output-dir",
        default="./models",
        help="Output directory for trained model",
    )
    parser.add_argument(
        "--fidelity",
        default=None,
        help="Fidelity level (e.g., pbe, hse06)",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=100,
        help="Number of training epochs",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=128,
        help="Training batch size",
    )
    parser.add_argument(
        "--learning-rate",
        type=float,
        default=3e-4,
        help="Initial learning rate",
    )
    parser.add_argument(
        "--robust",
        action="store_true",
        help="Use robust (uncertainty-aware) training",
    )
    parser.add_argument(
        "--device",
        default=None,
        help="Device to train on (cpu or cuda)",
    )
    parser.add_argument(
        "--patience",
        type=int,
        default=None,
        help="Early stopping patience",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for data splitting",
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    train_roost_model(
        data_path=args.data_path,
        target_column=args.target,
        property_name=args.property,
        output_dir=args.output_dir,
        fidelity=args.fidelity,
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.learning_rate,
        robust=args.robust,
        device=args.device,
        patience=args.patience,
        data_seed=args.seed,
    )


if __name__ == "__main__":
    main()
