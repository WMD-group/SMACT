"""Convert trained ROOST checkpoint to SMACT format.

This script converts checkpoints from the training script format to the
standardised SMACT model format (model.json, model.pt, state.pt).

Usage:
    python convert_checkpoint.py --input roost-mp-103k-v3-bugfix_best.pth --output Roost-MP-2024.12.0-band_gap
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import torch


def convert_checkpoint(
    input_path: Path,
    output_dir: Path,
    property_name: str = "band_gap",
    dataset: str = "MP",
    fidelity: str | None = None,
    description: str | None = None,
) -> None:
    """Convert training checkpoint to SMACT format.

    Args:
        input_path: Path to input .pth checkpoint.
        output_dir: Directory to save SMACT-format checkpoint.
        property_name: Name of the predicted property.
        dataset: Dataset the model was trained on.
        fidelity: Fidelity level (e.g., "pbe", "hse06") or None.
        description: Optional description of the model.
    """
    print(f"Loading checkpoint from {input_path}")
    checkpoint = torch.load(input_path, map_location="cpu", weights_only=False)

    # Extract components
    model_params = checkpoint.get("model_params", {})
    state_dict = checkpoint.get("model_state_dict", checkpoint.get("state_dict", {}))
    normalizer_mean = checkpoint.get("normalizer_mean")
    normalizer_std = checkpoint.get("normalizer_std")
    test_metrics = checkpoint.get("test_metrics", {})

    # Create normalizer dict with tensor values (required by aviary Normalizer)
    # Key must be the property name, not "target"
    normalizer_dict = {
        property_name: {
            "mean": torch.tensor(normalizer_mean) if normalizer_mean is not None else None,
            "std": torch.tensor(normalizer_std) if normalizer_std is not None else None,
        }
    }

    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save state dict
    print("Saving state.pt...")
    torch.save(state_dict, output_dir / "state.pt")

    # Build model params with required fields for Roost constructor
    # Convert lists to tuples where Roost expects Sequence[int]
    params_to_save = {
        "elem_embedding": model_params.get("elem_embedding", "matscholar200"),
        "elem_fea_len": model_params.get("elem_fea_len", 64),
        "n_graph": model_params.get("n_graph", 3),
        "elem_heads": model_params.get("elem_heads", 3),
        "elem_gate": tuple(model_params.get("elem_gate", [256])),
        "elem_msg": tuple(model_params.get("elem_msg", [256])),
        "cry_heads": model_params.get("cry_heads", 3),
        "cry_gate": tuple(model_params.get("cry_gate", [256])),
        "cry_msg": tuple(model_params.get("cry_msg", [256])),
        "trunk_hidden": tuple(model_params.get("trunk_hidden", [1024, 512])),
        "out_hidden": tuple(model_params.get("out_hidden", [256, 128, 64])),
        # Required by Roost constructor
        "robust": model_params.get("robust", True),
        "n_targets": [1],  # Single property prediction
        "task_dict": {property_name: "regression"},
        # Normalizer for denormalization during inference
        "normalizer_dict": normalizer_dict,
    }

    # Save model params
    print("Saving model.pt...")
    torch.save(params_to_save, output_dir / "model.pt")

    # Create metadata
    metadata = {
        "property": property_name,
        "dataset": dataset,
        "fidelity": fidelity,
        "unit": "eV" if property_name == "band_gap" else "",
        "training_date": datetime.now(tz=timezone.utc).isoformat(),
        "test_mae": test_metrics.get("mae"),
        "test_rmse": test_metrics.get("rmse"),
        "description": description or f"ROOST model for {property_name} prediction",
        "normalizer_mean": float(normalizer_mean) if normalizer_mean is not None else None,
        "normalizer_std": float(normalizer_std) if normalizer_std is not None else None,
    }

    # Create model.json
    model_json = {
        "@class": "Roost",
        "@module": "aviary.roost.model",
        "@model_version": 1,
        "metadata": metadata,
        "kwargs": model_params,
    }

    print("Saving model.json...")
    with open(output_dir / "model.json", "w") as f:
        json.dump(model_json, f, indent=4, default=str)

    print("\nCheckpoint converted successfully!")
    print(f"Output directory: {output_dir}")
    print("\nFiles created:")
    for f in output_dir.iterdir():
        size_kb = f.stat().st_size / 1024
        print(f"  {f.name}: {size_kb:.1f} KB")

    print("\nMetadata:")
    print(f"  Property: {property_name}")
    print(f"  Dataset: {dataset}")
    print(f"  Fidelity: {fidelity or 'default'}")
    if test_metrics:
        print(f"  Test MAE: {test_metrics.get('mae', 'N/A'):.4f}")
        print(f"  Test RMSE: {test_metrics.get('rmse', 'N/A'):.4f}")


def main() -> None:
    """CLI entry point for checkpoint conversion."""
    parser = argparse.ArgumentParser(description="Convert ROOST checkpoint to SMACT format")
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        required=True,
        help="Input checkpoint path (.pth file)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output directory name",
    )
    parser.add_argument(
        "--property",
        type=str,
        default="band_gap",
        help="Property name (default: band_gap)",
    )
    parser.add_argument(
        "--dataset",
        type=str,
        default="MP",
        help="Dataset name (default: MP)",
    )
    parser.add_argument(
        "--fidelity",
        type=str,
        default=None,
        help="Fidelity level (e.g., pbe, hse06)",
    )
    parser.add_argument(
        "--description",
        type=str,
        default=None,
        help="Model description",
    )

    args = parser.parse_args()

    convert_checkpoint(
        input_path=Path(args.input),
        output_dir=Path(args.output),
        property_name=args.property,
        dataset=args.dataset,
        fidelity=args.fidelity,
        description=args.description,
    )


if __name__ == "__main__":
    main()
