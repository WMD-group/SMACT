"""ROOST-based property predictor implementation.

This module provides the RoostPropertyPredictor class for predicting
material properties from chemical composition using pre-trained ROOST
(Representation Learning from Stoichiometry) models.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader

from smact.property_prediction.base_predictor import (
    BasePropertyPredictor,
    PredictionResult,
)
from smact.property_prediction.io import load_checkpoint
from smact.property_prediction.registry import (
    get_default_model,
    get_property_unit,
    get_supported_properties,
)

logger = logging.getLogger(__name__)


class RoostPropertyPredictor(BasePropertyPredictor):
    """Property predictor using pre-trained ROOST models.

    ROOST (Representation Learning from Stoichiometry) predicts material
    properties from composition alone using message passing on stoichiometric
    graphs. This class loads pre-trained checkpoints and performs inference.

    Example:
        >>> predictor = RoostPropertyPredictor(property_name="band_gap")
        >>> predictions = predictor.predict(["NaCl", "TiO2"])

        >>> # With fidelity selection
        >>> predictor = RoostPropertyPredictor(property_name="band_gap", fidelity="hse06")

        >>> # With uncertainty quantification
        >>> result = predictor.predict(["NaCl"], return_uncertainty=True)
        >>> print(result.uncertainties)

    Attributes:
        available_properties: Class attribute listing all supported properties.
        elem_embedding: Element embedding type used by the model.
        batch_size: Batch size for inference.
    """

    # Available properties - can be accessed as class attribute
    available_properties: ClassVar[list[str]] = get_supported_properties()

    @property
    def supported_properties(self) -> list[str]:
        """List of properties supported by ROOST predictor."""
        return self.available_properties

    def __init__(
        self,
        property_name: str,
        fidelity: str | None = None,
        model_name: str | None = None,
        model_path: str | Path | None = None,
        device: str = "cpu",
        elem_embedding: str = "matscholar200",
        batch_size: int = 128,
        **kwargs: Any,
    ):
        """Initialise the ROOST property predictor.

        Args:
            property_name: Name of the property to predict
                (e.g., "band_gap", "bulk_modulus").
            fidelity: Fidelity level (e.g., "pbe", "hse06" for band_gap).
                If None, uses the default fidelity for the property.
            model_name: Specific model version to use
                (e.g., "Roost-MP-2024.1.0-band_gap-pbe").
                If None, uses default for property/fidelity.
            model_path: Local path to model directory. Overrides model_name.
            device: Device to run the model on ("cpu" or "cuda").
            elem_embedding: Element embedding type (default: "matscholar200").
            batch_size: Batch size for inference.
            **kwargs: Additional arguments.
        """
        super().__init__(
            property_name=property_name,
            fidelity=fidelity,
            model_name=model_name,
            model_path=model_path,
            device=device,
            **kwargs,
        )

        self.elem_embedding = elem_embedding
        self.batch_size = batch_size
        self._normaliser: Any = None
        self._robust = False
        self._target_name: str = property_name

        # Resolve model to load
        if model_path is None and model_name is None:
            self.model_name = get_default_model(property_name, fidelity)

        self._load_model()

    def _load_model(self) -> None:
        """Load the pre-trained ROOST model from checkpoint."""
        from aviary.data import Normalizer
        from aviary.roost.model import Roost

        # Load checkpoint
        model_source = self.model_path if self.model_path else self.model_name
        if model_source is None:
            raise ValueError("No model specified. Provide either model_path or model_name.")

        checkpoint = load_checkpoint(model_source, device=self.device)

        model_params = checkpoint["model_params"]
        self._metadata = checkpoint.get("metadata", {})

        # Track if model is robust (has uncertainty estimation)
        self._robust = model_params.get("robust", False)

        # Get target name from task_dict
        task_dict = model_params.get("task_dict", {})
        if task_dict:
            self._target_name = next(iter(task_dict.keys()))

        # Reconstruct the model
        self.model = Roost(**model_params)
        self.model.load_state_dict(checkpoint["state_dict"])
        self.model.to(self.device)
        self.model.eval()

        # Load normaliser for denormalisation
        normaliser_dict = checkpoint.get("normalizer_dict", {})
        normaliser_state = normaliser_dict.get(self._target_name)
        if normaliser_state is not None:
            self._normaliser = Normalizer.from_state_dict(normaliser_state)

        logger.info(f"Loaded ROOST model for {self.property_name} (robust={self._robust}, device={self.device})")

    def predict(
        self,
        compositions: str | list[str],
        return_uncertainty: bool = False,
        validate_smact: bool = True,
    ) -> np.ndarray | PredictionResult:
        """Predict property values for given compositions.

        Args:
            compositions: Single composition string or list of composition
                strings (e.g., "NaCl", ["TiO2", "GaN"]).
            return_uncertainty: If True, return PredictionResult with
                uncertainty estimates (requires robust model).
            validate_smact: Whether to validate compositions using SMACT.
                Set to False to skip validation.

        Returns:
            Property predictions as numpy array, or PredictionResult
            object if return_uncertainty is True.
        """
        from aviary.roost.data import CompositionData, collate_batch

        compositions = self._validate_compositions(compositions, validate_smact)

        # Create dataframe for CompositionData
        comp_data = pd.DataFrame(
            {
                "material_id": list(range(len(compositions))),
                "composition": compositions,
                self._target_name: [0.0] * len(compositions),  # Dummy target
            }
        )

        # Create dataset and dataloader
        task_dict = {self._target_name: "regression"}
        dataset = CompositionData(
            df=comp_data,
            task_dict=task_dict,
            inputs="composition",
            identifiers=["material_id", "composition"],
        )

        data_loader = DataLoader(
            dataset,
            batch_size=self.batch_size,
            shuffle=False,
            collate_fn=collate_batch,  # type: ignore[arg-type]  # aviary annotation is too narrow
        )

        # Run inference
        all_preds: list[np.ndarray] = []
        all_uncertainties: list[np.ndarray] = []

        self.model.eval()
        with torch.no_grad():
            for batch in data_loader:
                inputs = batch[0]  # First element is inputs tuple

                # Move inputs to device
                inputs = tuple(t.to(self.device) if isinstance(t, torch.Tensor) else t for t in inputs)

                # Forward pass
                outputs = self.model(*inputs)
                output = outputs[0]  # First (and only) target

                if self._robust:
                    # Robust model outputs mean and log_std
                    mean, log_std = output.unbind(dim=1)
                    preds = mean.cpu()
                    uncertainties = torch.exp(log_std).cpu()
                else:
                    preds = output.squeeze(1).cpu()
                    uncertainties = None

                # Denormalise predictions
                if self._normaliser is not None:
                    preds = self._normaliser.denorm(preds)
                    if uncertainties is not None:
                        # Scale uncertainty by normaliser std
                        uncertainties = uncertainties * self._normaliser.std

                all_preds.append(preds.numpy())
                if uncertainties is not None:
                    all_uncertainties.append(uncertainties.numpy())

        predictions = np.concatenate(all_preds)

        if return_uncertainty:
            uncertainties_arr = None
            if all_uncertainties:
                uncertainties_arr = np.concatenate(all_uncertainties)

            return PredictionResult(
                predictions=predictions,
                uncertainties=uncertainties_arr,
                unit=get_property_unit(self.property_name),
                compositions=compositions,
            )

        return predictions

    @classmethod
    def from_checkpoint(
        cls,
        checkpoint_path: str | Path,
        property_name: str,
        device: str = "cpu",
        **kwargs: Any,
    ) -> RoostPropertyPredictor:
        """Create predictor from a local checkpoint directory.

        Args:
            checkpoint_path: Path to checkpoint directory containing
                model.json, model.pt, and state.pt files.
            property_name: Name of the property.
            device: Device to run on ("cpu" or "cuda").
            **kwargs: Additional arguments passed to constructor.

        Returns:
            Initialised RoostPropertyPredictor.
        """
        return cls(
            property_name=property_name,
            model_path=checkpoint_path,
            device=device,
            **kwargs,
        )

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"{type(self).__name__}("
            f"property='{self.property_name}', "
            f"fidelity='{self.fidelity}', "
            f"model='{self.model_name}', "
            f"robust={self._robust})"
        )
