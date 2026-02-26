"""Model I/O utilities for loading, saving, and caching pretrained models."""

from __future__ import annotations

import json
import logging
import os
import shutil
import tarfile
from pathlib import Path
from typing import Any

import requests
import torch

from smact.property_prediction.config import MODELS_CACHE, PRETRAINED_MODELS_BASE_URL, PRETRAINED_MODELS_DIR

logger = logging.getLogger(__name__)


class RemoteFile:
    """Handle download and caching of remote model files.

    Downloads model archives from a remote URL and extracts them to
    the local cache directory.

    Attributes:
        uri: The remote URL of the model archive.
        cache_location: Local directory for caching models.
        model_name: Name of the model (derived from URI).
        local_path: Path to the extracted model directory.
    """

    def __init__(
        self,
        uri: str,
        cache_location: Path = MODELS_CACHE,
        force_download: bool = False,
    ):
        """Initialise the RemoteFile handler.

        Args:
            uri: Remote URL to download from.
            cache_location: Local directory for caching.
            force_download: If True, re-download even if cached.
        """
        self.uri = uri
        self.cache_location = cache_location
        self.force_download = force_download

        # Parse model name from URI
        self.model_name = uri.split("/")[-1].replace(".tar.gz", "")
        self.local_path = self.cache_location / self.model_name

        if not self.local_path.exists() or force_download:
            logger.info(f"Downloading model from {uri}...")
            self._download()
        else:
            logger.debug(f"Using cached model at {self.local_path}")

    def _download(self) -> None:
        """Download and extract the model archive."""
        os.makedirs(self.cache_location, exist_ok=True)

        # Download the archive
        response = requests.get(self.uri, stream=True, timeout=120)
        if response.status_code != 200:
            raise requests.RequestException(
                f"Failed to download model from {self.uri}. Status code: {response.status_code}"
            )

        # Save to temporary tar.gz file
        tar_path = self.cache_location / f"{self.model_name}.tar.gz"
        with open(tar_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        logger.info(f"Downloaded to {tar_path}, extracting...")

        # Extract the archive safely using data filter for security
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(self.cache_location, filter="data")

        # Clean up the tar file
        tar_path.unlink()
        logger.info(f"Model extracted to {self.local_path}")


def load_model_files(
    model_name: str | Path,
    force_download: bool = False,
) -> dict[str, Path]:
    """Load model files from local path or download from remote.

    Checks for model files in the following order:
    1. If model_name is a local directory path with required files
    2. If model exists in pretrained_models/ directory (repo installs)
    3. If model exists in cache (~/.cache/smact/models/)
    4. Download from remote URL

    Args:
        model_name: Model name or local path to model directory.
        force_download: If True, re-download even if cached.

    Returns:
        Dictionary mapping filename to Path for model.json, model.pt, state.pt.

    Raises:
        FileNotFoundError: If model files cannot be found or downloaded.
    """
    path = Path(model_name)
    required_files = ("model.json", "model.pt", "state.pt")

    # Check if it's a local directory path
    if path.is_dir() and all((path / fn).exists() for fn in required_files):
        logger.debug(f"Loading model from local path: {path}")
        return {fn: path / fn for fn in required_files}

    # Check pretrained_models directory in the repository
    pretrained_path = PRETRAINED_MODELS_DIR / str(model_name)
    if pretrained_path.is_dir() and all((pretrained_path / fn).exists() for fn in required_files):
        logger.debug(f"Loading model from pretrained_models: {pretrained_path}")
        return {fn: pretrained_path / fn for fn in required_files}

    # Check cache
    cached_path = MODELS_CACHE / str(model_name)
    if cached_path.is_dir() and all((cached_path / fn).exists() for fn in required_files) and not force_download:
        logger.debug(f"Loading model from cache: {cached_path}")
        return {fn: cached_path / fn for fn in required_files}

    # Download from remote
    try:
        remote = RemoteFile(
            f"{PRETRAINED_MODELS_BASE_URL}{model_name}.tar.gz",
            force_download=force_download,
        )
        return {fn: remote.local_path / fn for fn in required_files}
    except requests.RequestException as e:
        raise FileNotFoundError(
            f"Could not find or download model '{model_name}'. "
            f"Check that the model name is correct and you have internet access. "
            f"Original error: {e}"
        ) from e


def load_checkpoint(
    model_name: str | Path,
    device: str = "cpu",
    force_download: bool = False,
) -> dict[str, Any]:
    """Load a complete checkpoint from model files.

    Args:
        model_name: Model name or local path to model directory.
        device: Device to load tensors to ("cpu" or "cuda").
        force_download: If True, re-download even if cached.

    Returns:
        Checkpoint dictionary containing:
            - model_params: Model hyperparameters for reconstruction
            - state_dict: Model weights
            - normalizer_dict: Normaliser states for denormalisation
            - metadata: Additional model metadata
    """
    fpaths = load_model_files(model_name, force_download)

    map_location = torch.device(device)

    # Load metadata from JSON
    with open(fpaths["model.json"]) as f:
        metadata = json.load(f)

    # Load model parameters and state dict
    model_params = torch.load(fpaths["model.pt"], map_location=map_location, weights_only=False)
    state_dict = torch.load(fpaths["state.pt"], map_location=map_location, weights_only=False)

    # Check version compatibility (warn if mismatch)
    model_version = metadata.get("@model_version", 0)
    if model_version > 1:
        logger.warning(
            f"Model version ({model_version}) is newer than supported. Some features may not work correctly."
        )

    return {
        "model_params": model_params,
        "state_dict": state_dict,
        "metadata": metadata,
        "normalizer_dict": model_params.get("normalizer_dict", {}),
    }


def save_checkpoint(
    model_params: dict[str, Any],
    state_dict: dict[str, Any],
    normalizer_dict: dict[str, Any],
    path: Path | str,
    metadata: dict[str, Any] | None = None,
    model_class: str = "Roost",
    model_module: str = "aviary.roost.model",
) -> None:
    """Save a checkpoint in the standard SMACT format.

    Creates three files in the specified directory:
    - model.json: Metadata and hyperparameters
    - model.pt: Model parameters (including normaliser dict)
    - state.pt: Model weights (state dict)

    Args:
        model_params: Model hyperparameters for reconstruction.
        state_dict: Model weights.
        normalizer_dict: Normaliser states for each target.
        path: Directory to save the checkpoint to.
        metadata: Additional metadata to include.
        model_class: Class name for deserialisation.
        model_module: Module path for deserialisation.
    """
    path = Path(path)
    os.makedirs(path, exist_ok=True)

    # Save model parameters (including normaliser dict)
    params_to_save = {**model_params, "normalizer_dict": normalizer_dict}
    torch.save(params_to_save, path / "model.pt")

    # Save state dict
    torch.save(state_dict, path / "state.pt")

    # Save metadata JSON
    model_json = {
        "@class": model_class,
        "@module": model_module,
        "@model_version": 1,
        "metadata": metadata or {},
        "kwargs": model_params,
    }
    with open(path / "model.json", "w") as f:
        json.dump(model_json, f, indent=4, default=str)

    logger.info(f"Checkpoint saved to {path}")


def clear_cache(confirm: bool = True) -> None:
    """Clear the model cache directory.

    Args:
        confirm: If True, ask for confirmation before deleting.
    """
    if not MODELS_CACHE.exists():
        print(f"Cache directory {MODELS_CACHE} does not exist.")
        return

    if confirm:
        answer = input(f"Delete all cached models in {MODELS_CACHE}? (y/n): ").lower()
        if answer != "y":
            print("Cancelled.")
            return

    shutil.rmtree(MODELS_CACHE)
    print(f"Cleared cache at {MODELS_CACHE}")


def get_cache_size() -> int:
    """Get the total size of the model cache in bytes.

    Returns:
        Total size of cached models in bytes.
    """
    if not MODELS_CACHE.exists():
        return 0

    total_size = 0
    for path in MODELS_CACHE.rglob("*"):
        if path.is_file():
            total_size += path.stat().st_size

    return total_size


def list_cached_models() -> list[str]:
    """List all models currently in the cache.

    Returns:
        List of cached model names.
    """
    if not MODELS_CACHE.exists():
        return []

    models = []
    for path in MODELS_CACHE.iterdir():
        if path.is_dir() and (path / "model.json").exists():
            models.append(path.name)

    return sorted(models)
