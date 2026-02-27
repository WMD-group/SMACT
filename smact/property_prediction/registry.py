"""Model registry for discovering and resolving pretrained models."""

from __future__ import annotations

import json
import logging
from pathlib import Path

from smact.property_prediction.config import (
    DEFAULT_MODELS,
    MODELS_CACHE,
    PRETRAINED_MODELS_BASE_URL,
    PRETRAINED_MODELS_DIR,
    PROPERTY_METADATA,
)

logger = logging.getLogger(__name__)


def _is_valid_model_dir(path: Path) -> bool:
    """Check if a directory contains a complete model."""
    return (
        path.is_dir()
        and (path / "model.json").exists()
        and (path / "model.pt").exists()
        and (path / "state.pt").exists()
    )


def get_available_models(include_cached: bool = True) -> list[str]:
    """Query for available pretrained models.

    Attempts to fetch the model manifest from the remote server.
    Falls back to listing locally cached models if remote is unavailable.

    Args:
        include_cached: Whether to include locally cached models in the list.

    Returns:
        List of available model names.
    """
    import requests

    models = set()

    # Try to fetch remote manifest
    try:
        manifest_url = f"{PRETRAINED_MODELS_BASE_URL}manifest.json"
        response = requests.get(manifest_url, timeout=10)
        if response.status_code == 200:
            manifest = json.loads(response.content)
            models.update(manifest.get("models", []))
    except (requests.RequestException, json.JSONDecodeError) as e:
        logger.debug(f"Could not fetch remote manifest: {e}")

    # Include models from pretrained_models directory in the repository
    if PRETRAINED_MODELS_DIR.exists():
        for path in PRETRAINED_MODELS_DIR.iterdir():
            if _is_valid_model_dir(path):
                models.add(path.name)

    # Include locally cached models
    if include_cached and MODELS_CACHE.exists():
        for path in MODELS_CACHE.iterdir():
            if _is_valid_model_dir(path):
                models.add(path.name)

    return sorted(models)


def get_default_model(
    property_name: str,
    fidelity: str | None = None,
) -> str:
    """Get the default model name for a property/fidelity combination.

    Args:
        property_name: Name of the property (e.g., "band_gap", "bulk_modulus").
        fidelity: Optional fidelity level (e.g., "pbe", "hse06").
            If None, uses the default fidelity for that property.

    Returns:
        Model name string (e.g., "Roost-MP-2024.1.0-band_gap-pbe").

    Raises:
        ValueError: If property or fidelity is not supported.
    """
    if property_name not in DEFAULT_MODELS:
        available = list(DEFAULT_MODELS.keys())
        raise ValueError(f"Unknown property: '{property_name}'. Available properties: {available}")

    prop_models = DEFAULT_MODELS[property_name]
    available_fidelities = [k for k in prop_models if k != "default"]

    # Use default fidelity if not specified
    if fidelity is None:
        default_value = prop_models.get("default")
        if default_value is None:
            raise ValueError(f"No default model configured for property '{property_name}'")

        # Check if default points to a fidelity key or directly to a model name
        if default_value in prop_models:
            # default is a fidelity key (e.g., "pbe" for band_gap)
            fidelity = default_value
        else:
            # default is the model name directly (properties without fidelity variants)
            return default_value

    if fidelity not in prop_models:
        raise ValueError(
            f"Fidelity '{fidelity}' not available for property '{property_name}'. "
            f"Available fidelities: {available_fidelities}"
        )

    return prop_models[fidelity]


def get_supported_properties() -> list[str]:
    """Get list of properties with available models.

    Returns:
        List of property names that have pretrained models.
    """
    return list(DEFAULT_MODELS.keys())


def get_property_fidelities(property_name: str) -> list[str] | None:
    """Get available fidelity levels for a property.

    Args:
        property_name: Name of the property.

    Returns:
        List of available fidelities, or None if property has no fidelity variants.

    Raises:
        ValueError: If property is not supported.
    """
    if property_name not in PROPERTY_METADATA:
        available = list(PROPERTY_METADATA.keys())
        raise ValueError(f"Unknown property: '{property_name}'. Available properties: {available}")

    fidelities = PROPERTY_METADATA[property_name].get("fidelities")
    if fidelities is None:
        return None
    return list(fidelities)


def get_property_unit(property_name: str) -> str:
    """Get the unit string for a property.

    Args:
        property_name: Name of the property.

    Returns:
        Unit string (e.g., "eV", "GPa").
    """
    unit = PROPERTY_METADATA.get(property_name, {}).get("unit", "")
    return str(unit) if unit else ""


def get_property_description(property_name: str) -> str:
    """Get the description for a property.

    Args:
        property_name: Name of the property.

    Returns:
        Description string.
    """
    desc = PROPERTY_METADATA.get(property_name, {}).get("description", "")
    return str(desc) if desc else ""


def model_exists(model_name: str) -> bool:
    """Check if a model exists locally or remotely.

    Args:
        model_name: Name of the model to check.

    Returns:
        True if model exists, False otherwise.
    """
    # Check pretrained_models directory in the repository
    pretrained_path = PRETRAINED_MODELS_DIR / model_name
    if _is_valid_model_dir(pretrained_path):
        return True

    # Check local cache
    cached_path = MODELS_CACHE / model_name
    if _is_valid_model_dir(cached_path):
        return True

    # Check remote availability
    import requests

    try:
        url = f"{PRETRAINED_MODELS_BASE_URL}{model_name}.tar.gz"
        response = requests.head(url, timeout=5, allow_redirects=True)
        return response.status_code == 200
    except requests.RequestException:
        return False


def parse_model_name(model_name: str) -> dict[str, str | None]:
    """Parse a model name into its components.

    Model naming convention: Roost-<dataset>-<version>-<property>[-<fidelity>]

    Args:
        model_name: Full model name string.

    Returns:
        Dictionary with keys: model_type, dataset, version, property, fidelity.
    """
    parts = model_name.split("-")

    if len(parts) < 4:
        return {
            "model_type": None,
            "dataset": None,
            "version": None,
            "property": None,
            "fidelity": None,
        }

    result = {
        "model_type": parts[0],  # e.g., "Roost"
        "dataset": parts[1],  # e.g., "MP"
        "version": parts[2],  # e.g., "2024.1.0"
        "property": parts[3],  # e.g., "band_gap"
        "fidelity": parts[4] if len(parts) > 4 else None,  # e.g., "pbe"
    }

    # Handle properties with underscores (e.g., "band_gap")
    # If we have more parts, they might be part of property name
    if len(parts) > 5:
        # Rejoin property parts
        result["property"] = "_".join(parts[3:-1])
        result["fidelity"] = parts[-1]
    elif len(parts) == 5:
        # Could be property_fidelity or property_part
        # Check if last part is a known fidelity
        all_fidelities: set[str] = set()
        for prop_meta in PROPERTY_METADATA.values():
            fidelities = prop_meta.get("fidelities")
            if fidelities and isinstance(fidelities, list):
                all_fidelities.update(fidelities)

        if parts[4] in all_fidelities:
            result["fidelity"] = parts[4]
        else:
            # It's part of the property name
            result["property"] = f"{parts[3]}_{parts[4]}"
            result["fidelity"] = None

    return result
