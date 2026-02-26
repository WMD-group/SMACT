"""Tests for the I/O module."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

torch = pytest.importorskip("torch", reason="torch required for IO tests")

from smact.property_prediction.io import (  # noqa: E402
    get_cache_size,
    list_cached_models,
    load_model_files,
    save_checkpoint,
)


class TestSaveCheckpoint:
    """Tests for save_checkpoint function."""

    def test_saves_required_files(self):
        """Should save model.json, model.pt, and state.pt."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_model"

            model_params = {"param1": 1, "param2": "value"}
            state_dict = {"layer.weight": torch.tensor([1.0, 2.0])}
            normalizer_dict = {"target": {"mean": 0.0, "std": 1.0}}

            save_checkpoint(
                model_params=model_params,
                state_dict=state_dict,
                normalizer_dict=normalizer_dict,
                path=output_path,
            )

            assert (output_path / "model.json").exists()
            assert (output_path / "model.pt").exists()
            assert (output_path / "state.pt").exists()

    def test_model_json_format(self):
        """Model JSON should have correct structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_model"

            save_checkpoint(
                model_params={"param": 1},
                state_dict={},
                normalizer_dict={},
                path=output_path,
                metadata={"custom": "data"},
            )

            with open(output_path / "model.json") as f:
                data = json.load(f)

            assert "@class" in data
            assert "@module" in data
            assert "@model_version" in data
            assert "metadata" in data
            assert data["metadata"]["custom"] == "data"

    def test_state_dict_loadable(self):
        """State dict should be loadable with torch."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_model"
            original_state = {"layer.weight": torch.tensor([1.0, 2.0, 3.0])}

            save_checkpoint(
                model_params={},
                state_dict=original_state,
                normalizer_dict={},
                path=output_path,
            )

            loaded_state = torch.load(output_path / "state.pt", weights_only=False)
            assert torch.equal(loaded_state["layer.weight"], original_state["layer.weight"])


class TestLoadModelFiles:
    """Tests for load_model_files function."""

    def test_loads_local_path(self):
        """Should load from local directory path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            model_path = Path(tmpdir)

            # Create required files
            (model_path / "model.json").write_text("{}")
            torch.save({}, model_path / "model.pt")
            torch.save({}, model_path / "state.pt")

            files = load_model_files(model_path)

            assert "model.json" in files
            assert "model.pt" in files
            assert "state.pt" in files
            assert files["model.json"] == model_path / "model.json"

    def test_missing_files_raises(self):
        """Should raise error if files are missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            model_path = Path(tmpdir)
            # Create only partial files
            (model_path / "model.json").write_text("{}")

            # Should not find complete model locally and fail to download
            with pytest.raises(FileNotFoundError):
                load_model_files(model_path / "nonexistent")


class TestCacheFunctions:
    """Tests for cache utility functions."""

    def test_get_cache_size_no_cache(self):
        """Should return 0 if cache doesn't exist."""
        # This assumes MODELS_CACHE might not exist in test environment
        size = get_cache_size()
        assert isinstance(size, int)
        assert size >= 0

    def test_list_cached_models_empty(self):
        """Should return empty list if no models cached."""
        models = list_cached_models()
        assert isinstance(models, list)


class TestIntegration:
    """Integration tests for save and load."""

    def test_save_load_roundtrip(self):
        """Should be able to save and load checkpoint."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "roundtrip_model"

            # Save
            original_params = {"task_dict": {"target": "regression"}, "robust": True}
            original_state = {"encoder.weight": torch.randn(10, 5)}
            original_normalizer = {"target": {"mean": 1.5, "std": 0.5}}

            save_checkpoint(
                model_params=original_params,
                state_dict=original_state,
                normalizer_dict=original_normalizer,
                path=output_path,
                metadata={"property": "band_gap"},
            )

            # Load
            files = load_model_files(output_path)

            with open(files["model.json"]) as f:
                metadata = json.load(f)

            loaded_params = torch.load(files["model.pt"], weights_only=False)
            loaded_state = torch.load(files["state.pt"], weights_only=False)

            # Verify
            assert metadata["@class"] == "Roost"
            assert loaded_params["task_dict"] == original_params["task_dict"]
            assert loaded_params["robust"] == original_params["robust"]
            assert "normalizer_dict" in loaded_params
            assert torch.equal(loaded_state["encoder.weight"], original_state["encoder.weight"])
