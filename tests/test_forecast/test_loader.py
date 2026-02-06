"""Tests for ForecastLoader orchestrator."""

from __future__ import annotations

import tempfile
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import xarray as xr

from ecmwf_downloader.bbox import BBox
from ecmwf_downloader.forecast import ForecastLoader


class TestForecastLoaderInit:
    """Test ForecastLoader initialization."""

    def test_init_with_tuple_bbox(self):
        """ForecastLoader should accept tuple bbox."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
            )
            assert isinstance(loader.bbox, BBox)
            assert loader.bbox.west == 7.7
            assert loader.bbox.south == 45.95

    def test_init_with_bbox_instance(self, sample_bbox: BBox):
        """ForecastLoader should accept BBox instance."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=sample_bbox,
                output_dir=tmpdir,
            )
            assert loader.bbox == sample_bbox

    def test_init_creates_output_dir(self):
        """ForecastLoader should create output directory if missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "new_dir" / "forecast"
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=str(output_path),
            )
            assert output_path.exists()

    def test_init_default_backend(self):
        """ForecastLoader should default to auto backend."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
            )
            assert loader.backend_name == "auto"

    def test_init_default_priority_strategy(self):
        """ForecastLoader should default to speed priority."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
            )
            assert loader.priority_strategy == "speed"

    def test_init_custom_priority_strategy(self):
        """ForecastLoader should accept reliability priority."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                priority_strategy="reliability",
            )
            assert loader.priority_strategy == "reliability"

    def test_init_default_pressure_levels(self):
        """ForecastLoader should have default pressure levels."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
            )
            assert loader.pressure_levels == [1000, 850, 700, 500, 300]

    def test_init_custom_pressure_levels(self):
        """ForecastLoader should accept custom pressure levels."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                pressure_levels=[850, 500],
            )
            assert loader.pressure_levels == [850, 500]

    def test_init_required_variables(self):
        """ForecastLoader should store required_variables as set."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                required_variables=["strd", "q"],
            )
            assert loader.required_variables == {"strd", "q"}


class TestBackendSelection:
    """Test backend selection logic."""

    def test_explicit_backend_returns_that_backend(self):
        """Explicit backend selection should return that backend."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="ecmwf_opendata",
            )
            assert loader._select_backend() == "ecmwf_opendata"

    def test_explicit_openmeteo_backend(self):
        """Explicit openmeteo backend selection."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="openmeteo_ifs",
            )
            assert loader._select_backend() == "openmeteo_ifs"

    def test_auto_speed_priority_selects_openmeteo(self):
        """Auto mode with speed priority should select openmeteo."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="auto",
                priority_strategy="speed",
            )
            assert loader._select_backend() == "openmeteo_ifs"

    def test_auto_reliability_priority_selects_ecmwf(self):
        """Auto mode with reliability priority should select ecmwf."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="auto",
                priority_strategy="reliability",
            )
            assert loader._select_backend() == "ecmwf_opendata"

    def test_blacklisted_var_strd_forces_ecmwf(self):
        """Required strd should force ecmwf backend."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="auto",
                priority_strategy="speed",  # Would normally select openmeteo
                required_variables=["strd"],
            )
            assert loader._select_backend() == "ecmwf_opendata"

    def test_blacklisted_var_q_forces_ecmwf(self):
        """Required q should force ecmwf backend."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="auto",
                priority_strategy="speed",
                required_variables=["q"],
            )
            assert loader._select_backend() == "ecmwf_opendata"

    def test_multiple_blacklisted_vars_force_ecmwf(self):
        """Multiple blacklisted variables should force ecmwf."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="auto",
                required_variables=["strd", "q", "t2m"],
            )
            assert loader._select_backend() == "ecmwf_opendata"

    def test_non_blacklisted_vars_allow_openmeteo(self):
        """Non-blacklisted required vars should allow openmeteo with speed priority."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="auto",
                priority_strategy="speed",
                required_variables=["t2m", "tp", "ssrd"],
            )
            assert loader._select_backend() == "openmeteo_ifs"


class TestGetBackend:
    """Test backend instantiation."""

    def test_get_backend_caches_instance(self):
        """_get_backend should cache the backend instance."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="openmeteo_ifs",
            )
            backend1 = loader._get_backend()
            backend2 = loader._get_backend()
            assert backend1 is backend2

    def test_get_backend_passes_pressure_levels(self):
        """Backend should receive pressure levels from loader."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="openmeteo_ifs",
                pressure_levels=[850, 500],
            )
            backend = loader._get_backend()
            assert backend.pressure_levels == [500, 850]  # Sorted

    def test_get_backend_passes_bbox(self, sample_bbox: BBox):
        """Backend should receive bbox from loader."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=sample_bbox,
                output_dir=tmpdir,
                backend="openmeteo_ifs",
            )
            backend = loader._get_backend()
            assert backend.bbox == sample_bbox


class TestOpen:
    """Test opening forecast data."""

    def test_open_raises_if_no_data(self):
        """open() should raise if no forecast data exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
            )
            with pytest.raises(FileNotFoundError, match="No forecast data"):
                loader.open()

    def test_open_returns_dataset_if_data_exists(
        self,
        sample_surface_dataset: xr.Dataset,
        sample_pressure_dataset: xr.Dataset,
    ):
        """open() should return dataset if forecast.zarr exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a mock zarr store
            zarr_path = Path(tmpdir) / "forecast.zarr"
            # Rename z in surface to z_surf to avoid conflict with pressure z
            ds_surf = sample_surface_dataset.rename({"z": "z_surf"})
            merged = xr.merge([ds_surf, sample_pressure_dataset])
            merged.to_zarr(str(zarr_path))

            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
            )
            ds = loader.open()
            assert isinstance(ds, xr.Dataset)
            assert "t2m" in ds.data_vars


class TestClose:
    """Test resource cleanup."""

    def test_close_clears_backend(self):
        """close() should clear the cached backend."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
                backend="openmeteo_ifs",
            )
            # Get backend to cache it
            _ = loader._get_backend()
            assert loader._backend is not None

            loader.close()
            assert loader._backend is None

    def test_close_is_safe_without_backend(self):
        """close() should be safe even if backend never created."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
            )
            # Should not raise
            loader.close()


class TestGetExistingTimes:
    """Test detection of already-downloaded forecasts."""

    def test_empty_dir_returns_empty_set(self):
        """Empty directory should return empty set."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
            )
            assert loader._get_existing_times() == set()

    def test_nonexistent_zarr_returns_empty_set(self):
        """Missing zarr store should return empty set."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = ForecastLoader(
                bbox=(7.7, 45.95, 7.85, 46.05),
                output_dir=tmpdir,
            )
            assert loader._get_existing_times() == set()
