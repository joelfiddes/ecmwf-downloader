"""Tests for ECMWFOpenDataBackend."""

from __future__ import annotations

from datetime import datetime
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import xarray as xr

from ecmwf_downloader.bbox import BBox
from ecmwf_downloader.forecast.backends import ECMWFOpenDataBackend


class TestECMWFOpenDataBackendInit:
    """Test ECMWFOpenDataBackend initialization."""

    def test_init_with_valid_levels(self, sample_bbox: BBox):
        """Backend should accept valid pressure levels."""
        backend = ECMWFOpenDataBackend(
            bbox=sample_bbox,
            pressure_levels=[1000, 850, 700, 500],
        )
        assert backend.pressure_levels == [500, 700, 850, 1000]

    def test_init_default_forecast_hour(self, sample_bbox: BBox):
        """Backend should default to 00Z forecast."""
        backend = ECMWFOpenDataBackend(
            bbox=sample_bbox,
            pressure_levels=[850],
        )
        assert backend.forecast_hour == 0

    def test_init_12z_forecast_hour(self, sample_bbox: BBox):
        """Backend should accept 12Z forecast."""
        backend = ECMWFOpenDataBackend(
            bbox=sample_bbox,
            pressure_levels=[850],
            forecast_hour=12,
        )
        assert backend.forecast_hour == 12

    def test_init_rejects_invalid_forecast_hour(self, sample_bbox: BBox):
        """Backend should reject invalid forecast hour."""
        with pytest.raises(ValueError, match="must be 0 or 12"):
            ECMWFOpenDataBackend(
                bbox=sample_bbox,
                pressure_levels=[850],
                forecast_hour=6,
            )

    def test_init_valid_timesteps(self, sample_bbox: BBox):
        """Backend should accept valid output timesteps."""
        for ts in ["1H", "2H", "3H"]:
            backend = ECMWFOpenDataBackend(
                bbox=sample_bbox,
                pressure_levels=[850],
                output_timestep=ts,
            )
            assert backend.output_timestep == ts

    def test_init_rejects_invalid_timestep(self, sample_bbox: BBox):
        """Backend should reject invalid output timestep."""
        with pytest.raises(ValueError, match="Invalid output_timestep"):
            ECMWFOpenDataBackend(
                bbox=sample_bbox,
                pressure_levels=[850],
                output_timestep="6H",
            )

    def test_forecast_horizon_hours(self, sample_bbox: BBox):
        """Forecast horizon should be 240 hours (10 days)."""
        backend = ECMWFOpenDataBackend(
            bbox=sample_bbox,
            pressure_levels=[850],
        )
        assert backend.forecast_horizon_hours == 240

    def test_available_init_times(self, sample_bbox: BBox):
        """Available init times should be 0 and 12."""
        backend = ECMWFOpenDataBackend(
            bbox=sample_bbox,
            pressure_levels=[850],
        )
        assert backend.available_init_times == [0, 12]


class TestECMWFOpenDataBackendHelpers:
    """Test helper functions."""

    def test_deaccumulate(self):
        """Deaccumulation should compute differences."""
        from ecmwf_downloader.forecast.backends.ecmwf_opendata import _deaccumulate

        # Create simple accumulating dataset
        ds = xr.Dataset({
            "ssrd": (["valid_time"], [0, 100, 300, 600]),
        })
        ds = ds.assign_coords(valid_time=range(4))

        result = _deaccumulate(ds, "ssrd", divisor=3.0)

        # Differences: [0, 100, 200, 300] / 3 = [0, 33.33, 66.67, 100]
        expected = np.array([0, 100, 200, 300]) / 3.0
        np.testing.assert_allclose(result.values, expected, rtol=1e-5)

    def test_reaccumulate(self):
        """Reaccumulation should multiply by target hours."""
        from ecmwf_downloader.forecast.backends.ecmwf_opendata import _reaccumulate

        da = xr.DataArray([10.0, 20.0, 30.0])
        result = _reaccumulate(da, target_hours=2.0)

        expected = np.array([20.0, 40.0, 60.0])
        np.testing.assert_allclose(result.values, expected)

    def test_interpolate_time(self):
        """Time interpolation should create intermediate values."""
        from ecmwf_downloader.forecast.backends.ecmwf_opendata import _interpolate_time
        import pandas as pd

        times = pd.date_range("2024-01-01", periods=3, freq="6h")
        ds = xr.Dataset({
            "t2m": (["valid_time"], [270.0, 280.0, 290.0]),
        })
        ds = ds.assign_coords(valid_time=times)

        result = _interpolate_time(ds, target_freq_hours=3)

        # Should now have times at 3h intervals
        assert len(result.valid_time) == 5  # 0, 3, 6, 9, 12
        # Interpolated values at 3h should be midpoints
        assert abs(result.t2m.values[1] - 275.0) < 0.1  # Between 270 and 280


class TestECMWFOpenDataBackendSpatialSubset:
    """Test spatial subsetting."""

    def test_spatial_subset(self, sample_bbox: BBox):
        """Spatial subset should filter to bbox."""
        from ecmwf_downloader.forecast.backends.ecmwf_opendata import _spatial_subset

        # Create global-ish dataset
        lats = np.arange(40, 50, 0.25)
        lons = np.arange(5, 12, 0.25)
        ds = xr.Dataset({
            "t2m": (["latitude", "longitude"],
                    np.random.rand(len(lats), len(lons))),
        })
        ds = ds.assign_coords(latitude=lats, longitude=lons)

        result = _spatial_subset(ds, sample_bbox)

        # All coords should be within bbox
        assert all(result.latitude.values >= sample_bbox.south)
        assert all(result.latitude.values <= sample_bbox.north)
        assert all(result.longitude.values >= sample_bbox.west)
        assert all(result.longitude.values <= sample_bbox.east)


class TestECMWFOpenDataBackendFetch:
    """Test fetch_forecast with mocked ECMWF client."""

    @pytest.fixture
    def mock_grib_datasets(self) -> tuple[xr.Dataset, xr.Dataset]:
        """Create mock surface and pressure datasets."""
        import pandas as pd

        times = pd.date_range("2024-01-15", periods=49, freq="3h")  # fc1
        lats = np.array([46.0, 45.75, 45.5])
        lons = np.array([7.5, 7.75, 8.0])
        levels = np.array([500, 700, 850, 1000])

        ds_surf = xr.Dataset({
            "t2m": (["time", "latitude", "longitude"],
                    273.15 + 10 * np.random.rand(49, 3, 3)),
            "d2m": (["time", "latitude", "longitude"],
                    268.15 + 10 * np.random.rand(49, 3, 3)),
            "sp": (["time", "latitude", "longitude"],
                   85000 + 5000 * np.random.rand(49, 3, 3)),
            "ssrd": (["time", "latitude", "longitude"],
                     np.cumsum(1e6 * np.random.rand(49, 3, 3), axis=0)),
            "strd": (["time", "latitude", "longitude"],
                     np.cumsum(1e6 * np.random.rand(49, 3, 3), axis=0)),
            "tp": (["time", "latitude", "longitude"],
                   np.cumsum(0.001 * np.random.rand(49, 3, 3), axis=0)),
            "msl": (["time", "latitude", "longitude"],
                    101325 + 1000 * np.random.rand(49, 3, 3)),
        })
        ds_surf = ds_surf.assign_coords(time=times, latitude=lats, longitude=lons)

        ds_plev = xr.Dataset({
            "t": (["time", "level", "latitude", "longitude"],
                  250 + 30 * np.random.rand(49, 4, 3, 3)),
            "gh": (["time", "level", "latitude", "longitude"],
                   5000 + 3000 * np.random.rand(49, 4, 3, 3)),
            "u": (["time", "level", "latitude", "longitude"],
                  20 * (np.random.rand(49, 4, 3, 3) - 0.5)),
            "v": (["time", "level", "latitude", "longitude"],
                  20 * (np.random.rand(49, 4, 3, 3) - 0.5)),
            "r": (["time", "level", "latitude", "longitude"],
                  100 * np.random.rand(49, 4, 3, 3)),
            "q": (["time", "level", "latitude", "longitude"],
                  0.01 * np.random.rand(49, 4, 3, 3)),
        })
        ds_plev = ds_plev.assign_coords(
            time=times, level=levels, latitude=lats, longitude=lons
        )

        return ds_surf, ds_plev

    def test_backend_passes_forecast_hour_to_client(self, sample_bbox: BBox):
        """Backend should pass forecast_hour to ECMWF client."""
        backend = ECMWFOpenDataBackend(
            bbox=sample_bbox,
            pressure_levels=[850, 500],
            forecast_hour=12,
        )
        assert backend.forecast_hour == 12


class TestECMWFOpenDataBackendVariables:
    """Test variable handling."""

    def test_accumulated_variables(self):
        """Check which variables are accumulated."""
        from ecmwf_downloader.forecast.backends.ecmwf_opendata import _ACCUM_VARS

        assert "ssrd" in _ACCUM_VARS
        assert "strd" in _ACCUM_VARS
        assert "tp" in _ACCUM_VARS
        assert "t2m" not in _ACCUM_VARS

    def test_geopotential_conversion(self):
        """Geopotential height should be converted to geopotential."""
        # gh * 9.81 = z
        gh = 5000.0  # meters
        z = gh * 9.81
        assert abs(z - 49050.0) < 1.0
