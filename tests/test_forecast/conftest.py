"""Shared fixtures for forecast module tests."""

from __future__ import annotations

from datetime import datetime
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from ecmwf_downloader.bbox import BBox


@pytest.fixture
def sample_bbox() -> BBox:
    """Sample bounding box for testing (Zermatt area, Switzerland)."""
    return BBox.from_tuple((7.7, 45.95, 7.85, 46.05))


@pytest.fixture
def sample_pressure_levels() -> list[int]:
    """Standard pressure levels for testing."""
    return [1000, 850, 700, 500, 300]


@pytest.fixture
def sample_init_time() -> datetime:
    """Sample forecast initialization time."""
    return datetime(2024, 1, 15, 0, 0, 0)


@pytest.fixture
def sample_surface_dataset(sample_init_time: datetime) -> xr.Dataset:
    """Create a sample surface forecast dataset."""
    n_times = 24
    n_lats = 3
    n_lons = 3

    valid_times = pd.date_range(
        sample_init_time, periods=n_times, freq="1h"
    )
    lats = np.array([46.05, 46.0, 45.95])
    lons = np.array([7.7, 7.775, 7.85])

    # Create sample data
    np.random.seed(42)

    ds = xr.Dataset(
        {
            "t2m": (["valid_time", "latitude", "longitude"],
                    273.15 + 10 * np.random.rand(n_times, n_lats, n_lons)),
            "d2m": (["valid_time", "latitude", "longitude"],
                    268.15 + 10 * np.random.rand(n_times, n_lats, n_lons)),
            "sp": (["valid_time", "latitude", "longitude"],
                   85000 + 5000 * np.random.rand(n_times, n_lats, n_lons)),
            "ssrd": (["valid_time", "latitude", "longitude"],
                     1e6 * np.random.rand(n_times, n_lats, n_lons)),
            "tp": (["valid_time", "latitude", "longitude"],
                   0.001 * np.random.rand(n_times, n_lats, n_lons)),
            "u10": (["valid_time", "latitude", "longitude"],
                    5 * (np.random.rand(n_times, n_lats, n_lons) - 0.5)),
            "v10": (["valid_time", "latitude", "longitude"],
                    5 * (np.random.rand(n_times, n_lats, n_lons) - 0.5)),
            "z": (["valid_time", "latitude", "longitude"],
                  9.81 * 1500 * np.ones((n_times, n_lats, n_lons))),
        },
        coords={
            "valid_time": valid_times,
            "latitude": lats,
            "longitude": lons,
            "init_time": sample_init_time,
        },
    )
    return ds


@pytest.fixture
def sample_pressure_dataset(
    sample_init_time: datetime,
    sample_pressure_levels: list[int]
) -> xr.Dataset:
    """Create a sample pressure-level forecast dataset."""
    n_times = 24
    n_levels = len(sample_pressure_levels)
    n_lats = 3
    n_lons = 3

    valid_times = pd.date_range(
        sample_init_time, periods=n_times, freq="1h"
    )
    lats = np.array([46.05, 46.0, 45.95])
    lons = np.array([7.7, 7.775, 7.85])
    levels = np.array(sorted(sample_pressure_levels))

    np.random.seed(42)

    ds = xr.Dataset(
        {
            "t": (["valid_time", "level", "latitude", "longitude"],
                  250 + 30 * np.random.rand(n_times, n_levels, n_lats, n_lons)),
            "z": (["valid_time", "level", "latitude", "longitude"],
                  9.81 * (5000 + 3000 * np.random.rand(n_times, n_levels, n_lats, n_lons))),
            "u": (["valid_time", "level", "latitude", "longitude"],
                  20 * (np.random.rand(n_times, n_levels, n_lats, n_lons) - 0.5)),
            "v": (["valid_time", "level", "latitude", "longitude"],
                  20 * (np.random.rand(n_times, n_levels, n_lats, n_lons) - 0.5)),
            "r": (["valid_time", "level", "latitude", "longitude"],
                  100 * np.random.rand(n_times, n_levels, n_lats, n_lons)),
        },
        coords={
            "valid_time": valid_times,
            "level": levels,
            "latitude": lats,
            "longitude": lons,
            "init_time": sample_init_time,
        },
    )
    return ds


@pytest.fixture
def mock_openmeteo_response() -> list[dict]:
    """Create a mock Open-Meteo API response for a single grid point."""
    # 24 hours of data
    base_time = 1705276800  # 2024-01-15 00:00:00 UTC
    times = [base_time + i * 3600 for i in range(24)]

    np.random.seed(42)

    return [{
        "latitude": 46.0,
        "longitude": 7.8,
        "elevation": 1500.0,
        "hourly": {
            "time": times,
            "temperature_2m": list(np.random.rand(24) * 10 - 5),  # -5 to 5 °C
            "dew_point_2m": list(np.random.rand(24) * 10 - 10),  # -10 to 0 °C
            "surface_pressure": list(850 + np.random.rand(24) * 20),  # hPa
            "shortwave_radiation": list(np.random.rand(24) * 500),  # W/m²
            "precipitation": list(np.random.rand(24) * 2),  # mm
            "wind_speed_10m": list(np.random.rand(24) * 10),  # m/s
            "wind_direction_10m": list(np.random.rand(24) * 360),  # degrees
            # Pressure level vars at 850 hPa
            "temperature_850hPa": list(np.random.rand(24) * 20 - 10),
            "geopotential_height_850hPa": list(1400 + np.random.rand(24) * 100),
            "relative_humidity_850hPa": list(np.random.rand(24) * 100),
            "wind_speed_850hPa": list(np.random.rand(24) * 30),
            "wind_direction_850hPa": list(np.random.rand(24) * 360),
        },
    }]
