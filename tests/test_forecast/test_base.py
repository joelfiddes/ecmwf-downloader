"""Tests for ForecastBackend abstract base class."""

from __future__ import annotations

from datetime import datetime

import pytest
import xarray as xr

from ecmwf_downloader.bbox import BBox
from ecmwf_downloader.forecast.base import ForecastBackend


class DummyBackend(ForecastBackend):
    """Concrete implementation for testing the ABC."""

    def fetch_forecast(self, init_time: datetime) -> tuple[xr.Dataset, xr.Dataset]:
        return xr.Dataset(), xr.Dataset()

    @property
    def forecast_horizon_hours(self) -> int:
        return 240

    @property
    def available_init_times(self) -> list[int]:
        return [0, 12]


class TestForecastBackendABC:
    """Test the ForecastBackend abstract base class."""

    def test_init_stores_bbox(self, sample_bbox: BBox):
        """Backend should store bbox."""
        backend = DummyBackend(
            bbox=sample_bbox,
            pressure_levels=[1000, 850, 700],
        )
        assert backend.bbox == sample_bbox

    def test_init_sorts_pressure_levels(self, sample_bbox: BBox):
        """Backend should sort pressure levels ascending."""
        backend = DummyBackend(
            bbox=sample_bbox,
            pressure_levels=[1000, 500, 850, 700],
        )
        assert backend.pressure_levels == [500, 700, 850, 1000]

    def test_init_default_timestep(self, sample_bbox: BBox):
        """Backend should default to 1H timestep."""
        backend = DummyBackend(
            bbox=sample_bbox,
            pressure_levels=[1000, 850],
        )
        assert backend.output_timestep == "1H"

    def test_init_custom_timestep(self, sample_bbox: BBox):
        """Backend should accept custom timestep."""
        backend = DummyBackend(
            bbox=sample_bbox,
            pressure_levels=[1000, 850],
            output_timestep="3H",
        )
        assert backend.output_timestep == "3H"

    def test_forecast_horizon_hours(self, sample_bbox: BBox):
        """Backend should report forecast horizon."""
        backend = DummyBackend(
            bbox=sample_bbox,
            pressure_levels=[1000, 850],
        )
        assert backend.forecast_horizon_hours == 240

    def test_available_init_times(self, sample_bbox: BBox):
        """Backend should report available init times."""
        backend = DummyBackend(
            bbox=sample_bbox,
            pressure_levels=[1000, 850],
        )
        assert backend.available_init_times == [0, 12]

    def test_close_is_callable(self, sample_bbox: BBox):
        """Backend close() should be callable (no-op by default)."""
        backend = DummyBackend(
            bbox=sample_bbox,
            pressure_levels=[1000, 850],
        )
        # Should not raise
        backend.close()


class TestForecastBackendABCEnforcement:
    """Test that ABC properly enforces interface."""

    def test_cannot_instantiate_abc_directly(self, sample_bbox: BBox):
        """Cannot instantiate ForecastBackend directly."""
        with pytest.raises(TypeError, match="abstract"):
            ForecastBackend(
                bbox=sample_bbox,
                pressure_levels=[1000, 850],
            )

    def test_must_implement_fetch_forecast(self, sample_bbox: BBox):
        """Subclass must implement fetch_forecast."""

        class IncompleteBackend(ForecastBackend):
            @property
            def forecast_horizon_hours(self) -> int:
                return 240

            @property
            def available_init_times(self) -> list[int]:
                return [0, 12]

        with pytest.raises(TypeError, match="fetch_forecast"):
            IncompleteBackend(
                bbox=sample_bbox,
                pressure_levels=[1000, 850],
            )

    def test_must_implement_forecast_horizon(self, sample_bbox: BBox):
        """Subclass must implement forecast_horizon_hours."""

        class IncompleteBackend(ForecastBackend):
            def fetch_forecast(self, init_time):
                return xr.Dataset(), xr.Dataset()

            @property
            def available_init_times(self) -> list[int]:
                return [0, 12]

        with pytest.raises(TypeError, match="forecast_horizon_hours"):
            IncompleteBackend(
                bbox=sample_bbox,
                pressure_levels=[1000, 850],
            )

    def test_must_implement_available_init_times(self, sample_bbox: BBox):
        """Subclass must implement available_init_times."""

        class IncompleteBackend(ForecastBackend):
            def fetch_forecast(self, init_time):
                return xr.Dataset(), xr.Dataset()

            @property
            def forecast_horizon_hours(self) -> int:
                return 240

        with pytest.raises(TypeError, match="available_init_times"):
            IncompleteBackend(
                bbox=sample_bbox,
                pressure_levels=[1000, 850],
            )
