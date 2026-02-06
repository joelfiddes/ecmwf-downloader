"""Abstract base class for forecast data backends."""

from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime

import xarray as xr

from ecmwf_downloader.bbox import BBox


class ForecastBackend(ABC):
    """Abstract base class for forecast data backends.

    Each backend must implement ``fetch_forecast()`` which returns a pair of
    xarray Datasets (surface, pressure-level) for a single forecast run.

    Standard surface variables: t2m, d2m, sp, ssrd, strd, tp, z, u10, v10
    Standard pressure variables: t, z, u, v, q, r
    Standard dimensions: valid_time, latitude, longitude, level (hPa ascending)

    Note: Forecast data uses ``valid_time`` (when the forecast is valid for),
    not ``time`` (when the forecast was issued). The ``init_time`` parameter
    specifies when the forecast was initialized.
    """

    def __init__(
        self,
        bbox: BBox,
        pressure_levels: list[int],
        output_timestep: str = "1H",
        **kwargs,
    ):
        """Initialize forecast backend.

        Args:
            bbox: Bounding box (west, south, east, north).
            pressure_levels: Pressure levels in hPa.
            output_timestep: Output time resolution: '1H', '2H', or '3H'.
        """
        self.bbox = bbox
        self.pressure_levels = sorted(pressure_levels)
        self.output_timestep = output_timestep

    @abstractmethod
    def fetch_forecast(
        self, init_time: datetime
    ) -> tuple[xr.Dataset, xr.Dataset]:
        """Fetch a single forecast run.

        Args:
            init_time: Forecast initialization time (base time).

        Returns:
            Tuple of (ds_surf, ds_plev) with:
            - ``valid_time`` dimension (forecast valid times)
            - Standardized variable names
            - Accumulated variables (ssrd, strd, tp) per output timestep
        """
        ...

    @property
    @abstractmethod
    def forecast_horizon_hours(self) -> int:
        """Maximum forecast lead time in hours."""
        ...

    @property
    @abstractmethod
    def available_init_times(self) -> list[int]:
        """Available forecast initialization hours (UTC).

        E.g., [0, 12] for 00Z and 12Z forecasts.
        """
        ...

    def close(self):
        """Clean up resources. Override in subclasses if needed."""
        pass
