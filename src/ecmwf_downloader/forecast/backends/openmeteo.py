"""Open-Meteo backend for IFS historical forecasts.

Uses the Open-Meteo historical forecast API to fetch IFS forecast data.
Fast, free, no credentials required. Available from 2022-present.

Note: This backend fetches "historical forecasts" - archived IFS forecasts
that can be queried by valid_time. It does not provide real-time operational
forecasts (use ECMWFOpenDataBackend for that).

Limitations vs ECMWF OpenData:
- No strd (longwave radiation)
- No q (specific humidity) - only r (relative humidity)
- u/v wind has resolution artifacts (see backend-priority.md blacklist)
"""

from __future__ import annotations

import logging
import math
import threading
import time as _time
from collections import deque
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import xarray as xr

from ecmwf_downloader.bbox import BBox
from ecmwf_downloader.forecast.base import ForecastBackend

logger = logging.getLogger(__name__)

# API endpoint
IFS_FORECAST_URL = "https://historical-forecast-api.open-meteo.com/v1/forecast"

# Unit conversions
_C_TO_K = 273.15
_G = 9.80665

# Surface variable mapping: Open-Meteo param -> (ERA5 name, conversion)
SURFACE_MAP = {
    "temperature_2m": ("t2m", lambda x: x + _C_TO_K),
    "dew_point_2m": ("d2m", lambda x: x + _C_TO_K),
    "surface_pressure": ("sp", lambda x: x * 100.0),  # hPa -> Pa
    "shortwave_radiation": ("ssrd", lambda x: x * 3600.0),  # W/m² -> J/m²/h
    "precipitation": ("tp", lambda x: x / 1000.0),  # mm -> m
    "wind_speed_10m": ("ws10", lambda x: x),  # For u10/v10 computation
    "wind_direction_10m": ("wd10", lambda x: x),
}

# Pressure level variable mapping
PLEV_MAP = {
    "temperature": ("t", lambda x: x + _C_TO_K),
    "geopotential_height": ("z", lambda x: x * _G),  # gpm -> m²/s²
    "relative_humidity": ("r", lambda x: x),
}

# Available pressure levels on Open-Meteo IFS
AVAILABLE_LEVELS = [
    1000, 975, 950, 925, 900, 875, 850, 800, 750, 700,
    650, 600, 550, 500, 450, 400, 350, 300, 250, 200,
    150, 100, 70, 50, 30,
]

# Accumulated variables
_ACCUM_VARS = {"ssrd", "tp"}

# Rate limiting
_MAX_REQUESTS_PER_MINUTE = 500


class _RateLimiter:
    """Sliding-window rate limiter."""

    def __init__(self, max_per_minute: int = _MAX_REQUESTS_PER_MINUTE):
        self._window: deque[float] = deque()
        self._max = max_per_minute
        self._lock = threading.Lock()

    def wait(self) -> None:
        with self._lock:
            now = _time.monotonic()
            while self._window and self._window[0] < now - 60:
                self._window.popleft()
            if len(self._window) >= self._max:
                sleep_for = 60.0 - (now - self._window[0])
                if sleep_for > 0:
                    logger.debug("Rate-limiting: sleeping %.1f s", sleep_for)
                    _time.sleep(sleep_for)
            self._window.append(_time.monotonic())


class OpenMeteoIFSBackend(ForecastBackend):
    """Open-Meteo backend for IFS historical forecasts.

    Args:
        bbox: Bounding box (west, south, east, north).
        pressure_levels: Pressure levels in hPa.
        output_timestep: Output time resolution: '1H', '2H', or '3H'.
        forecast_days: Number of forecast days to fetch (default 10).
        batch_size: Number of grid points per API request.

    Note:
        This backend does NOT provide strd (longwave radiation) or q
        (specific humidity). Use ECMWFOpenDataBackend if these are required.
    """

    def __init__(
        self,
        bbox: BBox,
        pressure_levels: list[int] | None = None,
        output_timestep: str = "1H",
        forecast_days: int = 10,
        batch_size: int = 20,
        **kwargs,
    ):
        levels = pressure_levels or [1000, 850, 700, 500, 300]
        super().__init__(bbox, levels, output_timestep, **kwargs)

        # Validate levels
        bad = [l for l in self.pressure_levels if l not in AVAILABLE_LEVELS]
        if bad:
            raise ValueError(
                f"Pressure levels {bad} not available on Open-Meteo IFS. "
                f"Available: {AVAILABLE_LEVELS}"
            )

        self.forecast_days = forecast_days
        self._batch_size = batch_size
        self._rate_limiter = _RateLimiter()

        self._timestep_hours = {"1H": 1, "2H": 2, "3H": 3}.get(output_timestep, 1)

        # Warn about missing variables
        logger.warning(
            "OpenMeteoIFSBackend does not provide strd (longwave) or q (specific humidity). "
            "Use ECMWFOpenDataBackend if these are required."
        )

    @property
    def forecast_horizon_hours(self) -> int:
        return self.forecast_days * 24

    @property
    def available_init_times(self) -> list[int]:
        # Open-Meteo historical API doesn't expose init times directly
        # Returns data by valid_time, so we simulate 00Z
        return [0]

    def _build_grid(self) -> tuple[np.ndarray, np.ndarray]:
        """Build a regular 0.25° grid from the bounding box."""
        step = 0.25
        lat_min = math.ceil(self.bbox.south / step) * step
        lat_max = math.floor(self.bbox.north / step) * step
        lon_min = math.ceil(self.bbox.west / step) * step
        lon_max = math.floor(self.bbox.east / step) * step

        lats = np.arange(lat_min, lat_max + step / 2, step)
        lons = np.arange(lon_min, lon_max + step / 2, step)
        return lats, lons

    def _query_api(
        self,
        lats: list[float],
        lons: list[float],
        start_date: str,
        end_date: str,
    ) -> list[dict]:
        """Query Open-Meteo API for a batch of grid points."""
        try:
            import requests
        except ImportError:
            raise ImportError(
                "The 'requests' package is required for OpenMeteoIFSBackend. "
                "Install with: pip install requests"
            )
        from tenacity import retry, stop_after_attempt, wait_exponential

        # Build parameter list
        hourly_params = list(SURFACE_MAP.keys())

        for level in self.pressure_levels:
            for base_param in PLEV_MAP:
                hourly_params.append(f"{base_param}_{level}hPa")
            # Wind at pressure levels (for u/v computation)
            hourly_params.append(f"wind_speed_{level}hPa")
            hourly_params.append(f"wind_direction_{level}hPa")

        params = {
            "latitude": ",".join(str(lat) for lat in lats),
            "longitude": ",".join(str(lon) for lon in lons),
            "hourly": ",".join(hourly_params),
            "start_date": start_date,
            "end_date": end_date,
            "timeformat": "unixtime",
        }

        @retry(
            stop=stop_after_attempt(3),
            wait=wait_exponential(min=2, max=30),
            reraise=True,
        )
        def _do_request():
            self._rate_limiter.wait()
            resp = requests.get(IFS_FORECAST_URL, params=params, timeout=120)
            resp.raise_for_status()
            return resp.json()

        data = _do_request()

        # Handle single vs multiple points
        if isinstance(data, dict) and "latitude" in data:
            if not isinstance(data["latitude"], list):
                return [data]
        if isinstance(data, list):
            return data
        return [data]

    def _assemble_grid(
        self,
        responses: list[dict],
        lats: np.ndarray,
        lons: np.ndarray,
    ) -> tuple[xr.Dataset, xr.Dataset]:
        """Assemble API responses into gridded xarray Datasets."""
        # Parse time from first response
        first = responses[0]
        times_unix = first["hourly"]["time"]
        times = pd.to_datetime(times_unix, unit="s", utc=True).tz_localize(None)
        n_times = len(times)
        n_lats = len(lats)
        n_lons = len(lons)

        # Surface variables
        surf_arrays: dict[str, np.ndarray] = {}
        for om_param, (era5_name, convert) in SURFACE_MAP.items():
            if era5_name in ("ws10", "wd10"):
                continue  # Handle wind separately
            arr = np.full((n_times, n_lats, n_lons), np.nan, dtype="float32")
            for idx, resp in enumerate(responses):
                ilat = idx // n_lons
                ilon = idx % n_lons
                values = resp["hourly"].get(om_param)
                if values is not None:
                    col = np.array(values, dtype="float32")
                    col = np.where(np.isnan(col), np.nan, convert(col))
                    arr[:, ilat, ilon] = col
            surf_arrays[era5_name] = arr

        # 10m wind: speed + direction -> u10, v10
        u10_arr = np.full((n_times, n_lats, n_lons), np.nan, dtype="float32")
        v10_arr = np.full((n_times, n_lats, n_lons), np.nan, dtype="float32")
        for idx, resp in enumerate(responses):
            ilat = idx // n_lons
            ilon = idx % n_lons
            speed = resp["hourly"].get("wind_speed_10m")
            direction = resp["hourly"].get("wind_direction_10m")
            if speed is not None and direction is not None:
                speed = np.array(speed, dtype="float32")
                direction = np.array(direction, dtype="float32")
                dir_rad = np.deg2rad(direction)
                u10_arr[:, ilat, ilon] = -speed * np.sin(dir_rad)
                v10_arr[:, ilat, ilon] = -speed * np.cos(dir_rad)
        surf_arrays["u10"] = u10_arr
        surf_arrays["v10"] = v10_arr

        # Surface geopotential from elevation
        z_surf = np.full((n_times, n_lats, n_lons), np.nan, dtype="float32")
        for idx, resp in enumerate(responses):
            ilat = idx // n_lons
            ilon = idx % n_lons
            elev = resp.get("elevation")
            if elev is not None:
                z_surf[:, ilat, ilon] = float(elev) * _G
        surf_arrays["z"] = z_surf

        ds_surf = xr.Dataset(
            {
                name: (["valid_time", "latitude", "longitude"], data)
                for name, data in surf_arrays.items()
            },
            coords={
                "valid_time": times,
                "latitude": lats,
                "longitude": lons,
            },
        )

        # Pressure level variables
        levels = np.array(sorted(self.pressure_levels), dtype="int32")
        n_levels = len(levels)

        plev_arrays: dict[str, np.ndarray] = {}

        for base_param, (era5_name, convert) in PLEV_MAP.items():
            arr = np.full(
                (n_times, n_levels, n_lats, n_lons), np.nan, dtype="float32"
            )
            for li, level in enumerate(levels):
                om_key = f"{base_param}_{level}hPa"
                for idx, resp in enumerate(responses):
                    ilat = idx // n_lons
                    ilon = idx % n_lons
                    values = resp["hourly"].get(om_key)
                    if values is not None:
                        col = np.array(values, dtype="float32")
                        col = np.where(np.isnan(col), np.nan, convert(col))
                        arr[:, li, ilat, ilon] = col
            plev_arrays[era5_name] = arr

        # Pressure level wind (note: blacklisted due to resolution issues)
        u_arr = np.full((n_times, n_levels, n_lats, n_lons), np.nan, dtype="float32")
        v_arr = np.full((n_times, n_levels, n_lats, n_lons), np.nan, dtype="float32")
        for li, level in enumerate(levels):
            speed_key = f"wind_speed_{level}hPa"
            dir_key = f"wind_direction_{level}hPa"
            for idx, resp in enumerate(responses):
                ilat = idx // n_lons
                ilon = idx % n_lons
                speed = resp["hourly"].get(speed_key)
                direction = resp["hourly"].get(dir_key)
                if speed is not None and direction is not None:
                    speed = np.array(speed, dtype="float32")
                    direction = np.array(direction, dtype="float32")
                    dir_rad = np.deg2rad(direction)
                    u_arr[:, li, ilat, ilon] = -speed * np.sin(dir_rad)
                    v_arr[:, li, ilat, ilon] = -speed * np.cos(dir_rad)
        plev_arrays["u"] = u_arr
        plev_arrays["v"] = v_arr

        ds_plev = xr.Dataset(
            {
                name: (["valid_time", "level", "latitude", "longitude"], data)
                for name, data in plev_arrays.items()
            },
            coords={
                "valid_time": times,
                "level": levels,
                "latitude": lats,
                "longitude": lons,
            },
        )

        return ds_surf, ds_plev

    def fetch_forecast(
        self, init_time: datetime
    ) -> tuple[xr.Dataset, xr.Dataset]:
        """Fetch IFS forecast data from Open-Meteo.

        Args:
            init_time: Forecast initialization time (used as start date).

        Returns:
            Tuple of (ds_surf, ds_plev) with valid_time dimension.
        """
        lats, lons = self._build_grid()

        start_date = init_time.strftime("%Y-%m-%d")
        end_date = (init_time + timedelta(days=self.forecast_days)).strftime("%Y-%m-%d")

        logger.info(
            "Open-Meteo IFS: fetching %s to %s for %d×%d grid",
            start_date,
            end_date,
            len(lats),
            len(lons),
        )

        # Build all (lat, lon) pairs
        grid_lats = []
        grid_lons = []
        for lat in lats:
            for lon in lons:
                grid_lats.append(float(lat))
                grid_lons.append(float(lon))

        # Batch API calls
        all_responses: list[dict] = []
        for i in range(0, len(grid_lats), self._batch_size):
            batch_lats = grid_lats[i : i + self._batch_size]
            batch_lons = grid_lons[i : i + self._batch_size]
            logger.debug(
                "Fetching batch %d/%d (%d points)",
                i // self._batch_size + 1,
                math.ceil(len(grid_lats) / self._batch_size),
                len(batch_lats),
            )
            responses = self._query_api(batch_lats, batch_lons, start_date, end_date)
            all_responses.extend(responses)

        # Assemble into gridded datasets
        ds_surf, ds_plev = self._assemble_grid(all_responses, lats, lons)

        # Sum accumulated variables over timestep window
        if self._timestep_hours > 1:
            for var in _ACCUM_VARS:
                if var in ds_surf:
                    ds_surf[var] = (
                        ds_surf[var]
                        .rolling(valid_time=self._timestep_hours, min_periods=1)
                        .sum()
                    )

        # Subsample to target timestep
        target_hours = list(range(0, 24, self._timestep_hours))
        ds_surf = ds_surf.sel(
            valid_time=ds_surf.valid_time.dt.hour.isin(target_hours)
        )
        ds_plev = ds_plev.sel(
            valid_time=ds_plev.valid_time.dt.hour.isin(target_hours)
        )

        # Add init_time as coordinate
        ds_surf = ds_surf.assign_coords(init_time=init_time)
        ds_plev = ds_plev.assign_coords(init_time=init_time)

        logger.info(
            "Open-Meteo IFS: fetched SURF vars=%s, PLEV vars=%s",
            list(ds_surf.data_vars),
            list(ds_plev.data_vars),
        )

        return ds_surf.load(), ds_plev.load()
