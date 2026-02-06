"""Unified forecast loader with backend selection.

Orchestrates forecast downloads from multiple backends (ECMWF OpenData,
Open-Meteo IFS) with automatic priority selection and fallback.
"""

from __future__ import annotations

import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Literal

import pandas as pd
import xarray as xr

from ecmwf_downloader.bbox import BBox
from ecmwf_downloader.forecast.base import ForecastBackend
from ecmwf_downloader.forecast.backends import ECMWFOpenDataBackend, OpenMeteoIFSBackend
from ecmwf_downloader.writer import ZarrWriter

logger = logging.getLogger(__name__)

# Variables that are blacklisted from OpenMeteo (see backend-priority.md)
_OPENMETEO_BLACKLIST = {"strd", "q"}

# Backend registry
_BACKENDS = {
    "ecmwf_opendata": ECMWFOpenDataBackend,
    "openmeteo_ifs": OpenMeteoIFSBackend,
}


class ForecastLoader:
    """Unified forecast loader with backend selection.

    Supports automatic backend selection based on priority strategy,
    with fallback when preferred backend fails.

    Args:
        bbox: Bounding box as (W, S, E, N) tuple or BBox instance.
        output_dir: Directory for output Zarr store.
        backend: Backend to use: 'ecmwf_opendata', 'openmeteo_ifs', or 'auto'.
        priority_strategy: For auto mode: 'speed' or 'reliability'.
        output_timestep: Output time resolution: '1H', '2H', or '3H'.
        pressure_levels: Pressure levels in hPa.
        required_variables: Variables that must be present (affects backend selection).
        forecast_hour: Forecast initialization hour (0 or 12) for ECMWF.
        backend_kwargs: Additional kwargs passed to the backend.

    Example::

        loader = ForecastLoader(
            bbox=(7.7, 45.95, 7.85, 46.05),
            output_dir="./forecast/",
            backend="auto",
            priority_strategy="speed",
        )
        loader.download()
        ds = loader.open()
    """

    def __init__(
        self,
        bbox: tuple | BBox,
        output_dir: str = "./forecast/",
        backend: Literal["ecmwf_opendata", "openmeteo_ifs", "auto"] = "auto",
        priority_strategy: Literal["speed", "reliability"] = "speed",
        output_timestep: str = "1H",
        pressure_levels: list[int] | None = None,
        required_variables: list[str] | None = None,
        forecast_hour: int = 0,
        backend_kwargs: dict | None = None,
    ):
        if isinstance(bbox, (tuple, list)):
            self.bbox = BBox.from_tuple(tuple(bbox))
        else:
            self.bbox = bbox

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.backend_name = backend
        self.priority_strategy = priority_strategy
        self.output_timestep = output_timestep
        self.pressure_levels = pressure_levels or [1000, 850, 700, 500, 300]
        self.required_variables = set(required_variables or [])
        self.forecast_hour = forecast_hour
        self.backend_kwargs = backend_kwargs or {}

        self._backend: ForecastBackend | None = None
        self._writer: ZarrWriter | None = None

    def _select_backend(self) -> str:
        """Select backend based on priority strategy and requirements."""
        if self.backend_name != "auto":
            return self.backend_name

        # Check if any required variables are blacklisted from OpenMeteo
        blacklisted = self.required_variables & _OPENMETEO_BLACKLIST
        if blacklisted:
            logger.info(
                "Required variables %s not available on OpenMeteo, using ecmwf_opendata",
                blacklisted,
            )
            return "ecmwf_opendata"

        # Apply priority strategy
        if self.priority_strategy == "speed":
            # OpenMeteo is faster
            return "openmeteo_ifs"
        else:
            # ECMWF is more reliable/complete
            return "ecmwf_opendata"

    def _get_backend(self) -> ForecastBackend:
        """Get or create the backend instance."""
        if self._backend is not None:
            return self._backend

        backend_name = self._select_backend()
        logger.info("Using forecast backend: %s", backend_name)

        backend_cls = _BACKENDS[backend_name]

        kwargs = {
            "bbox": self.bbox,
            "pressure_levels": self.pressure_levels,
            "output_timestep": self.output_timestep,
            **self.backend_kwargs,
        }

        if backend_name == "ecmwf_opendata":
            kwargs["forecast_hour"] = self.forecast_hour

        self._backend = backend_cls(**kwargs)
        return self._backend

    def _get_writer(self) -> ZarrWriter:
        """Get or create the Zarr writer."""
        if self._writer is None:
            self._writer = ZarrWriter(
                self.output_dir,
                chunks={
                    "valid_time": 240,  # 10 days at hourly
                    "latitude": 10,
                    "longitude": 10,
                    "level": -1,
                },
            )
        return self._writer

    def _get_existing_times(self) -> set[datetime]:
        """Get forecast init times already in the Zarr store."""
        zarr_path = self.output_dir / "forecast.zarr"
        if not zarr_path.exists():
            return set()

        try:
            ds = xr.open_zarr(str(zarr_path))
            if "init_time" in ds.coords:
                times = set(pd.Timestamp(t).to_pydatetime() for t in ds.init_time.values)
                ds.close()
                return times
        except Exception:
            pass

        return set()

    def download(
        self,
        init_time: datetime | None = None,
        backfill_days: int = 3,
    ) -> None:
        """Download forecast data.

        Args:
            init_time: Specific initialization time to download. If None,
                downloads today's forecast.
            backfill_days: Number of historical days to backfill (default 3,
                matching ECMWF open data retention).
        """
        backend = self._get_backend()

        # Determine which forecasts to download
        if init_time is not None:
            init_times = [init_time]
        else:
            today = datetime.now().replace(
                hour=self.forecast_hour, minute=0, second=0, microsecond=0
            )
            init_times = []

            # Backfill historical
            for days_ago in range(backfill_days, 0, -1):
                init_times.append(today - timedelta(days=days_ago))

            # Today's forecast
            init_times.append(today)

        # Filter out already downloaded
        existing = self._get_existing_times()
        to_download = [t for t in init_times if t not in existing]

        if not to_download:
            logger.info("All requested forecasts already downloaded")
            return

        logger.info("Downloading %d forecast(s)", len(to_download))

        for init in to_download:
            logger.info("Fetching forecast for %s", init.strftime("%Y-%m-%d %H:%M"))
            try:
                ds_surf, ds_plev = backend.fetch_forecast(init)
                self._write_forecast(ds_surf, ds_plev, init)
            except Exception as e:
                logger.error("Failed to fetch forecast for %s: %s", init, e)
                # Try fallback if using auto
                if self.backend_name == "auto":
                    self._try_fallback(init)

    def _try_fallback(self, init_time: datetime) -> None:
        """Try fallback backend after primary failure."""
        current = self._select_backend()
        fallback = "ecmwf_opendata" if current == "openmeteo_ifs" else "openmeteo_ifs"

        # Check if fallback can provide required variables
        if fallback == "openmeteo_ifs" and (self.required_variables & _OPENMETEO_BLACKLIST):
            logger.warning("Fallback to openmeteo_ifs not possible due to required variables")
            return

        logger.info("Trying fallback backend: %s", fallback)

        try:
            fallback_cls = _BACKENDS[fallback]
            kwargs = {
                "bbox": self.bbox,
                "pressure_levels": self.pressure_levels,
                "output_timestep": self.output_timestep,
                **self.backend_kwargs,
            }
            if fallback == "ecmwf_opendata":
                kwargs["forecast_hour"] = self.forecast_hour

            fallback_backend = fallback_cls(**kwargs)
            ds_surf, ds_plev = fallback_backend.fetch_forecast(init_time)
            self._write_forecast(ds_surf, ds_plev, init_time)
            fallback_backend.close()
        except Exception as e:
            logger.error("Fallback also failed: %s", e)

    def _write_forecast(
        self,
        ds_surf: xr.Dataset,
        ds_plev: xr.Dataset,
        init_time: datetime,
    ) -> None:
        """Write forecast to Zarr store."""
        zarr_path = self.output_dir / "forecast.zarr"

        # Merge surface and pressure level data
        # Rename z in surface to z_surf to avoid clash
        if "z" in ds_surf:
            ds_surf = ds_surf.rename({"z": "z_surf"})

        ds = xr.merge([ds_surf, ds_plev])

        if not zarr_path.exists():
            # Create new store
            ds.to_zarr(str(zarr_path), mode="w")
            logger.info("Created forecast Zarr store: %s", zarr_path)
        else:
            # Append to existing
            ds.to_zarr(str(zarr_path), mode="a", append_dim="valid_time")
            logger.info("Appended forecast to Zarr store")

    def open(self) -> xr.Dataset:
        """Open the forecast Zarr store as a lazy Dataset.

        Returns:
            xarray Dataset with lazy loading from Zarr.
        """
        zarr_path = self.output_dir / "forecast.zarr"
        if not zarr_path.exists():
            raise FileNotFoundError(
                f"No forecast data found at {zarr_path}. Run download() first."
            )
        return xr.open_zarr(str(zarr_path))

    def close(self):
        """Clean up resources."""
        if self._backend is not None:
            self._backend.close()
            self._backend = None
