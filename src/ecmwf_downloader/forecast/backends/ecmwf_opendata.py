"""ECMWF Open Data backend for IFS forecasts.

Downloads IFS forecasts from the official ECMWF open data service using the
ecmwf-opendata client. Data is in GRIB format, requires cfgrib.

IFS open data structure:
    - For times 00z & 12z: steps 0-144 by 3h, steps 150-240 by 6h
    - Surface params: 2t, sp, 2d, ssrd, strd, tp, msl
    - Pressure params: gh, u, v, r, q, t
    - Default levels: 1000, 925, 850, 700, 600, 500, 400, 300
"""

from __future__ import annotations

import logging
import shutil
import tempfile
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

from ecmwf_downloader.bbox import BBox
from ecmwf_downloader.derived import compute_surface_geopotential
from ecmwf_downloader.forecast.base import ForecastBackend
from ecmwf_downloader.variables import (
    IFS_DEFAULT_LEVELS,
    IFS_PLEV_PARAMS,
    IFS_RENAME_MAP,
    IFS_SURF_PARAMS,
)

logger = logging.getLogger(__name__)

# Accumulated variables that need deaccumulation/re-accumulation
_ACCUM_VARS = {"ssrd", "strd", "tp"}

# Output timestep in hours
_TIMESTEP_HOURS = {"1H": 1, "2H": 2, "3H": 3}


def _open_grib_surface(path: str) -> xr.Dataset:
    """Open a surface GRIB file with cfgrib, merging all variables."""
    import cfgrib

    datasets = cfgrib.open_datasets(path)
    if len(datasets) == 1:
        return datasets[0]
    return xr.merge(datasets, compat="override")


def _open_grib_pressure(path: str) -> xr.Dataset:
    """Open a pressure-level GRIB file with cfgrib."""
    import cfgrib

    datasets = cfgrib.open_datasets(path)
    if len(datasets) == 1:
        return datasets[0]
    return xr.merge(datasets, compat="override")


def _spatial_subset(ds: xr.Dataset, bbox: BBox) -> xr.Dataset:
    """Subset dataset to bounding box."""
    lat_name = "latitude" if "latitude" in ds.dims else "lat"
    lon_name = "longitude" if "longitude" in ds.dims else "lon"

    lat = ds[lat_name].values
    lon = ds[lon_name].values

    lat_mask = (lat >= bbox.south) & (lat <= bbox.north)
    lon_mask = (lon >= bbox.west) & (lon <= bbox.east)

    return ds.sel({lat_name: lat[lat_mask], lon_name: lon[lon_mask]})


def _deaccumulate(ds: xr.Dataset, var_name: str, divisor: float) -> xr.DataArray:
    """Deaccumulate a forecast variable: diff between consecutive steps."""
    acc = ds[var_name]
    deacc = acc - acc.shift(valid_time=1, fill_value=0)
    return deacc / divisor


def _interpolate_time(ds: xr.Dataset, target_freq_hours: int) -> xr.Dataset:
    """Interpolate dataset to a finer time resolution."""
    times = pd.DatetimeIndex(ds.valid_time.values)
    new_times = pd.date_range(times[0], times[-1], freq=f"{target_freq_hours}h")
    return ds.interp(valid_time=new_times, method="linear")


def _reaccumulate(da: xr.DataArray, target_hours: float) -> xr.DataArray:
    """Re-accumulate a rate variable to match target timestep."""
    return da * target_hours


class ECMWFOpenDataBackend(ForecastBackend):
    """ECMWF Open Data backend for IFS forecasts.

    Uses the official ecmwf-opendata client to download GRIB files.

    Args:
        bbox: Bounding box (west, south, east, north).
        pressure_levels: Pressure levels in hPa.
        output_timestep: Output time resolution: '1H', '2H', or '3H'.
        forecast_hour: Forecast initialization hour: 0 or 12 (UTC).
    """

    def __init__(
        self,
        bbox: BBox,
        pressure_levels: list[int] | None = None,
        output_timestep: str = "1H",
        forecast_hour: int = 0,
        show_progress: bool = True,
        **kwargs,
    ):
        levels = pressure_levels or IFS_DEFAULT_LEVELS
        super().__init__(bbox, levels, output_timestep, **kwargs)

        if output_timestep not in _TIMESTEP_HOURS:
            raise ValueError(
                f"Invalid output_timestep '{output_timestep}'. "
                f"Must be one of: {list(_TIMESTEP_HOURS.keys())}"
            )
        self._target_hours = _TIMESTEP_HOURS[output_timestep]

        if forecast_hour not in (0, 12):
            raise ValueError("forecast_hour must be 0 or 12")
        self.forecast_hour = forecast_hour
        self.show_progress = show_progress

    @property
    def forecast_horizon_hours(self) -> int:
        return 240  # 10 days

    @property
    def available_init_times(self) -> list[int]:
        return [0, 12]

    def _download_grib(self, init_time: datetime, tmp_dir: Path) -> None:
        """Download GRIB files for a forecast run."""
        from ecmwf.opendata import Client

        client = Client()

        # Calculate date offset from today
        today = datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)
        delta_days = (init_time.date() - today.date()).days

        fc1_steps = list(range(0, 147, 3))  # 0 to 144 by 3h
        fc2_steps = list(range(150, 241, 6))  # 150 to 240 by 6h

        logger.info(
            "Downloading IFS forecast GRIB for %s %02dZ",
            init_time.strftime("%Y-%m-%d"),
            self.forecast_hour,
        )

        downloads = [
            ("SURF fc1", dict(
                time=self.forecast_hour,
                date=delta_days,
                step=fc1_steps,
                type="fc",
                param=IFS_SURF_PARAMS,
                target=str(tmp_dir / "SURF_fc1.grib2"),
            )),
            ("PLEV fc1", dict(
                time=self.forecast_hour,
                date=delta_days,
                step=fc1_steps,
                type="fc",
                param=IFS_PLEV_PARAMS,
                levelist=self.pressure_levels,
                target=str(tmp_dir / "PLEV_fc1.grib2"),
            )),
            ("SURF fc2", dict(
                time=self.forecast_hour,
                date=delta_days,
                step=fc2_steps,
                type="fc",
                param=IFS_SURF_PARAMS,
                target=str(tmp_dir / "SURF_fc2.grib2"),
            )),
            ("PLEV fc2", dict(
                time=self.forecast_hour,
                date=delta_days,
                step=fc2_steps,
                type="fc",
                param=IFS_PLEV_PARAMS,
                levelist=self.pressure_levels,
                target=str(tmp_dir / "PLEV_fc2.grib2"),
            )),
        ]

        iterator = downloads
        if self.show_progress:
            iterator = tqdm(
                downloads,
                desc=f"IFS {init_time.strftime('%Y-%m-%d')} {self.forecast_hour:02d}Z",
                unit="file",
                leave=False,
            )

        for name, params in iterator:
            client.retrieve(**params)

    def _process_surface(
        self, fc1_path: Path, fc2_path: Path
    ) -> xr.Dataset:
        """Process surface GRIB files."""
        ds1 = _spatial_subset(_open_grib_surface(str(fc1_path)), self.bbox)
        ds2 = _spatial_subset(_open_grib_surface(str(fc2_path)), self.bbox)

        # Rename time -> valid_time for forecast convention
        if "time" in ds1.dims:
            ds1 = ds1.rename({"time": "valid_time"})
        if "time" in ds2.dims:
            ds2 = ds2.rename({"time": "valid_time"})

        # Handle tp variable name quirk
        for ds in [ds1, ds2]:
            if "tp" not in ds and "param193.1.0" in ds:
                ds["tp"] = ds["param193.1.0"]

        # Save last values of fc1 accumulators for fc2 boundary correction
        fc1_last = {}
        for var in _ACCUM_VARS:
            if var in ds1:
                fc1_last[var] = ds1[var].isel(valid_time=-1)

        # Deaccumulate fc1 (3h steps)
        for var in _ACCUM_VARS:
            if var in ds1:
                ds1[var] = _deaccumulate(ds1, var, 3.0)

        # fc2 accumulations continue from fc1 end
        for var in _ACCUM_VARS:
            if var in ds2 and var in fc1_last:
                ds2[var] = ds2[var] - fc1_last[var]
                ds2[var] = _deaccumulate(ds2, var, 6.0)

        # Interpolate fc2 (6h) to 3h resolution
        ds2_3h = _interpolate_time(ds2, 3)

        # Concatenate fc1 + fc2
        ds_cat = xr.concat([ds1, ds2_3h], dim="valid_time")

        # Drop duplicate times
        _, unique_idx = np.unique(ds_cat.valid_time.values, return_index=True)
        ds_cat = ds_cat.isel(valid_time=sorted(unique_idx))

        # Rename variables to ERA5 convention
        rename = {}
        for old, new in IFS_RENAME_MAP.items():
            if old in ds_cat or old in ds_cat.dims:
                rename[old] = new
        if rename:
            ds_cat = ds_cat.rename(rename)

        # Compute surface geopotential from msl, sp, t2m
        if "msl" in ds_cat and "sp" in ds_cat and "t2m" in ds_cat:
            ds_cat["z"] = compute_surface_geopotential(
                ds_cat["sp"], ds_cat["msl"], ds_cat["t2m"]
            )
            ds_cat = ds_cat.drop_vars("msl")

        # Drop height dimension if present
        if "height" in ds_cat.dims:
            ds_cat = ds_cat.squeeze("height", drop=True)
        for coord in list(ds_cat.coords):
            if coord == "height":
                ds_cat = ds_cat.drop_vars("height")

        # Interpolate to target timestep
        if self._target_hours < 3:
            ds_cat = _interpolate_time(ds_cat, self._target_hours)

        # Re-accumulate to match output timestep convention
        for var in _ACCUM_VARS:
            if var in ds_cat:
                ds_cat[var] = _reaccumulate(ds_cat[var], self._target_hours)

        return ds_cat

    def _process_pressure(
        self, fc1_path: Path, fc2_path: Path
    ) -> xr.Dataset:
        """Process pressure-level GRIB files."""
        ds1 = _spatial_subset(_open_grib_pressure(str(fc1_path)), self.bbox)
        ds2 = _spatial_subset(_open_grib_pressure(str(fc2_path)), self.bbox)

        # Rename time -> valid_time
        if "time" in ds1.dims:
            ds1 = ds1.rename({"time": "valid_time"})
        if "time" in ds2.dims:
            ds2 = ds2.rename({"time": "valid_time"})

        # Rename dimensions
        for ds_ref in [ds1, ds2]:
            rename = {k: v for k, v in IFS_RENAME_MAP.items() if k in ds_ref or k in ds_ref.dims}
            if rename:
                ds_ref = ds_ref.rename(rename)

        # Convert geopotential height to geopotential (gh * g)
        for ds_ref in [ds1, ds2]:
            if "gh" in ds_ref:
                ds_ref["z"] = ds_ref["gh"] * 9.81
                ds_ref["z"].attrs = {
                    "long_name": "Geopotential",
                    "units": "m**2 s**-2",
                }

        if "gh" in ds1:
            ds1 = ds1.drop_vars("gh")
        if "gh" in ds2:
            ds2 = ds2.drop_vars("gh")

        # Convert level from Pa to hPa if needed
        for ds_ref in [ds1, ds2]:
            if "level" in ds_ref.dims:
                levels = ds_ref.level.values
                if levels.max() > 2000:  # Likely in Pa
                    ds_ref = ds_ref.assign_coords(level=levels / 100.0)

        # Interpolate fc2 (6h) to 3h
        ds2_3h = _interpolate_time(ds2, 3)

        # Concatenate
        ds_cat = xr.concat([ds1, ds2_3h], dim="valid_time")

        # Drop duplicate times
        _, unique_idx = np.unique(ds_cat.valid_time.values, return_index=True)
        ds_cat = ds_cat.isel(valid_time=sorted(unique_idx))

        # Sort levels ascending
        if "level" in ds_cat.dims:
            ds_cat = ds_cat.sortby("level", ascending=True)

        # Interpolate to target timestep
        if self._target_hours < 3:
            ds_cat = _interpolate_time(ds_cat, self._target_hours)

        return ds_cat

    def fetch_forecast(
        self, init_time: datetime
    ) -> tuple[xr.Dataset, xr.Dataset]:
        """Fetch a single IFS forecast run.

        Args:
            init_time: Forecast initialization time.

        Returns:
            Tuple of (ds_surf, ds_plev) with valid_time dimension.
        """
        tmp_dir = Path(tempfile.mkdtemp(prefix="ifs_"))

        try:
            self._download_grib(init_time, tmp_dir)

            logger.info("Processing surface forecast")
            ds_surf = self._process_surface(
                tmp_dir / "SURF_fc1.grib2",
                tmp_dir / "SURF_fc2.grib2",
            )

            logger.info("Processing pressure forecast")
            ds_plev = self._process_pressure(
                tmp_dir / "PLEV_fc1.grib2",
                tmp_dir / "PLEV_fc2.grib2",
            )

            # Add init_time as coordinate
            ds_surf = ds_surf.assign_coords(init_time=init_time)
            ds_plev = ds_plev.assign_coords(init_time=init_time)

            return ds_surf, ds_plev

        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)
