"""IFS Forecast Loader.

Downloads and processes ECMWF IFS open data forecasts using cfgrib + xarray.
No CDO dependency. Outputs match ERA5 conventions (accumulated radiation/precip
per timestep, standard variable names and dimensions).

Ported and redesigned from: snowmapper/fetch_ifs_forecast.py

IFS open data structure:
    - For times 00z & 12z: steps 0-144 by 3h, steps 150-240 by 6h
    - Surface params: 2t, sp, 2d, ssrd, strd, tp, msl
    - Pressure params: gh, u, v, r, q, t
    - Default levels: 1000, 925, 850, 700, 600, 500, 400, 300

Output files:
    {output_dir}/SURF_FC.nc     — full merged forecast
    {output_dir}/PLEV_FC.nc     — full merged forecast
    {output_dir}/SURF_FC_YYYY-MM-DD.nc  — daily hindcast
    {output_dir}/PLEV_FC_YYYY-MM-DD.nc  — daily hindcast
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

from ecmwf_downloader.bbox import BBox
from ecmwf_downloader.derived import compute_surface_geopotential
from ecmwf_downloader.variables import (
    IFS_DEFAULT_LEVELS,
    IFS_PLEV_PARAMS,
    IFS_RENAME_MAP,
    IFS_SURF_PARAMS,
)

logger = logging.getLogger(__name__)

# Accumulated variables that need deaccumulation/re-accumulation
_ACCUM_VARS = {"ssrd", "strd", "tp"}

# Output timestep in hours for interpolation targets
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
    """Subset dataset to bounding box.

    Handles both 'latitude'/'longitude' and 'lat'/'lon' dimension names.
    """
    lat_name = "latitude" if "latitude" in ds.dims else "lat"
    lon_name = "longitude" if "longitude" in ds.dims else "lon"

    lat = ds[lat_name].values
    lon = ds[lon_name].values

    lat_mask = (lat >= bbox.south) & (lat <= bbox.north)
    lon_mask = (lon >= bbox.west) & (lon <= bbox.east)

    return ds.sel({lat_name: lat[lat_mask], lon_name: lon[lon_mask]})


def _deaccumulate(ds: xr.Dataset, var_name: str, divisor: float) -> xr.DataArray:
    """Deaccumulate a forecast variable: diff between consecutive steps.

    Returns instantaneous rate (value / divisor) per timestep.
    """
    acc = ds[var_name]
    deacc = acc - acc.shift(time=1, fill_value=0)
    return deacc / divisor


def _interpolate_time(ds: xr.Dataset, target_freq_hours: int) -> xr.Dataset:
    """Interpolate dataset to a finer time resolution using xr.interp.

    Args:
        ds: Dataset with 'time' dimension.
        target_freq_hours: Target timestep in hours (1, 2, or 3).

    Returns:
        Dataset interpolated to the target frequency.
    """
    times = pd.DatetimeIndex(ds.time.values)
    new_times = pd.date_range(times[0], times[-1], freq=f"{target_freq_hours}h")
    return ds.interp(time=new_times, method="linear")


def _reaccumulate(da: xr.DataArray, source_hours: float, target_hours: float) -> xr.DataArray:
    """Re-accumulate a rate variable to match target timestep accumulation.

    Converts from instantaneous rate (value/source_hours) to accumulated
    per target_hours: rate * source_hours gives back original accumulation,
    then we need to express it as accumulation per target timestep.

    For interpolated data: the interpolated rate * target_hours gives the
    accumulation per target timestep.
    """
    return da * target_hours


class IFSForecastLoader:
    """Download and process ECMWF IFS open data forecasts.

    Args:
        bbox: Bounding box as (W, S, E, N) tuple or BBox instance.
        output_dir: Directory for output forecast files.
        output_timestep: Output time resolution: '1H', '2H', or '3H'.
        forecast_time: Forecast base time: 0 or 12 (UTC).
        pressure_levels: Pressure levels to extract.

    Example::

        fc = IFSForecastLoader(
            bbox=(59, 32, 81, 45),
            output_dir="./inputs/climate/forecast/",
            output_timestep="1H",
        )
        fc.download(backfill=True)
    """

    def __init__(
        self,
        bbox: tuple | BBox,
        output_dir: str = "./inputs/climate/forecast/",
        output_timestep: str = "1H",
        forecast_time: int = 0,
        pressure_levels: list[int] | None = None,
    ):
        if isinstance(bbox, (tuple, list)):
            self.bbox = BBox.from_tuple(tuple(bbox))
        else:
            self.bbox = bbox

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if output_timestep not in _TIMESTEP_HOURS:
            raise ValueError(
                f"Invalid output_timestep '{output_timestep}'. "
                f"Must be one of: {list(_TIMESTEP_HOURS.keys())}"
            )
        self.output_timestep = output_timestep
        self._target_hours = _TIMESTEP_HOURS[output_timestep]

        self.forecast_time = forecast_time
        self.pressure_levels = pressure_levels or IFS_DEFAULT_LEVELS

    def _get_existing_forecast_dates(self) -> set[datetime]:
        """Find existing daily hindcast files."""
        dates = set()
        for f in self.output_dir.glob("SURF_FC_20*.nc"):
            name = f.stem  # e.g. SURF_FC_2026-01-14
            if len(name) >= 18:
                date_str = name[8:18]
                try:
                    dates.add(datetime.strptime(date_str, "%Y-%m-%d"))
                except ValueError:
                    pass
        return dates

    def _get_missing_dates(self) -> list[tuple[datetime, int]]:
        """Identify missing forecast gap-fill dates.

        IFS open data keeps ~3 days of historical forecasts.

        Returns:
            List of (date, mydate_offset) tuples.
        """
        existing = self._get_existing_forecast_dates()
        today = datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)

        missing = []
        for days_ago in range(3, 0, -1):
            date = today - timedelta(days=days_ago)
            if date not in existing:
                missing.append((date, -days_ago))

        return missing

    def _download_grib(self, mydate: int, tmp_dir: Path) -> None:
        """Download 4 GRIB files (surf fc1, plev fc1, surf fc2, plev fc2).

        Args:
            mydate: Date offset (0=today, -1=yesterday, etc.)
            tmp_dir: Temporary directory for GRIB files.
        """
        from ecmwf.opendata import Client

        client = Client()

        fc1_steps = list(range(0, 147, 3))  # 0 to 144 by 3h
        fc2_steps = list(range(150, 241, 6))  # 150 to 240 by 6h

        logger.info("Downloading IFS forecast GRIB (mydate=%d)", mydate)

        # Surface fc1
        client.retrieve(
            time=self.forecast_time,
            date=mydate,
            step=fc1_steps,
            type="fc",
            param=IFS_SURF_PARAMS,
            target=str(tmp_dir / "SURF_fc1.grib2"),
        )

        # Pressure levels fc1
        client.retrieve(
            time=self.forecast_time,
            date=mydate,
            step=fc1_steps,
            type="fc",
            param=IFS_PLEV_PARAMS,
            levelist=self.pressure_levels,
            target=str(tmp_dir / "PLEV_fc1.grib2"),
        )

        # Surface fc2
        client.retrieve(
            time=self.forecast_time,
            date=mydate,
            step=fc2_steps,
            type="fc",
            param=IFS_SURF_PARAMS,
            target=str(tmp_dir / "SURF_fc2.grib2"),
        )

        # Pressure levels fc2
        client.retrieve(
            time=self.forecast_time,
            date=mydate,
            step=fc2_steps,
            type="fc",
            param=IFS_PLEV_PARAMS,
            levelist=self.pressure_levels,
            target=str(tmp_dir / "PLEV_fc2.grib2"),
        )

    def _process_surface(
        self, fc1_path: Path, fc2_path: Path
    ) -> xr.Dataset:
        """Process surface GRIB files: open, subset, deaccumulate, rename.

        Handles accumulation continuity at fc1/fc2 boundary.
        """
        ds1 = _spatial_subset(_open_grib_surface(str(fc1_path)), self.bbox)
        ds2 = _spatial_subset(_open_grib_surface(str(fc2_path)), self.bbox)

        # Handle tp variable name quirk
        for ds in [ds1, ds2]:
            if "tp" not in ds and "param193.1.0" in ds:
                ds["tp"] = ds["param193.1.0"]
                ds = ds.drop_vars("param193.1.0", errors="ignore")

        # Save last values of fc1 accumulators for fc2 boundary correction
        fc1_last = {}
        for var in _ACCUM_VARS:
            if var in ds1:
                fc1_last[var] = ds1[var].isel(time=-1)

        # Deaccumulate fc1 (3h steps)
        for var in _ACCUM_VARS:
            if var in ds1:
                ds1[var] = _deaccumulate(ds1, var, 3.0)

        # fc2 accumulations continue from fc1 end — subtract fc1 last value
        for var in _ACCUM_VARS:
            if var in ds2 and var in fc1_last:
                ds2[var] = ds2[var] - fc1_last[var]
                ds2[var] = _deaccumulate(ds2, var, 6.0)

        # Interpolate fc2 (6h) to 3h resolution
        ds2_3h = _interpolate_time(ds2, 3)

        # Concatenate fc1 + fc2
        ds_cat = xr.concat([ds1, ds2_3h], dim="time")

        # Drop duplicate times
        _, unique_idx = np.unique(ds_cat.time.values, return_index=True)
        ds_cat = ds_cat.isel(time=sorted(unique_idx))

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

        # Re-accumulate radiation and precip to match ERA5 convention
        # (J/m² and m accumulated per output timestep)
        for var in _ACCUM_VARS:
            if var in ds_cat:
                ds_cat[var] = _reaccumulate(
                    ds_cat[var], source_hours=3.0, target_hours=self._target_hours
                )

        return ds_cat

    def _process_pressure(
        self, fc1_path: Path, fc2_path: Path
    ) -> xr.Dataset:
        """Process pressure-level GRIB files: open, subset, rename."""
        ds1 = _spatial_subset(_open_grib_pressure(str(fc1_path)), self.bbox)
        ds2 = _spatial_subset(_open_grib_pressure(str(fc2_path)), self.bbox)

        # Rename dimensions
        for ds in [ds1, ds2]:
            rename = {}
            for old, new in IFS_RENAME_MAP.items():
                if old in ds or old in ds.dims:
                    rename[old] = new
            if rename:
                ds_tmp = ds.rename(rename)
                # We need to reassign since ds is in a list
                if ds is ds1:
                    ds1 = ds_tmp
                else:
                    ds2 = ds_tmp

        # Actually re-do renaming properly
        rename1 = {k: v for k, v in IFS_RENAME_MAP.items() if k in ds1 or k in ds1.dims}
        rename2 = {k: v for k, v in IFS_RENAME_MAP.items() if k in ds2 or k in ds2.dims}
        if rename1:
            ds1 = ds1.rename(rename1)
        if rename2:
            ds2 = ds2.rename(rename2)

        # Convert geopotential height to geopotential (gh * g)
        for ds_ref, label in [(ds1, "fc1"), (ds2, "fc2")]:
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
        ds_cat = xr.concat([ds1, ds2_3h], dim="time")

        # Drop duplicate times
        _, unique_idx = np.unique(ds_cat.time.values, return_index=True)
        ds_cat = ds_cat.isel(time=sorted(unique_idx))

        # Sort levels ascending
        if "level" in ds_cat.dims:
            ds_cat = ds_cat.sortby("level", ascending=True)

        # Interpolate to target timestep
        if self._target_hours < 3:
            ds_cat = _interpolate_time(ds_cat, self._target_hours)

        return ds_cat

    def _download_and_process(self, mydate: int) -> str:
        """Download and process one forecast.

        Args:
            mydate: Date offset (0=today, -1=yesterday, etc.)

        Returns:
            Date string (YYYY-MM-DD) of the forecast base date.
        """
        tmp_dir = Path(tempfile.mkdtemp(prefix="ifs_"))

        try:
            self._download_grib(mydate, tmp_dir)

            logger.info("Processing surface forecast (mydate=%d)", mydate)
            ds_surf = self._process_surface(
                tmp_dir / "SURF_fc1.grib2",
                tmp_dir / "SURF_fc2.grib2",
            )

            logger.info("Processing pressure forecast (mydate=%d)", mydate)
            ds_plev = self._process_pressure(
                tmp_dir / "PLEV_fc1.grib2",
                tmp_dir / "PLEV_fc2.grib2",
            )

            # Determine the date string from the first time step
            first_time = pd.Timestamp(ds_surf.time.values[0])
            day_str = first_time.strftime("%Y-%m-%d")

            # Write full forecast files
            ds_surf.to_netcdf(self.output_dir / "SURF_FC.nc", mode="w")
            ds_plev.to_netcdf(self.output_dir / "PLEV_FC.nc", mode="w")

            # Write daily hindcast (first 24 hours)
            n_steps = 24 // self._target_hours
            ds_surf_day = ds_surf.isel(time=slice(0, n_steps))
            ds_plev_day = ds_plev.isel(time=slice(0, n_steps))

            ds_surf_day.to_netcdf(
                self.output_dir / f"SURF_FC_{day_str}.nc", mode="w"
            )
            ds_plev_day.to_netcdf(
                self.output_dir / f"PLEV_FC_{day_str}.nc", mode="w"
            )

            logger.info("Saved forecast products for %s", day_str)
            return day_str

        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    def _merge_all_forecasts(self) -> None:
        """Merge all hindcast + forecast files into final SURF_FC.nc / PLEV_FC.nc."""
        for prefix in ["SURF", "PLEV"]:
            fc_files = sorted(self.output_dir.glob(f"{prefix}_FC_*.nc"))
            full_fc = self.output_dir / f"{prefix}_FC.nc"

            all_files = list(fc_files)
            if full_fc.exists():
                all_files.append(full_fc)

            if not all_files:
                continue

            logger.info("Merging %d %s forecast files", len(all_files), prefix)
            ds = xr.open_mfdataset(
                [str(f) for f in all_files],
                combine="by_coords",
            ).load()

            # Remove duplicate timestamps
            _, unique_idx = np.unique(ds.time.values, return_index=True)
            ds = ds.isel(time=sorted(unique_idx))

            ds.to_netcdf(self.output_dir / f"{prefix}_FC.nc", mode="w")
            ds.close()

            logger.info("Merged %s forecast saved", prefix)

    def download(self, backfill: bool = True) -> None:
        """Run the full IFS forecast download pipeline.

        Args:
            backfill: If True, download missing historical forecasts (up to 3 days).
        """
        logger.info("Starting IFS forecast download")

        if backfill:
            missing = self._get_missing_dates()
            if missing:
                logger.info("Backfilling %d missing dates", len(missing))
                for date, offset in missing:
                    logger.info(
                        "  %s (mydate=%d)", date.strftime("%Y-%m-%d"), offset
                    )
                    try:
                        self._download_and_process(offset)
                    except Exception:
                        logger.error(
                            "Failed to backfill %s",
                            date.strftime("%Y-%m-%d"),
                            exc_info=True,
                        )
            else:
                logger.info("No missing dates to backfill")

        # Download today's forecast
        today = datetime.now()
        today_surf = self.output_dir / f"SURF_FC_{today:%Y-%m-%d}.nc"
        today_plev = self.output_dir / f"PLEV_FC_{today:%Y-%m-%d}.nc"

        if today_surf.exists() and today_plev.exists():
            logger.info("Today's forecast already exists, skipping")
        else:
            logger.info("Downloading today's forecast (mydate=0)")
            self._download_and_process(0)

        # Merge all forecasts
        self._merge_all_forecasts()

        logger.info("IFS forecast download complete")
