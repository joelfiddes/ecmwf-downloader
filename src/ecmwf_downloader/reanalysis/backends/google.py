"""Google ARCO-ERA5 backend.

Downloads ERA5 data from Google Cloud Storage's public ARCO-ERA5 dataset.
Anonymous access, no credentials required. Requires ``gcsfs`` and ``h5netcdf``.

URI pattern:
    gs://gcp-public-data-arco-era5/raw/date-variable-{single_level|pressure_level}/
        {YYYY}/{MM}/{DD}/{variable}/{level|surface}.nc

Ported from: era5google/era5_downloader/core.py (Luke Gregor)
"""

from __future__ import annotations

import json
import logging
import pathlib
import tempfile
from typing import Iterator

import fsspec
import pandas as pd
import xarray as xr
from tenacity import retry, stop_after_attempt, wait_exponential
from tqdm import tqdm

from ecmwf_downloader.base import ERA5Backend
from ecmwf_downloader.bbox import BBox
from ecmwf_downloader.variables import (
    DEFAULT_PLEV_VARS,
    DEFAULT_SURF_VARS,
    PLEV_VARS_GOOGLE,
    SURF_VARS_GOOGLE,
    TIME_RESOLUTION_HOURS,
    VARS_BOTH_SURFACE_AND_LEVEL,
)

logger = logging.getLogger(__name__)

URI_LEVELS = (
    "filecache::gs://gcp-public-data-arco-era5/raw/"
    "date-variable-pressure_level/{t:%Y}/{t:%m}/{t:%d}/{variable}/{level}.nc"
)
URI_SURFACE = (
    "filecache::gs://gcp-public-data-arco-era5/raw/"
    "date-variable-single_level/{t:%Y}/{t:%m}/{t:%d}/{variable}/surface.nc"
)


class GoogleCloudBackend(ERA5Backend):
    """ERA5 backend using Google ARCO-ERA5 public netCDF files.

    Args:
        bbox: Bounding box (west, south, east, north) in -180:180 coords.
        pressure_levels: List of pressure levels in hPa.
        time_resolution: '1H', '2H', '3H', or '6H'.
        cache_dir: Directory for fsspec file cache. Cleaned after each day.
        surf_vars: Surface variable long names to download.
        plev_vars: Pressure level variable long names to download.
    """

    def __init__(
        self,
        bbox: BBox,
        pressure_levels: list[int],
        time_resolution: str = "1H",
        cache_dir: str = "./era5_cache/",
        surf_vars: list[str] | None = None,
        plev_vars: list[str] | None = None,
        show_progress: bool = True,
        **kwargs,
    ):
        super().__init__(bbox, pressure_levels, time_resolution, **kwargs)

        # Convert bbox to 0:360 for Google's data
        self._bbox_360 = self.bbox.to_0_360()

        self.surf_vars = surf_vars or DEFAULT_SURF_VARS
        self.plev_vars = plev_vars or DEFAULT_PLEV_VARS
        self.show_progress = show_progress

        self._time_steps = TIME_RESOLUTION_HOURS.get(time_resolution)
        if self._time_steps is None:
            raise ValueError(
                f"Invalid time_resolution '{time_resolution}'. "
                f"Must be one of: {list(TIME_RESOLUTION_HOURS.keys())}"
            )

        # Set up cache directory
        self._cache_dir = pathlib.Path(cache_dir).expanduser().resolve()
        self._cache_dir.mkdir(parents=True, exist_ok=True)

        # Discover valid variables from the store
        self._valid_surface_vars: list[str] | None = None
        self._valid_level_vars: list[str] | None = None
        self._valid_levels: list[int] | None = None

    def _ensure_valid_vars(self):
        """Lazily discover available variables from the Google store."""
        if self._valid_surface_vars is not None:
            return

        fs = fsspec.filesystem("gs", token="anon")
        t = pd.Timestamp("2000-01-01")

        base_surface = (
            "gcp-public-data-arco-era5/raw/date-variable-single_level/"
            f"{t:%Y}/{t:%m}/{t:%d}"
        )
        base_levels = (
            "gcp-public-data-arco-era5/raw/date-variable-pressure_level/"
            f"{t:%Y}/{t:%m}/{t:%d}"
        )

        def get_name(f):
            return f.split("/")[-1].replace(".nc", "")

        self._valid_surface_vars = [get_name(f) for f in fs.ls(base_surface)]
        self._valid_level_vars = [get_name(f) for f in fs.ls(base_levels)]

        level_path = f"{base_levels}/temperature"
        self._valid_levels = sorted([int(get_name(f)) for f in fs.ls(level_path)])

        logger.debug("Valid surface vars: %s", self._valid_surface_vars)
        logger.debug("Valid level vars: %s", self._valid_level_vars)
        logger.debug("Valid levels: %s", self._valid_levels)

    def _surface_or_level(self, variable: str) -> str:
        """Classify a variable as 'surface', 'level', or 'both'."""
        self._ensure_valid_vars()
        in_surf = variable in self._valid_surface_vars
        in_level = variable in self._valid_level_vars
        if in_surf and in_level:
            return "both"
        if in_surf:
            return "surface"
        if in_level:
            return "level"
        raise ValueError(
            f"Variable '{variable}' not found in Google ARCO-ERA5. "
            f"Surface: {self._valid_surface_vars}, Level: {self._valid_level_vars}"
        )

    def _make_uri_list(self, t: pd.Timestamp) -> Iterator[str]:
        """Generate URIs for all variables and levels for a given day."""
        all_vars = set(self.surf_vars) | set(self.plev_vars)
        for variable in all_vars:
            var_type = self._surface_or_level(variable)
            if var_type in ("surface", "both"):
                yield URI_SURFACE.format(t=t, variable=variable)
            if var_type in ("level", "both"):
                for level in self.pressure_levels:
                    yield URI_LEVELS.format(t=t, variable=variable, level=level)

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10), reraise=True)
    def _download_files(self, uri_list: list[str]) -> tuple[list[str], pathlib.Path, dict]:
        """Download files using fsspec with filecache.

        Returns:
            Tuple of (file_list, tmp_dir, fsspec_cache) - all local to this call,
            not stored as instance variables to avoid race conditions in parallel downloads.
        """
        tmp_dir = pathlib.Path(tempfile.mkdtemp(dir=str(self._cache_dir)))
        logger.debug("Downloading %d files to %s", len(uri_list), tmp_dir)

        flist = fsspec.open_local(
            url=uri_list,
            filecache=dict(cache_storage=str(tmp_dir)),
            gs=dict(token="anon"),
        )

        # Read fsspec cache for URI mapping
        cache_file = tmp_dir / "cache"
        if cache_file.exists():
            with open(cache_file) as f:
                cache = json.load(f)
            fsspec_cache = {cache[k]["fn"]: cache[k] for k in cache}
        else:
            fsspec_cache = {}

        return list(flist), tmp_dir, fsspec_cache

    def _get_original_uri(self, ds: xr.Dataset, fsspec_cache: dict) -> str:
        """Look up the original URI from fsspec cache.

        Args:
            ds: Dataset with encoding["source"] pointing to cached file.
            fsspec_cache: Cache dict mapping file hashes to metadata (local to fetch_day call).
        """
        fname_hash = pathlib.Path(ds.encoding["source"]).name
        entry = fsspec_cache.get(fname_hash, {})
        return entry.get("original", "")

    def _preprocess(self, ds: xr.Dataset, fsspec_cache: dict) -> xr.Dataset:
        """Preprocess a single file: subset, rename, add level dim.

        Args:
            ds: Raw dataset from a single netCDF file.
            fsspec_cache: Cache dict for URI lookup (local to fetch_day call).
        """
        uri = self._get_original_uri(ds, fsspec_cache)
        is_level = "pressure_level" in uri

        # Variable name from URI path
        parts = uri.split("/")
        if len(parts) >= 2:
            request_variable = parts[-2]
        else:
            request_variable = ""

        var_type = "level" if is_level else "surface"
        if request_variable in VARS_BOTH_SURFACE_AND_LEVEL:
            var_type = "level" if is_level else "both_surface"

        # Spatial + temporal subset
        bb = self._bbox_360
        ds = ds.sel(
            latitude=slice(bb.north, bb.south),
            longitude=slice(bb.west, bb.east),
            time=ds.time.dt.hour.isin(self._time_steps),
        )

        # Get the data variable key
        key = list(ds.data_vars)[0]

        # Rename surface geopotential to avoid clash with pressure-level geopotential
        if var_type == "both_surface":
            ds = ds.rename({key: key + "_surf"})
            key = key + "_surf"

        # Add level dimension for pressure-level files
        if is_level:
            level_str = uri.split("/")[-1].replace(".nc", "")
            if level_str.isdigit():
                ds = ds.expand_dims(level=[int(level_str)])

        ds = ds.astype("float32")
        return ds

    def _clear_cache(self, tmp_dir: pathlib.Path | None = None):
        """Remove temporary cache files.

        Args:
            tmp_dir: Temporary directory to clean up. If None, does nothing.
                     Each fetch_day call passes its own tmp_dir to avoid race conditions.
        """
        if tmp_dir is None or not tmp_dir.exists():
            return

        for f in tmp_dir.iterdir():
            if f.is_file():
                try:
                    f.unlink()
                except FileNotFoundError:
                    pass
        try:
            tmp_dir.rmdir()
        except OSError:
            pass

    def fetch_day(self, date: pd.Timestamp) -> tuple[xr.Dataset, xr.Dataset]:
        """Fetch one day of ERA5 data from Google ARCO-ERA5.

        This method is thread-safe: each call uses its own temporary directory
        and cache, avoiding race conditions when downloading multiple days in parallel.

        Args:
            date: The date to fetch.

        Returns:
            Tuple of (ds_surf, ds_plev) with standardised names and dims.
        """
        date = pd.Timestamp(date)
        logger.info("Fetching Google ARCO-ERA5 data for %s", date.strftime("%Y-%m-%d"))

        uri_list = list(self._make_uri_list(date))
        flist, tmp_dir, fsspec_cache = self._download_files(uri_list)

        logger.debug("Opening and preprocessing %d files", len(flist))
        # Open each file individually, preprocess, load into memory, then merge
        # (open_mfdataset with combine_by_coords fails when files have
        # different variables but same coordinates - they need merge, not concat)
        datasets = []
        iterator = flist
        if self.show_progress:
            iterator = tqdm(
                flist,
                desc=f"Google {date.strftime('%Y-%m-%d')}",
                unit="file",
                leave=False,
            )
        for f in iterator:
            ds_single = xr.open_dataset(f, engine="scipy")
            ds_single = self._preprocess(ds_single, fsspec_cache).load()  # Load immediately
            datasets.append(ds_single)

        # Merge all datasets (handles different variables with same coords)
        ds = xr.merge(datasets, compat="override", join="outer")

        # Clean up this call's temporary directory (thread-safe: each call has its own)
        self._clear_cache(tmp_dir)

        # Split into surface and pressure-level datasets
        surf_vars = []
        plev_vars = []
        for var in ds.data_vars:
            if "level" in ds[var].dims:
                plev_vars.append(var)
            else:
                surf_vars.append(var)

        ds_surf = ds[surf_vars] if surf_vars else xr.Dataset()
        ds_plev = ds[plev_vars] if plev_vars else xr.Dataset()

        # Sort pressure levels ascending
        if "level" in ds_plev.dims:
            ds_plev = ds_plev.sortby("level", ascending=True)

        # Rename short names from Google's variable names
        surf_rename = {}
        for long_name, short_name in SURF_VARS_GOOGLE.items():
            # Surface geopotential was renamed to z_surf in preprocess
            if long_name in VARS_BOTH_SURFACE_AND_LEVEL:
                if short_name + "_surf" in ds_surf:
                    surf_rename[short_name + "_surf"] = short_name
            # Other variables keep their short name from the netCDF
            # (Google files already use short names as data var keys)

        if surf_rename:
            ds_surf = ds_surf.rename(surf_rename)

        logger.info(
            "Fetched: SURF vars=%s, PLEV vars=%s, levels=%s",
            list(ds_surf.data_vars),
            list(ds_plev.data_vars),
            list(ds_plev.level.values) if "level" in ds_plev.dims else [],
        )

        return ds_surf, ds_plev

    def close(self):
        """Clean up any remaining cache.

        Note: With the thread-safe design, each fetch_day() cleans its own cache,
        so this is typically a no-op. Kept for interface compatibility.
        """
        pass
