"""Hybrid backend — fetches from multiple sources per variable group.

The hybrid backend optimizes for both speed and completeness by routing
different variable groups to different backends:

- Surface variables (t2m, d2m, sp, ssrd, u10, v10): OpenMeteo S3
- Precipitation (tp): IFS 9km (2022+) or OpenMeteo S3 (pre-2022)
- Pressure levels + longwave (strd): Google ARCO

Special case: If bbox is within Central Asia and time ≤ 2023, uses s3zarr
for everything (fastest + complete).

See docs/backend-priority.md for the full priority logic.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd
import xarray as xr

from ecmwf_downloader.base import ERA5Backend
from ecmwf_downloader.bbox import BBox

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)

# ── Regional coverage ────────────────────────────────────────────────────────

# s3zarr Central Asia coverage
S3ZARR_BBOX = BBox.from_tuple((43.0, 24.0, 90.0, 58.0))
S3ZARR_END_DATE = pd.Timestamp("2023-12-31")

# ── Variable groups ──────────────────────────────────────────────────────────

SURFACE_VARS = {"t2m", "d2m", "sp", "ssrd", "u10", "v10", "msl"}
PRECIP_VARS = {"tp"}
PLEV_VARS = {"t", "z", "u", "v", "q", "r"}
STRD_VARS = {"strd"}

# Variables that must come from Google/CDS (OpenMeteo blacklisted)
GOOGLE_ONLY_VARS = {"strd", "q"}  # strd missing, q missing (only r available)
GOOGLE_ONLY_PLEV = {"u", "v"}    # 43m/s and 28m/s errors on OpenMeteo


def _bbox_within(inner: BBox, outer: BBox) -> bool:
    """Check if inner bbox is fully contained within outer bbox."""
    return (
        inner.west >= outer.west
        and inner.east <= outer.east
        and inner.south >= outer.south
        and inner.north <= outer.north
    )


class HybridBackend(ERA5Backend):
    """Hybrid backend fetching from multiple sources per variable group.

    Routes different variables to optimal backends:
    - Surface: OpenMeteo S3 (fast, exact grid)
    - Precip: IFS 9km (2022+) or S3 (pre-2022)
    - Pressure levels + strd: Google ARCO

    Args:
        bbox: Bounding box (west, south, east, north).
        pressure_levels: Pressure levels to fetch (hPa).
        time_resolution: Time resolution ('1H', '3H', '6H').
        cache_dir: Cache directory for backends.
        include_strd: Whether to fetch longwave radiation (requires Google).
        **kwargs: Additional arguments passed to sub-backends.

    Note:
        If bbox is within Central Asia and end_date ≤ 2023, automatically
        uses s3zarr for all variables (fastest + complete).
    """

    def __init__(
        self,
        bbox: BBox,
        pressure_levels: list[int] | None = None,
        time_resolution: str = "1H",
        cache_dir: str = "/tmp/hybrid_cache",
        include_strd: bool = True,
        start_date: str | pd.Timestamp | None = None,
        end_date: str | pd.Timestamp | None = None,
        **kwargs,
    ):
        super().__init__(bbox, pressure_levels or [], time_resolution, **kwargs)

        self.cache_dir = cache_dir
        self.include_strd = include_strd
        self._start_date = pd.Timestamp(start_date) if start_date else None
        self._end_date = pd.Timestamp(end_date) if end_date else None

        # Determine mode: regional (s3zarr) vs hybrid
        self._use_regional = self._check_regional_coverage()

        # Initialize sub-backends lazily
        self._surface_backend: ERA5Backend | None = None
        self._precip_backend: ERA5Backend | None = None
        self._plev_backend: ERA5Backend | None = None
        self._regional_backend: ERA5Backend | None = None

        logger.info(
            "HybridBackend: bbox=%s, plev=%s, regional=%s",
            (bbox.west, bbox.south, bbox.east, bbox.north),
            pressure_levels,
            self._use_regional,
        )

    def _check_regional_coverage(self) -> bool:
        """Check if we can use s3zarr for everything."""
        if not _bbox_within(self.bbox, S3ZARR_BBOX):
            return False
        if self._end_date and self._end_date > S3ZARR_END_DATE:
            return False
        # Check if s3zarr is available
        try:
            from ecmwf_downloader.reanalysis.backends import _check_available
            return _check_available("s3zarr")
        except ImportError:
            return False

    def _get_surface_backend(self) -> ERA5Backend:
        """Get or create the surface variables backend."""
        if self._surface_backend is None:
            from ecmwf_downloader.reanalysis.backends import get_backend, _check_available

            # Try backends in order
            for name in ["openmeteo_s3", "openmeteo", "google", "cds"]:
                if _check_available(name):
                    try:
                        kwargs = {
                            "bbox": self.bbox,
                            "cache_dir": self.cache_dir,
                            "pressure_levels": [],  # Surface only
                        }
                        if name == "openmeteo":
                            kwargs["model"] = "era5"
                            # OpenMeteo needs date range
                            if self._start_date and self._end_date:
                                kwargs["start_date"] = str(self._start_date.date())
                                kwargs["end_date"] = str(self._end_date.date())
                        if name == "openmeteo_s3":
                            # S3 backend doesn't need pressure_levels in kwargs
                            del kwargs["pressure_levels"]
                        self._surface_backend = get_backend(name, **kwargs)
                        logger.info("HybridBackend: surface backend = %s", name)
                        break
                    except Exception as e:
                        logger.warning("Failed to init %s for surface: %s", name, e)

            if self._surface_backend is None:
                raise RuntimeError("No surface backend available")

        return self._surface_backend

    def _get_precip_backend(self, date: pd.Timestamp) -> ERA5Backend:
        """Get or create the precipitation backend."""
        # For 2022+, prefer IFS 9km
        use_ifs = date >= pd.Timestamp("2022-01-01")

        if self._precip_backend is None or (use_ifs and "ifs" not in str(type(self._precip_backend)).lower()):
            from ecmwf_downloader.reanalysis.backends import get_backend, _check_available

            if use_ifs:
                candidates = ["openmeteo", "openmeteo_s3", "google", "cds"]
            else:
                candidates = ["openmeteo_s3", "openmeteo", "google", "cds"]

            for name in candidates:
                if _check_available(name):
                    try:
                        kwargs = {
                            "bbox": self.bbox,
                            "cache_dir": self.cache_dir,
                            "pressure_levels": [],
                        }
                        if name == "openmeteo" and use_ifs:
                            kwargs["model"] = "ifs"  # 9km precipitation
                            # OpenMeteo needs date range
                            if self._start_date and self._end_date:
                                kwargs["start_date"] = str(self._start_date.date())
                                kwargs["end_date"] = str(self._end_date.date())
                        elif name == "openmeteo":
                            kwargs["model"] = "era5"
                            if self._start_date and self._end_date:
                                kwargs["start_date"] = str(self._start_date.date())
                                kwargs["end_date"] = str(self._end_date.date())
                        if name == "openmeteo_s3":
                            del kwargs["pressure_levels"]
                        self._precip_backend = get_backend(name, **kwargs)
                        logger.info("HybridBackend: precip backend = %s (ifs=%s)", name, use_ifs)
                        break
                    except Exception as e:
                        logger.warning("Failed to init %s for precip: %s", name, e)

            if self._precip_backend is None:
                raise RuntimeError("No precip backend available")

        return self._precip_backend

    def _get_plev_backend(self) -> ERA5Backend:
        """Get or create the pressure level + strd backend."""
        if self._plev_backend is None:
            from ecmwf_downloader.reanalysis.backends import get_backend, _check_available

            # Only Google/CDS have complete plev + strd
            # Use explicit variable lists to avoid relative_humidity issue
            plev_vars = [
                "temperature",
                "geopotential",
                "u_component_of_wind",
                "v_component_of_wind",
                "specific_humidity",
            ]

            for name in ["google", "cds"]:
                if _check_available(name):
                    try:
                        self._plev_backend = get_backend(
                            name,
                            bbox=self.bbox,
                            pressure_levels=self.pressure_levels,
                            cache_dir=self.cache_dir,
                            plev_vars=plev_vars,
                        )
                        logger.info("HybridBackend: plev backend = %s", name)
                        break
                    except Exception as e:
                        logger.warning("Failed to init %s for plev: %s", name, e)

            if self._plev_backend is None:
                raise RuntimeError("No plev backend available (need Google or CDS)")

        return self._plev_backend

    def _get_regional_backend(self) -> ERA5Backend:
        """Get or create the regional (s3zarr) backend."""
        if self._regional_backend is None:
            from ecmwf_downloader.reanalysis.backends import get_backend

            self._regional_backend = get_backend(
                "s3zarr",
                bbox=self.bbox,
                pressure_levels=self.pressure_levels,
                cache_dir=self.cache_dir,
            )
            logger.info("HybridBackend: using regional s3zarr backend")

        return self._regional_backend

    def fetch_day(self, date: pd.Timestamp) -> tuple[xr.Dataset, xr.Dataset]:
        """Fetch one day of ERA5 data from multiple backends.

        Args:
            date: The date to fetch.

        Returns:
            Tuple of (ds_surf, ds_plev).
        """
        date = pd.Timestamp(date).normalize()

        # ─── Regional mode: s3zarr for everything ────────────────────────
        if self._use_regional:
            backend = self._get_regional_backend()
            return backend.fetch_day(date)

        # ─── Hybrid mode: multiple backends ──────────────────────────────

        # 1. Surface variables from OpenMeteo S3
        surface_backend = self._get_surface_backend()
        ds_surf_base, _ = surface_backend.fetch_day(date)

        # 2. Precipitation (potentially from different backend for 9km)
        precip_backend = self._get_precip_backend(date)
        if precip_backend is not surface_backend:
            ds_precip, _ = precip_backend.fetch_day(date)
            # Extract just tp and merge
            if "tp" in ds_precip:
                ds_surf_base["tp"] = ds_precip["tp"]
        # else: tp already in ds_surf_base

        # 3. Pressure levels + strd from Google
        ds_plev = xr.Dataset()
        if self.pressure_levels or self.include_strd:
            plev_backend = self._get_plev_backend()
            ds_surf_google, ds_plev = plev_backend.fetch_day(date)

            # Add strd from Google to surface dataset
            if self.include_strd and "strd" in ds_surf_google:
                # Align coordinates before merging
                strd = ds_surf_google["strd"]
                # Interpolate to match surface backend grid if needed
                if not _coords_match(ds_surf_base, ds_surf_google):
                    strd = strd.interp(
                        latitude=ds_surf_base.latitude,
                        longitude=ds_surf_base.longitude,
                        method="linear",
                    )
                ds_surf_base["strd"] = strd

        logger.info(
            "HybridBackend: fetched %s — SURF vars=%s, PLEV vars=%s",
            date.strftime("%Y-%m-%d"),
            list(ds_surf_base.data_vars),
            list(ds_plev.data_vars) if ds_plev else [],
        )

        return ds_surf_base, ds_plev

    def close(self):
        """Clean up all sub-backends."""
        for backend in [
            self._surface_backend,
            self._precip_backend,
            self._plev_backend,
            self._regional_backend,
        ]:
            if backend is not None:
                try:
                    backend.close()
                except Exception:
                    pass


def _coords_match(ds1: xr.Dataset, ds2: xr.Dataset, tol: float = 0.01) -> bool:
    """Check if two datasets have matching lat/lon coordinates."""
    try:
        lat1 = ds1.latitude.values
        lat2 = ds2.latitude.values
        lon1 = ds1.longitude.values
        lon2 = ds2.longitude.values

        if len(lat1) != len(lat2) or len(lon1) != len(lon2):
            return False

        return (
            abs(lat1 - lat2).max() < tol
            and abs(lon1 - lon2).max() < tol
        )
    except Exception:
        return False
