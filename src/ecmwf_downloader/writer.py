"""Output writers for ERA5 data: NetCDF and Zarr formats."""

from __future__ import annotations

import logging
import shutil
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path

import xarray as xr

logger = logging.getLogger(__name__)


class BaseWriter(ABC):
    """Abstract writer interface."""

    @abstractmethod
    def write_day(
        self, ds_surf: xr.Dataset, ds_plev: xr.Dataset, date_str: str
    ) -> None:
        """Write one day of data."""
        ...

    @abstractmethod
    def finalize(self) -> None:
        """Finalize output (merge, consolidate, etc.)."""
        ...


class NetCDFWriter(BaseWriter):
    """Write ERA5 data as daily netCDF files, merged to yearly files.

    Output structure:
        {output_dir}/daily/dSURF_YYYYMMDD.nc
        {output_dir}/daily/dPLEV_YYYYMMDD.nc
        {output_dir}/yearly/SURF_YYYY.nc  (after finalize)
        {output_dir}/yearly/PLEV_YYYY.nc  (after finalize)

    Args:
        output_dir: Base output directory.
        complevel: zlib compression level (1-9). Default 1.
    """

    def __init__(self, output_dir: str | Path, complevel: int = 1):
        self.output_dir = Path(output_dir)
        self.complevel = complevel

        self.daily_dir = self.output_dir / "daily"
        self.yearly_dir = self.output_dir / "yearly"
        self.daily_dir.mkdir(parents=True, exist_ok=True)
        self.yearly_dir.mkdir(parents=True, exist_ok=True)

        self._written_dates: list[str] = []

    def _encoding(self, ds: xr.Dataset) -> dict:
        """Build per-variable encoding dict for netCDF compression."""
        return {
            var: {"zlib": True, "complevel": self.complevel, "dtype": "float32"}
            for var in ds.data_vars
        }

    def _atomic_write(self, ds: xr.Dataset, target: Path) -> None:
        """Write dataset atomically via temp file + rename."""
        tmp_fd, tmp_path = tempfile.mkstemp(
            suffix=".nc.tmp", dir=str(target.parent)
        )
        try:
            import os
            os.close(tmp_fd)
            ds.to_netcdf(tmp_path, encoding=self._encoding(ds), engine="netcdf4")
            shutil.move(tmp_path, str(target))
        except Exception:
            Path(tmp_path).unlink(missing_ok=True)
            raise

    def write_day(
        self, ds_surf: xr.Dataset, ds_plev: xr.Dataset, date_str: str
    ) -> None:
        """Write daily SURF and PLEV netCDF files.

        Args:
            ds_surf: Surface dataset.
            ds_plev: Pressure level dataset.
            date_str: Date string in YYYYMMDD format.
        """
        surf_path = self.daily_dir / f"dSURF_{date_str}.nc"
        plev_path = self.daily_dir / f"dPLEV_{date_str}.nc"

        if ds_surf.data_vars:
            self._atomic_write(ds_surf, surf_path)
            logger.debug("Wrote %s", surf_path.name)

        if ds_plev.data_vars:
            self._atomic_write(ds_plev, plev_path)
            logger.debug("Wrote %s", plev_path.name)

        self._written_dates.append(date_str)

    def finalize(self) -> None:
        """Merge daily files into yearly files using xr.open_mfdataset.

        Replaces CDO mergetime with pure xarray.
        """
        # Collect years from written dates
        years = set()
        for d in self._written_dates:
            years.add(d[:4])

        # Also discover years from existing daily files
        for f in self.daily_dir.glob("dSURF_*.nc"):
            years.add(f.stem.split("_")[1][:4])

        for year in sorted(years):
            for prefix, out_prefix in [("dSURF", "SURF"), ("dPLEV", "PLEV")]:
                pattern = f"{prefix}_{year}*.nc"
                files = sorted(self.daily_dir.glob(pattern))
                if not files:
                    continue

                out_path = self.yearly_dir / f"{out_prefix}_{year}.nc"
                logger.info("Merging %d files → %s", len(files), out_path.name)

                ds = xr.open_mfdataset(
                    [str(f) for f in files],
                    combine="by_coords",
                    parallel=True,
                )
                self._atomic_write(ds.load(), out_path)
                ds.close()

        logger.info("NetCDF yearly merge complete")


class ZarrWriter(BaseWriter):
    """Write ERA5 data to a single Zarr store.

    Output structure:
        {output_dir}/ERA5.zarr
            Surface vars: dims (time, latitude, longitude)
            Pressure vars: dims (time, level, latitude, longitude)
            Surface geopotential stored as 'z_surf' to avoid name clash.

    Uses Blosc/lz4 compression with bitshuffle.

    Args:
        output_dir: Base output directory.
        chunks: Chunk sizes for the zarr store.
    """

    def __init__(
        self,
        output_dir: str | Path,
        chunks: dict | None = None,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.zarr_path = self.output_dir / "ERA5.zarr"
        self.chunks = chunks or {
            "time": 365 * 24,
            "latitude": 3,
            "longitude": 3,
            "level": -1,
        }
        self._initialized = False

    def _get_compressor(self):
        """Get Blosc compressor, handling different zarr versions."""
        try:
            from zarr.codecs import BloscCodec
            return BloscCodec(cname="lz4", clevel=5, shuffle="bitshuffle")
        except ImportError:
            try:
                from zarr.codecs._blosc import BloscCodec
                return BloscCodec(cname="lz4", clevel=5, shuffle="bitshuffle")
            except ImportError:
                from zarr import Blosc
                return Blosc(cname="lz4", clevel=5, shuffle=Blosc.BITSHUFFLE)

    def write_day(
        self, ds_surf: xr.Dataset, ds_plev: xr.Dataset, date_str: str
    ) -> None:
        """Append one day of data to the zarr store.

        Surface geopotential 'z' is renamed to 'z_surf' to avoid clash
        with pressure-level geopotential.
        """
        # Rename surface z to z_surf
        if "z" in ds_surf:
            ds_surf = ds_surf.rename({"z": "z_surf"})

        # Merge surface and plev into one dataset
        ds = xr.merge([ds_surf, ds_plev])

        if not self._initialized or not self.zarr_path.exists():
            # Create the zarr store
            compressor = self._get_compressor()
            encoding = {
                var: {"compressors": compressor} for var in ds.data_vars
            }
            ds.chunk(self.chunks).to_zarr(
                str(self.zarr_path), mode="w", encoding=encoding
            )
            self._initialized = True
            logger.debug("Created zarr store: %s", self.zarr_path)
        else:
            # Append with region='auto'
            ds.to_zarr(str(self.zarr_path), mode="a", region="auto", align_chunks=True)
            logger.debug("Appended %s to zarr store", date_str)

    def finalize(self) -> None:
        """Consolidate zarr metadata."""
        import zarr

        if self.zarr_path.exists():
            zarr.consolidate_metadata(str(self.zarr_path))
            logger.info("Zarr metadata consolidated: %s", self.zarr_path)
