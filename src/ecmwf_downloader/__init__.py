"""ecmwf-downloader: Download ECMWF climate data from multiple sources."""

from ecmwf_downloader._version import __version__
from ecmwf_downloader.bbox import BBox
from ecmwf_downloader.derived import compute_relative_humidity, compute_surface_geopotential
from ecmwf_downloader.loader import ERA5Loader
from ecmwf_downloader.ifs_forecast import IFSForecastLoader

__all__ = [
    "__version__",
    "BBox",
    "ERA5Loader",
    "IFSForecastLoader",
    "compute_relative_humidity",
    "compute_surface_geopotential",
]
