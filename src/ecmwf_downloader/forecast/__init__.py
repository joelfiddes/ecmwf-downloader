"""Forecast data module.

Provides unified access to IFS forecast data from multiple backends:
- ECMWF OpenData: Official ECMWF open data (GRIB format)
- Open-Meteo IFS: Historical IFS forecasts via REST API

Example::

    from ecmwf_downloader.forecast import ForecastLoader

    loader = ForecastLoader(
        bbox=(7.7, 45.95, 7.85, 46.05),
        output_dir="./forecast/",
        backend="auto",  # or "ecmwf_opendata", "openmeteo_ifs"
    )
    loader.download()
    ds = loader.open()
"""

from ecmwf_downloader.forecast.base import ForecastBackend
from ecmwf_downloader.forecast.loader import ForecastLoader
from ecmwf_downloader.forecast.backends import (
    ECMWFOpenDataBackend,
    OpenMeteoIFSBackend,
)

__all__ = [
    "ForecastLoader",
    "ForecastBackend",
    "ECMWFOpenDataBackend",
    "OpenMeteoIFSBackend",
]
