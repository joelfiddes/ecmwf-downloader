"""Forecast data backends."""

from ecmwf_downloader.forecast.backends.ecmwf_opendata import ECMWFOpenDataBackend
from ecmwf_downloader.forecast.backends.openmeteo import OpenMeteoIFSBackend

__all__ = ["ECMWFOpenDataBackend", "OpenMeteoIFSBackend"]
