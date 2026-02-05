"""Command-line interface for ecmwf-downloader.

Two subcommands:
    ecmwf-downloader era5      — download ERA5 reanalysis data
    ecmwf-downloader forecast  — download IFS open data forecasts
"""

from __future__ import annotations

import argparse
import logging
import sys


def _setup_logging(verbose: bool = False) -> None:
    """Configure root logger for CLI use."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _add_era5_parser(subparsers) -> None:
    """Add the 'era5' subcommand."""
    p = subparsers.add_parser("era5", help="Download ERA5 reanalysis data")

    p.add_argument(
        "--backend",
        choices=["google", "cds", "s3zarr", "openmeteo"],
        default="google",
        help="Data source backend (default: google)",
    )
    p.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        required=True,
        metavar=("W", "S", "E", "N"),
        help="Bounding box: west south east north",
    )
    p.add_argument("--start", required=True, help="Start date (YYYY-MM-DD)")
    p.add_argument("--end", required=True, help="End date (YYYY-MM-DD)")
    p.add_argument(
        "--levels",
        nargs="+",
        type=int,
        default=[300, 500, 700, 850, 1000],
        help="Pressure levels in hPa (default: 300 500 700 850 1000)",
    )
    p.add_argument(
        "--output-format",
        choices=["netcdf", "zarr"],
        default="netcdf",
        help="Output format (default: netcdf)",
    )
    p.add_argument(
        "--output-dir",
        default="./inputs/climate/",
        help="Output directory (default: ./inputs/climate/)",
    )
    p.add_argument(
        "--time-resolution",
        choices=["1H", "2H", "3H", "6H"],
        default="1H",
        help="Time resolution (default: 1H)",
    )
    p.add_argument(
        "--max-workers",
        type=int,
        default=4,
        help="Number of parallel download threads (default: 4)",
    )
    p.add_argument(
        "--no-rh",
        action="store_true",
        help="Skip computing relative humidity from specific humidity",
    )
    p.add_argument(
        "--no-skip",
        action="store_true",
        help="Re-download even if files already exist",
    )
    p.add_argument(
        "--zarr-url",
        default="",
        help="Zarr store URL (required for s3zarr backend)",
    )
    p.add_argument(
        "--openmeteo-model",
        choices=["era5", "ifs"],
        default="era5",
        help="Open-Meteo model: era5 (surface only, 1940+) or ifs (surface+plev, 2022+)",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")


def _add_forecast_parser(subparsers) -> None:
    """Add the 'forecast' subcommand."""
    p = subparsers.add_parser("forecast", help="Download IFS open data forecasts")

    p.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        required=True,
        metavar=("W", "S", "E", "N"),
        help="Bounding box: west south east north",
    )
    p.add_argument(
        "--output-dir",
        default="./inputs/climate/forecast/",
        help="Output directory (default: ./inputs/climate/forecast/)",
    )
    p.add_argument(
        "--output-timestep",
        choices=["1H", "2H", "3H"],
        default="1H",
        help="Output timestep (default: 1H)",
    )
    p.add_argument(
        "--forecast-time",
        type=int,
        choices=[0, 12],
        default=0,
        help="Forecast base time in UTC (default: 0)",
    )
    p.add_argument(
        "--levels",
        nargs="+",
        type=int,
        default=[1000, 925, 850, 700, 600, 500, 400, 300],
        help="Pressure levels in hPa",
    )
    p.add_argument(
        "--no-backfill",
        action="store_true",
        help="Skip backfilling missing historical forecasts",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")


def _run_era5(args) -> None:
    """Execute the era5 subcommand."""
    from ecmwf_downloader.loader import ERA5Loader

    backend_kwargs = {}
    if args.backend == "s3zarr":
        if not args.zarr_url:
            print("ERROR: --zarr-url is required for s3zarr backend", file=sys.stderr)
            sys.exit(1)
        backend_kwargs["zarr_url"] = args.zarr_url
    elif args.backend == "openmeteo":
        backend_kwargs["start_date"] = args.start
        backend_kwargs["end_date"] = args.end
        backend_kwargs["model"] = args.openmeteo_model

    loader = ERA5Loader(
        backend=args.backend,
        bbox=tuple(args.bbox),
        start_date=args.start,
        end_date=args.end,
        pressure_levels=args.levels,
        output_format=args.output_format,
        output_dir=args.output_dir,
        time_resolution=args.time_resolution,
        max_workers=args.max_workers,
        compute_rh=not args.no_rh,
        skip_existing=not args.no_skip,
        backend_kwargs=backend_kwargs,
    )
    loader.download()


def _run_forecast(args) -> None:
    """Execute the forecast subcommand."""
    from ecmwf_downloader.ifs_forecast import IFSForecastLoader

    fc = IFSForecastLoader(
        bbox=tuple(args.bbox),
        output_dir=args.output_dir,
        output_timestep=args.output_timestep,
        forecast_time=args.forecast_time,
        pressure_levels=args.levels,
    )
    fc.download(backfill=not args.no_backfill)


def main(argv: list[str] | None = None) -> None:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="ecmwf-downloader",
        description="Download ECMWF climate data (ERA5 + IFS forecasts)",
    )
    subparsers = parser.add_subparsers(dest="command")

    _add_era5_parser(subparsers)
    _add_forecast_parser(subparsers)

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    _setup_logging(verbose=args.verbose)

    if args.command == "era5":
        _run_era5(args)
    elif args.command == "forecast":
        _run_forecast(args)


if __name__ == "__main__":
    main()
