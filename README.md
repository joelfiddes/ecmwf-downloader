# ecmwf-downloader

Standalone package for downloading ECMWF climate data (ERA5 reanalysis + IFS forecasts) from multiple sources. Pure Python — no CDO dependency.

## Install

```bash
pip install ecmwf-downloader[all]
```

Or install only the backends you need:

```bash
pip install ecmwf-downloader[google]   # Google ARCO-ERA5 (anonymous)
pip install ecmwf-downloader[cds]      # CDS API (requires ~/.cdsapirc)
pip install ecmwf-downloader[s3]       # S3 zarr stores
pip install ecmwf-downloader[forecast] # IFS open data forecasts
```

## Quick Start

### Python API

```python
from ecmwf_downloader import ERA5Loader

loader = ERA5Loader(
    backend="google",
    bbox=(6.5, 60.0, 8.5, 61.5),
    start_date="2020-01-01",
    end_date="2020-12-31",
    pressure_levels=[300, 500, 700, 850, 1000],
    output_format="netcdf",
    output_dir="./inputs/climate/",
    time_resolution="1H",
    max_workers=4,
)
loader.download()
```

### IFS Forecasts

```python
from ecmwf_downloader import IFSForecastLoader

fc = IFSForecastLoader(
    bbox=(59, 32, 81, 45),
    output_dir="./inputs/climate/forecast/",
    output_timestep="1H",
)
fc.download(backfill=True)
```

### CLI

```bash
# ERA5 reanalysis
ecmwf-downloader era5 \
    --backend google --bbox 6.5 60.0 8.5 61.5 \
    --start 2020-01-01 --end 2020-12-31 \
    --levels 300 500 700 850 1000 \
    --output-format netcdf --output-dir ./inputs/climate/ \
    --time-resolution 1H --max-workers 8

# IFS forecast
ecmwf-downloader forecast \
    --bbox 59 32 81 45 --output-dir ./inputs/climate/forecast/ \
    --output-timestep 1H
```

## Backends

| Backend | Source | Auth | Install |
|---------|--------|------|---------|
| `google` | Google ARCO-ERA5 netCDF | Anonymous | `pip install .[google]` |
| `cds` | Copernicus CDS API | `~/.cdsapirc` | `pip install .[cds]` |
| `s3zarr` | S3-compatible zarr store | Env vars | `pip install .[s3]` |

## Output Conventions

All output uses ERA5 native conventions:

| Variable | Convention | Units |
|----------|-----------|-------|
| ssrd, strd | Accumulated per timestep | J/m² |
| tp | Accumulated per timestep | m |
| t, t2m, d2m | Instantaneous | K |
| u, v | Instantaneous | m/s |
| q | Instantaneous | kg/kg |
| r | Instantaneous | % (0-100) |
| sp | Instantaneous | Pa |
| z | Geopotential | m² s⁻² |

Standard dimensions: `time`, `latitude`, `longitude`, `level` (hPa ascending).
