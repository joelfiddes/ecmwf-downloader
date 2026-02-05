# Open-Meteo Backend

Free, credential-free backend using [Open-Meteo](https://open-meteo.com/) APIs. Two modes controlled by the `model` parameter:

- **`model="era5"`** (default): ERA5 reanalysis via the archive API. Surface variables only, 1940-present, ~0.25 deg resolution.
- **`model="ifs"`**: IFS historical forecasts via the historical forecast API. Surface + pressure levels, 2022-present, ~9km native resolution.

## Installation

```bash
pip install ecmwf-downloader[openmeteo]
```

## Usage

### Python API

```python
from ecmwf_downloader import ERA5Loader

# ERA5 mode (surface only, 1940+)
loader = ERA5Loader(
    backend="openmeteo",
    bbox=(68.0, 39.0, 69.0, 39.8),
    start_date="2020-01-01",
    end_date="2020-01-03",
    pressure_levels=[500, 700, 850, 1000],
    output_dir="./inputs/climate/",
    time_resolution="3H",
    backend_kwargs={
        "start_date": "2020-01-01",
        "end_date": "2020-01-03",
    },
)
loader.download()

# IFS mode (surface + pressure levels, 2022+)
loader = ERA5Loader(
    backend="openmeteo",
    bbox=(68.0, 39.0, 69.0, 39.8),
    start_date="2023-01-01",
    end_date="2023-01-03",
    pressure_levels=[500, 700, 850, 1000],
    output_dir="./inputs/climate/",
    time_resolution="3H",
    backend_kwargs={
        "start_date": "2023-01-01",
        "end_date": "2023-01-03",
        "model": "ifs",
    },
)
loader.download()
```

### CLI

```bash
# ERA5 mode
ecmwf-downloader era5 --backend openmeteo --bbox 68.0 39.0 69.0 39.8 \
    --start 2020-01-01 --end 2020-01-03 --time-resolution 3H

# IFS mode
ecmwf-downloader era5 --backend openmeteo --openmeteo-model ifs \
    --bbox 68.0 39.0 69.0 39.8 \
    --start 2023-01-01 --end 2023-01-03 --time-resolution 3H
```

## Variable Availability

### Surface Variables

| Variable | ERA5 short name | ERA5 mode | IFS mode | Notes |
|----------|----------------|-----------|----------|-------|
| 2m temperature | `t2m` | yes | yes | C to K |
| 2m dewpoint | `d2m` | yes | yes | C to K |
| Surface pressure | `sp` | yes | yes | hPa to Pa |
| Shortwave radiation | `ssrd` | yes | yes | W/m2 to J/m2 per timestep |
| Total precipitation | `tp` | yes | yes | mm to m per timestep |
| Surface geopotential | `z` | yes | yes | From DEM elevation x g |
| Longwave radiation | `strd` | **no** | **no** | Not available on Open-Meteo |

### Pressure-Level Variables (IFS mode only)

| Variable | ERA5 short name | Notes |
|----------|----------------|-------|
| Temperature | `t` | C to K |
| Geopotential | `z` | gpm to m2/s2 (x g) |
| Relative humidity | `r` | % (no conversion) |
| U-wind | `u` | From speed+direction decomposition |
| V-wind | `v` | From speed+direction decomposition |
| Specific humidity | `q` | **Not available** (r is present) |

### Available Pressure Levels (IFS mode)

1000, 975, 950, 925, 900, 875, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 30 hPa

## Architecture

The backend uses a **pre-fetch strategy**: Open-Meteo returns the full date range in one request per grid point. On the first `fetch_day()` call, the backend:

1. Builds a 0.25 deg grid from the bounding box (matching Google/CDS grid snapping)
2. Pre-fetches the entire date range for all grid points in batches of 20
3. Caches the result in memory as xarray Datasets
4. Subsequent `fetch_day()` calls slice from the cache

This makes it extremely fast for multi-day downloads.

### Accumulated Variables

For `ssrd` and `tp`, the hourly values are summed over the time-resolution window using a rolling sum. For example, at 3H resolution the output represents the 3-hour accumulation ending at each timestamp, matching ERA5/IFS convention.

### Rate Limiting

A sliding-window rate limiter throttles requests to 500/minute to respect Open-Meteo's API limits.

## Known Limitations

- No longwave radiation (`strd`) available from Open-Meteo
- No specific humidity (`q`) on pressure levels (relative humidity `r` is available)
- Surface geopotential `z` is derived from Open-Meteo's DEM elevation (native ~9km), not the ERA5 model terrain (~31km). This creates differences in surface temperature and pressure compared to ERA5 backends in mountainous regions.
- ERA5 mode has no pressure levels (use IFS mode for surface + plev)
- IFS mode only available from 2022 onwards
