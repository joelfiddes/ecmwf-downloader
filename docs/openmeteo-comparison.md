# Open-Meteo Backend Comparison Results

Comparison tests performed on the Fan Mountains, Tajikistan bbox (68.0, 39.0, 69.0, 39.8) with pressure levels [500, 700, 850, 1000] at 3H resolution.

## 1. Speed: Open-Meteo ERA5 vs Google ARCO-ERA5

Date range: 2020-01-01 to 2020-01-03, surface only.

| Backend | Total (s) | Per day (s) | Speedup |
|---------|----------|------------|---------|
| **openmeteo** | **1.2** | **0.4** | **34x** |
| google | 41.3 | 13.8 | 1x |

Open-Meteo pre-fetches all days in one batch of API calls, then slices from cache.

## 2. Speed: Open-Meteo IFS vs Google ARCO-ERA5

Date range: 2023-01-01 to 2023-01-03, surface + pressure levels.

| Backend | Total (s) | Per day (s) | Speedup |
|---------|----------|------------|---------|
| **openmeteo** | **3.0** | **1.0** | **14x** |
| google | 40.9 | 13.6 | 1x |

## 3. Open-Meteo ERA5 vs Google ARCO-ERA5 (Surface)

Both sourcing ERA5 reanalysis, same 4x5 grid (0.25 deg), 2020-01-01 to 2020-01-03.

| Variable | MaxAbsDiff | MeanAbsDiff | Notes |
|----------|-----------|-------------|-------|
| d2m | 8.2 K | 2.5 K | Terrain-driven (different DEM) |
| sp | 12258 Pa | 3552 Pa | Terrain-driven |
| ssrd | 1794 J/m2 | 471 J/m2 | Excellent (~0.1% relative) |
| t2m | 8.2 K | 2.5 K | Terrain-driven |
| tp | 0.0001 m | 0.00002 m | Negligible |
| z | 12370 m2/s2 | 3815 m2/s2 | Different DEM source |
| strd | - | - | Only available from Google |

The surface variable differences are driven by Open-Meteo using its own DEM (~9km) for surface geopotential rather than the ERA5 model terrain (~31km). This cascades to temperature and pressure in mountainous terrain. Radiation (ssrd) and precipitation (tp) agree very closely.

## 4. Open-Meteo IFS vs ECMWF Open Data (IFS Forecast)

Direct comparison of IFS forecast data from Open-Meteo historical forecast API vs ecmwf-opendata package, for 2026-02-03 at the same 0.25 deg grid points.

**Important context**: Open-Meteo serves IFS at its native ~9km resolution, while ecmwf-opendata serves the 0.25 deg (~25km) open data product. Even when querying at the same lat/lon, Open-Meteo returns the nearest 9km grid point value.

### Surface - Instantaneous Variables

| Variable | MaxAbsDiff | MeanAbsDiff | RelDiff% |
|----------|-----------|-------------|----------|
| t2m | 17.8 K | 4.9 K | 1.9% |
| d2m | 9.5 K | 3.0 K | 1.2% |
| sp | 7850 Pa | 2880 Pa | 3.9% |

### Surface - Accumulated Variables (per 3h timestep)

| Variable | ECMWF mean | OM mean | Ratio |
|----------|-----------|---------|-------|
| ssrd | 1,625,469 J/m2 | 1,685,745 J/m2 | 0.96 |
| tp | 0.0003 m | 0.0002 m | 1.48 |

The ssrd accumulation convention is correct (ratio ~0.96, within 4%). The tp ratio is higher but absolute values are tiny (0.3 mm/3h vs 0.2 mm/3h).

### Pressure Levels

| Variable | MaxAbsDiff | MeanAbsDiff | RelDiff% | Assessment |
|----------|-----------|-------------|----------|------------|
| t | 3.0 K | 0.7 K | **0.27%** | Excellent |
| z | 289 m2/s2 | 56 m2/s2 | **0.22%** | Excellent |
| r | 57% | 13% | 20.4% | Moderate (resolution effect) |
| u | 43 m/s | 11 m/s | large | Resolution-driven |
| v | 28 m/s | 10 m/s | large | Resolution-driven |

**Temperature and geopotential on pressure levels agree within 0.3%** - essentially the same data. The wind differences are driven by the 9km vs 25km resolution difference (wind is very sensitive to orography). Humidity shows moderate divergence.

### Variables only in one source

- **ecmwf-opendata only**: strd (longwave radiation), q (specific humidity)
- **Open-Meteo only**: z (surface geopotential from DEM)

## Summary (API)

- Open-Meteo is **14-34x faster** than Google ARCO-ERA5
- Pressure-level temperature and geopotential are virtually identical (0.2-0.3% difference)
- Surface variables differ due to terrain resolution (9km vs 25km DEM)
- Radiation accumulation convention is correct after fix
- Missing variables: strd (longwave), q (specific humidity on pressure levels)
- Free, no credentials required

---

## 5. Open-Meteo S3 Direct Access

An alternative to the API is direct access to `.om` files on `s3://openmeteo`. This section documents the tradeoffs.

### S3 vs API vs Google Overview

| Aspect | Open-Meteo API | Open-Meteo S3 | Google Cloud |
|--------|----------------|---------------|--------------|
| **Grid alignment** | Off-grid (~0.03-0.08°) | Exact 0.25° | Exact 0.25° |
| **Rate limits** | 500 req/min | None | None |
| **Surface variables** | All | Most (no strd) | All |
| **Pressure levels** | IFS mode only (2022+) | None | All |
| **Historical data** | 1940+ (ERA5) | 1950+ (year files) | 1940+ |
| **Credentials** | None | None | None |

### Grid Alignment Issue

**Critical finding**: The API returns coordinates that do not match the standard ERA5 0.25° grid.

Test request for `(60.25, 7.25)`:

| Source | Returned Coordinates | Offset |
|--------|---------------------|--------|
| Open-Meteo S3 | (60.250000, 7.250000) | 0, 0 |
| Google Cloud | (60.250000, 7.250000) | 0, 0 |
| Open-Meteo API | (60.281193, 7.166276) | +0.031°, -0.084° |

This ~3-8 km offset can cause issues when mixing backends or expecting standard grid alignment.

### S3 Performance Benchmarks

#### Small Region (Norway: 7×9 grid points, 4 variables)

**S3 vs Google Cloud** (store opened once for Google):

| Time Range | S3 Direct | Google Cloud | S3 Speedup |
|------------|-----------|--------------|------------|
| 1 day | 17.6s | 42s (+107s startup) | 2.4x |
| 7 days | 16.1s | 71s | 4.4x |
| 21 days | 15.3s | 163s | 10.7x |

**S3 parallel fetching** (4 workers vs sequential):

| Mode | Time | Speedup |
|------|------|---------|
| Sequential | 17.5s | 1x |
| Parallel (4 workers) | 7.2s | 2.4x |

**API vs S3** (small region, under rate limit):

| Time Range | API | S3 (parallel) | Winner |
|------------|-----|---------------|--------|
| 6 hours | 0.13s | 5.4s | API |
| 1 day | 0.13s | 4.1s | API |
| 21 days | 0.60s | 4.2s | API |

For small regions, the API is faster due to S3's ~4s per-file connection overhead.

#### Large Region (Scandinavia: 6,213 grid points)

| Backend | Result |
|---------|--------|
| S3 Direct | 923s (4 vars, 24h) — completed |
| API | **429 Too Many Requests** — rate limited |

For large regions, S3 is the only viable option.

#### Rate Limit Calculation (Central Asia)

- Bbox: 50°E-90°E, 35°N-55°N
- Grid points: 81 × 161 = **13,041 points**
- API requests needed: ~653 (at 20 points/batch)
- Time at rate limit: ~78 seconds minimum
- **Verdict**: Tight; S3 recommended for reliability

### S3 Bucket Structure

**Bucket**: `s3://openmeteo` (us-west-2, public, anonymous access)

```
data/
├── copernicus_era5/              # ERA5 reanalysis
│   ├── temperature_2m/
│   │   ├── chunk_904.om          # Rolling ~4 years (2022+)
│   │   └── ...chunk_975.om
│   └── static/meta.json
│
├── copernicus_era5_land/         # ERA5-Land
│   ├── temperature_2m/
│   │   ├── year_1950.om          # Historical (9 GB/year)
│   │   ├── ...year_2024.om
│   │   ├── chunk_904.om          # Recent data
│   │   └── ...
```

| File Type | Pattern | Coverage | Size |
|-----------|---------|----------|------|
| Chunk | `chunk_XXX.om` | Rolling ~4 years | ~230 MB |
| Year | `year_YYYY.om` | 1950-present | ~9 GB |

Chunk timing: 504 hours (21 days) per chunk. Reference: chunk 975 starts 2026-01-11.

### S3 Available Variables

**copernicus_era5** (chunks only, 2022+):
- temperature_2m, dew_point_2m, pressure_msl, shortwave_radiation, precipitation
- cloud_cover, wind components, soil variables

**copernicus_era5_land** (year files 1950+):
- temperature_2m, dew_point_2m (limited variables in historical files)

**Not available in S3**:
- Surface geopotential (z) — static file has nodata issues
- Longwave radiation (strd)
- Specific humidity (q)
- **Any pressure level data**

### Reading .om Files

```python
import fsspec
from omfiles import OmFileReader

backend = fsspec.open(
    "blockcache::s3://openmeteo/data/copernicus_era5/temperature_2m/chunk_975.om",
    mode="rb",
    s3={"anon": True},
    blockcache={"cache_storage": "/tmp/cache"},
)

with backend as f:
    reader = OmFileReader(f)
    # Shape is (lat, lon, time) - lat descending 90° to -90°
    data = reader[114:121, 746:755, 0:24]  # bbox subset, 24 hours
```

### Recommendations

| Use Case | Recommended Backend |
|----------|---------------------|
| Small region, quick results | API |
| Large region (>1000 points) | S3 |
| Grid consistency required | S3 (not API) |
| Pressure levels needed | Google/CDS |
| Historical data (pre-2022) | S3 year files or Google |
| Surface + pressure combined | Google (or hybrid) |

### Prototype S3 Backend

A prototype implementation exists at `src/OM/openmeteo_s3_backend.py` with:
- Parallel variable fetching (configurable workers)
- Year file support for historical data
- ERA5Backend-compatible interface

To integrate, the backend needs:
- Missing variable handling (z, strd, q)
- Cross-chunk reads for date ranges spanning boundaries
- Registration in backends registry
