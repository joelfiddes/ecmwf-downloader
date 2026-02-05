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

## Summary

- Open-Meteo is **14-34x faster** than Google ARCO-ERA5
- Pressure-level temperature and geopotential are virtually identical (0.2-0.3% difference)
- Surface variables differ due to terrain resolution (9km vs 25km DEM)
- Radiation accumulation convention is correct after fix
- Missing variables: strd (longwave), q (specific humidity on pressure levels)
- Free, no credentials required
