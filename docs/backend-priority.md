# Backend Priority and Source Selection

Strategy: try the fastest source first, fall back if a required variable or time range is unavailable. Longwave radiation (`strd`) can be computed from temperature and humidity when not available from the source.

## Priority Table

| Time range | Need plev? | Priority 1 | Priority 2 | Priority 3 |
|------------|-----------|------------|------------|------------|
| 2022+ | yes | **OM IFS** (~1s/day) | Google (~14s/day) | CDS (queued) |
| 2022+ | no | **OM ERA5** (~0.4s/day) | Google (~14s/day) | CDS (queued) |
| 1940-2021 | yes | **Google** (~14s/day) | CDS (queued) | - |
| 1940-2021 | no | **OM ERA5** (~0.4s/day) | Google (~14s/day) | CDS (queued) |

## Variable Source Map

Each variable lists the preferred source in priority order.

### Surface

| Variable | Description | Source priority | Notes |
|----------|-------------|----------------|-------|
| t2m | 2m temperature | OM > Google > CDS | All sources, OM fastest |
| d2m | 2m dewpoint | OM > Google > CDS | All sources |
| sp | Surface pressure | OM > Google > CDS | All sources |
| ssrd | Shortwave radiation | OM > Google > CDS | All sources; OM sums hourly values per timestep |
| tp | Precipitation | OM > Google > CDS | All sources |
| z | Surface geopotential | Google > CDS > OM | OM uses DEM elevation not model terrain; use Google/CDS if model terrain needed |
| strd | Longwave radiation | Google > CDS > **compute** | Not on OM; compute from t2m, d2m, cloud cover |

### Pressure levels

| Variable | Description | Source priority (2022+) | Source priority (pre-2022) |
|----------|-------------|------------------------|---------------------------|
| t | Temperature | OM IFS > Google > CDS | Google > CDS |
| z | Geopotential | OM IFS > Google > CDS | Google > CDS |
| u | U-wind | Google > CDS > OM IFS | Google > CDS |
| v | V-wind | Google > CDS > OM IFS | Google > CDS |
| r | Relative humidity | OM IFS > Google > CDS | Google > CDS |
| q | Specific humidity | Google > CDS | Google > CDS |

Notes on pressure-level wind: OM IFS wind shows large differences vs ECMWF open data (~300% relative for u) due to 9km vs 25km resolution. For applications sensitive to wind accuracy at 0.25 deg, prefer Google/CDS. For applications that benefit from higher resolution wind fields, OM IFS may be preferable.

## Decision Logic

```
given: start_date, end_date, required_variables

1. Determine if pressure levels are needed
   plev_needed = any variable in {t, z, u, v, r, q} is required

2. Determine time range category
   if start_date >= 2022-01-01:
       era = "post2022"
   else:
       era = "pre2022"

3. Check if strd is required
   strd_needed = "strd" in required_variables

4. Check if q is required (vs r being sufficient)
   q_needed = "q" in required_variables and "r" not sufficient

5. Select backend

   if era == "post2022" and plev_needed:
       if q_needed:
           use google  (OM IFS has r but not q)
       else:
           try openmeteo model=ifs
           fallback: google
           fallback: cds

   elif era == "post2022" and not plev_needed:
       try openmeteo model=era5
       fallback: google
       fallback: cds

   elif era == "pre2022" and plev_needed:
       try google
       fallback: cds

   elif era == "pre2022" and not plev_needed:
       try openmeteo model=era5
       fallback: google
       fallback: cds

6. If strd_needed and backend is openmeteo:
       compute strd from available variables

7. If surface geopotential z accuracy matters:
       prefer google/cds over openmeteo
       (OM uses 9km DEM, Google/CDS use ERA5 model terrain)
```

## Computed Variables

### strd (longwave radiation)

When `strd` is not available from the source, it can be estimated from:
- `t2m` (2m temperature)
- `d2m` (2m dewpoint) -> vapour pressure
- Optional: cloud fraction (improves accuracy)

Common parameterisations:
- **Dilley & O'Brien (1998)**: clear-sky from t2m and vapour pressure
- **Prata (1996)**: clear-sky emissivity from precipitable water
- **Konzelmann et al. (1994)**: includes cloud correction

These are already available in TopoPyScale's `meteo_util.py`.

### q (specific humidity)

When `q` is unavailable but `r` (relative humidity) and `t` (temperature) are present on pressure levels, `q` can be computed:
```
e_sat = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))
e = r/100 * e_sat
q = 0.622 * e / (p - 0.378 * e)
```
The loader already computes `r` from `q` when missing; the reverse is equally straightforward.

## Speed Benchmarks

Measured on Fan Mountains bbox (68.0-69.0E, 39.0-39.8N), 3 days, 3H resolution, 4x5 grid:

| Backend | Surface only | Surface + plev | Notes |
|---------|-------------|---------------|-------|
| OM ERA5 | 0.4 s/day | - | No plev available |
| OM IFS | - | 1.0 s/day | 2022+ only |
| Google | 14 s/day | 14 s/day | Includes plev |
| CDS | ~30-60 s/day | ~30-60 s/day | Queue dependent |

Open-Meteo pre-fetches the full date range on first call, so per-day cost is amortised. Longer ranges have even better per-day performance.
