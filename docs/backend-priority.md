# Backend Priority and Source Selection

> **This document governs backend/source selection logic in nwp-downloader.**
> Code should implement these rules. Update this doc as priorities evolve.

---

## Source Definitions

Formal definitions for each NWP source, organized by module.

---

### Reanalysis Sources (`nwp_downloader.era5`)

#### Global

| Source | ID | Resolution | Temporal Coverage | Endpoint | Auth |
|--------|-----|------------|-------------------|----------|------|
| Google ARCO-ERA5 | `google` | 31 km | 1940-present | `gs://gcp-public-data-arco-era5` | Anonymous |
| CDS (Copernicus) | `cds` | 31 km | 1940-present | `cds.climate.copernicus.eu` | API key |
| OpenMeteo ERA5 | `openmeteo` | 31 km | 1940-present | `archive-api.open-meteo.com` | Anonymous |
| OpenMeteo ERA5-Land | `openmeteo_era5land` | 9 km* | 1950-present | `archive-api.open-meteo.com` | Anonymous |

*ERA5-Land precip is NOT genuinely 9km — it's 31km ERA5 precip fed to a 9km LSM rerun. Use IFS for true 9km precip.

#### Regional

| Source | ID | Resolution | Spatial Coverage | Temporal | Variables | S3 URL | Auth |
|--------|-----|------------|------------------|----------|-----------|--------|------|
| Switch S3 (Central Asia) | `s3zarr` | 31 km | 43-90°E, 24-58°N | 1940-2023 | 15 (10 surf + 5 plev) | `s3://spi-pamir-c7-sdsc/...` | S3 keys |

**Variables in S3 store:** t2m, d2m, sp, z_surf, ssrd, strd, tisr, tp, u10, v10 (surface) + t, z, u, v, q (pressure levels @ 8 levels)

**Regional sources are preferred when bbox fits within coverage** - faster, pre-subset.

---

### Forecast Sources (`nwp_downloader.forecast`)

| Source | ID | Resolution | Temporal Coverage | Endpoint | Auth |
|--------|-----|------------|-------------------|----------|------|
| ECMWF OpenData | `ecmwf_opendata` | 9 km | present-3d → present+10d | `data.ecmwf.int` | Anonymous |
| OpenMeteo IFS | `openmeteo_ifs` | 9 km | present → present+16d | `api.open-meteo.com` | Anonymous |

**Note:** OpenMeteo IFS also provides historical data back to 2022, but for historical reanalysis use the ERA5 module instead.

---

### Climate Projection Sources (`nwp_downloader.cmip6`) — Future

| Source | ID | Resolution | Temporal Coverage | Endpoint | Auth | Status |
|--------|-----|------------|-------------------|----------|------|--------|
| ESGF | `esgf` | ~100 km | 1850-2100 | `esgf-node.llnl.gov` | OpenID | Planned |
| Google CMIP6 | `google_cmip6` | ~100 km | 1850-2100 | `gs://cmip6` | Anonymous | Planned |

**Variables:** tas, pr, hurs, rsds, rlds, sfcWind, ps, ...

**Scenarios:** historical, ssp126, ssp245, ssp370, ssp585

**Models:** Priority list TBD (e.g., EC-Earth3, MPI-ESM1-2-HR, UKESM1-0-LL, GFDL-ESM4)

---

## Detailed Source Specifications

### Reanalysis (`nwp_downloader.era5`)

#### `google` — Google ARCO-ERA5

```yaml
id: google
type: reanalysis
product: ERA5
endpoint: gs://gcp-public-data-arco-era5/raw/
format: NetCDF (via GCS)
resolution_km: 31
bbox: [-180, -90, 180, 90]  # Global
time_range: [1940-01-01, present - 5 days]
variables:
  surface: [t2m, d2m, sp, z, ssrd, strd, tp, u10, v10, msl]
  pressure: [t, z, u, v, q, r]
  levels: [300, 500, 700, 850, 1000]
auth: anonymous
latency: ~14s/day
notes: Full ERA5 mirror, reliable, complete variable set
```

#### `cds` — Copernicus Climate Data Store

```yaml
id: cds
type: reanalysis
product: ERA5
endpoint: https://cds.climate.copernicus.eu/api/v2
format: NetCDF (download)
resolution_km: 31
bbox: [-180, -90, 180, 90]  # Global
time_range: [1940-01-01, present - 5 days]
variables:
  surface: [all ERA5 variables]
  pressure: [all ERA5 variables]
  levels: [all 37 levels]
auth:
  type: api_key
  env_var: CDSAPI_KEY
  config_file: ~/.cdsapirc
latency: 30-60s/day (queue dependent)
notes: Authoritative source, full variable set, but slow due to queues
```

#### `s3zarr` — Switch S3 Central Asia Archive

```yaml
id: s3zarr
type: reanalysis
product: ERA5 (pre-subset)
endpoint: https://os.zhdk.cloud.switch.ch
bucket: spi-pamir-c7-sdsc
zarr_path: era5_data/central_asia.zarr/
full_url: s3://spi-pamir-c7-sdsc/era5_data/central_asia.zarr/
format: Zarr
resolution_km: 31
bbox: [43, 24, 90, 58]  # Central Asia: west, south, east, north
time_range: [1940-01-01, 2023-12-31]
variables:
  # All 15 variables in the store:
  surface: [t2m, d2m, sp, z_surf, ssrd, strd, tisr, tp, u10, v10]
  pressure: [t, z, u, v, q]
  # Note: z is on pressure levels, z_surf is surface geopotential
levels: [300, 500, 600, 700, 800, 850, 900, 1000]  # 8 levels
auth:
  type: s3
  env_vars: [AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY]
  endpoint_url: https://os.zhdk.cloud.switch.ch
latency: ~5s/day (pre-subset, fast)
notes: |
  Regional archive covering Central Asia + surrounding regions.
  Full ERA5 time range (1940-2023). Preferred when bbox fits.
  Requires S3 credentials.
```

#### `openmeteo` — Open-Meteo ERA5 API

```yaml
id: openmeteo
type: reanalysis
product: ERA5
endpoint: https://archive-api.open-meteo.com/v1/archive
format: REST/JSON
resolution_km: 31
bbox: [-180, -90, 180, 90]  # Global
time_range: [1940-01-01, present - 5 days]
variables:
  surface: [t2m, d2m, sp, ssrd, tp, u10, v10]  # No strd!
  pressure: [t, z, u, v, r]  # No q - has r instead
  levels: [300, 500, 700, 850, 1000]
auth: anonymous
latency: ~0.4s/day (surface only)
notes: Very fast, but missing strd and q. Uses DEM elevation not model terrain for z.
```

#### `openmeteo_era5land` — Open-Meteo ERA5-Land API

```yaml
id: openmeteo_era5land
type: reanalysis
product: ERA5-Land
endpoint: https://archive-api.open-meteo.com/v1/archive
format: REST/JSON
resolution_km: 9  # BUT SEE WARNING BELOW
bbox: [-180, -90, 180, 90]  # Global
time_range: [1950-01-01, present - 5 days]
variables:
  surface: [t2m, tp]  # Limited to land surface variables
  pressure: []  # No pressure levels
  levels: []
auth: anonymous
latency: ~0.5s/day
notes: |
  ⚠️ NOT RECOMMENDED FOR PRECIPITATION
  ERA5-Land is a land surface model rerun at 9km, BUT precipitation forcing
  comes from ERA5 at 31km. The "9km" precip is just interpolated 31km data.
  For genuinely higher-resolution precip, use OpenMeteo IFS (2022+).
  ERA5-Land t2m is similarly just interpolated — TPS2 lapse correction is better.
```

---

### Forecast (`nwp_downloader.forecast`)

#### `openmeteo_ifs` — Open-Meteo IFS API

```yaml
id: openmeteo_ifs
type: forecast
product: IFS (derivative)
endpoint: https://api.open-meteo.com/v1/forecast
format: REST/JSON
resolution_km: 9
bbox: [-180, -90, 180, 90]  # Global
time_range: [2022-01-01, present + 16 days]  # Historical + forecast
variables:
  surface: [t2m, d2m, sp, ssrd, tp, u10, v10]
  pressure: [t, z, u, v, r]  # No q
  levels: [300, 500, 700, 850, 1000]
auth: anonymous
latency: ~1s/day
notes: Fast, 2022+ only. Wind differs significantly from Google/CDS at coarse resolution.
```

#### `ecmwf_opendata` — ECMWF Open Data (IFS Forecasts)

```yaml
id: ecmwf_opendata
type: forecast
product: IFS HRES
endpoint: https://data.ecmwf.int/forecasts/
format: GRIB
resolution_km: 9
bbox: [-180, -90, 180, 90]  # Global
time_range: [present - 3 days, present + 10 days]
variables:
  surface: [t2m, d2m, sp, ssrd, strd, tp, u10, v10]
  pressure: [t, z, u, v, q]
  levels: [300, 500, 700, 850, 1000]
auth: anonymous
latency: ~30s/forecast
notes: Official IFS forecasts, GRIB format requires cfgrib
```

---

### Climate Projections (`nwp_downloader.cmip6`) — Planned

#### `esgf` — Earth System Grid Federation

```yaml
id: esgf
type: projection
product: CMIP6
endpoint: https://esgf-node.llnl.gov/esg-search/search
format: NetCDF
resolution_km: ~100 (model dependent)
bbox: [-180, -90, 180, 90]  # Global
time_range: [1850-01-01, 2100-12-31]  # Scenario dependent
variables:
  # CMIP6 variable names differ from ERA5
  surface: [tas, pr, hurs, rsds, rlds, sfcWind, ps]
  pressure: [ta, zg, ua, va, hur]
  levels: [85000, 70000, 50000, 30000]  # Pa, not hPa
auth:
  type: openid
  provider: esgf
scenarios: [historical, ssp126, ssp245, ssp370, ssp585]
models: TBD  # Priority models to support
status: planned
notes: Requires variable name mapping to ERA5 conventions
```

#### `google_cmip6` — Google Cloud CMIP6 Mirror

```yaml
id: google_cmip6
type: projection
product: CMIP6
endpoint: gs://cmip6/
format: Zarr
resolution_km: ~100 (model dependent)
bbox: [-180, -90, 180, 90]  # Global
time_range: [1850-01-01, 2100-12-31]
variables: [same as esgf]
auth: anonymous
status: planned
notes: Cloud-optimized Zarr, faster than ESGF for supported models
```

---

## Priority Logic

### Configuration: `priority_strategy`

```python
ERA5Loader(
    priority_strategy="speed",  # or "reliability"
    ...
)
```

| Strategy | Priority Order | Use Case |
|----------|---------------|----------|
| `speed` (default) | openmeteo > google > cds | Operational, large jobs |
| `reliability` | google > s3zarr > cds > openmeteo | Research, variable completeness |

**Default: `speed`** — 10-30x faster for most jobs, with automatic fallback when variables missing (strd, q, u/v on plev).

```yaml
# Config file example
nwp:
  priority_strategy: speed  # or "reliability"

  # Override for specific jobs needing all variables
  # priority_strategy: reliability
```

### Priority Table (Speed-First)

| Time range | Need plev? | Priority 1 | Priority 2 | Priority 3 |
|------------|-----------|------------|------------|------------|
| 2022+ | yes | **OM IFS** (~1s/day) | Google (~14s/day) | CDS (queued) |
| 2022+ | no | **OM ERA5** (~0.4s/day) | Google (~14s/day) | CDS (queued) |
| 1940-2021 | yes | **Google** (~14s/day) | CDS (queued) | - |
| 1940-2021 | no | **OM ERA5** (~0.4s/day) | Google (~14s/day) | CDS (queued) |

### Regional Override

When bbox fits within a regional source's coverage, prefer the regional source:

```
if bbox_within(request.bbox, s3zarr.bbox) and time_within(request.time, s3zarr.time_range):
    priority.insert(0, "s3zarr")  # Pre-subset = fastest
```

---

## Variable Source Map

Each variable lists the preferred source in priority order.

### Surface Variables

| Variable | Description | Source priority | Notes |
|----------|-------------|----------------|-------|
| t2m | 2m temperature | OM > Google > CDS | All sources |
| d2m | 2m dewpoint | OM > Google > CDS | All sources |
| sp | Surface pressure | OM > Google > CDS | All sources |
| ssrd | Shortwave radiation | OM > Google > CDS | OM sums hourly values per timestep |
| tp | Precipitation | **OM IFS (9km, 2022+)** > Google > CDS | Always use 9km IFS when available |
| z | Surface geopotential | Google > CDS > OM | OM uses DEM not model terrain |
| strd | Longwave radiation | Google > CDS > **compute** | Not on OM; compute from t2m, d2m |
| u10 | 10m U wind | OM > Google > CDS | All sources |
| v10 | 10m V wind | OM > Google > CDS | All sources |

### Pressure Level Variables

| Variable | Description | Source priority (2022+) | Source priority (pre-2022) | Notes |
|----------|-------------|------------------------|---------------------------|-------|
| t | Temperature | OM IFS > Google > CDS | Google > CDS | OM excellent (0.27% diff) |
| z | Geopotential | OM IFS > Google > CDS | Google > CDS | OM excellent (0.22% diff) |
| u | U-wind | **Google > CDS only** | Google > CDS | ⚠️ OM blacklisted (43 m/s error) |
| v | V-wind | **Google > CDS only** | Google > CDS | ⚠️ OM blacklisted (28 m/s error) |
| r | Relative humidity | Google > CDS > OM IFS | Google > CDS | OM moderate (20% diff) |
| q | Specific humidity | **Google > CDS only** | Google > CDS | ⚠️ Not available on OM |

---

## Variable Blacklist (OpenMeteo)

Based on comparison testing (see `openmeteo-comparison.md`), these variables should **NEVER** be sourced from OpenMeteo:

| Variable | Source | Issue | MaxAbsDiff |
|----------|--------|-------|------------|
| u (plev) | OM IFS | Resolution mismatch (9km vs 25km) | **43 m/s** |
| v (plev) | OM IFS | Resolution mismatch (9km vs 25km) | **28 m/s** |
| strd | OM * | Not available | - |
| q | OM * | Not available (only r) | - |
| z (surface) | OM * | Uses different DEM, not model terrain | 12370 m²/s² |

**Always use Google/CDS for:** u, v, strd, q, surface z

**Safe to use from OpenMeteo:**
- t, z on pressure levels (0.2-0.3% difference — excellent)
- ssrd (~4% difference — good)
- tp (negligible difference)
- t2m, d2m, sp (if terrain differences acceptable)

---

## Variable-Specific Source Preferences (Resolution)

### Precipitation (`tp`) — 9km Available

| Source | Resolution | Temporal | Recommendation |
|--------|------------|----------|----------------|
| ERA5 (Google/CDS) | 31 km | 1940-present | Default for pre-2022 |
| OpenMeteo IFS | **9 km** | 2022-present | **Always preferred when available** |
| ERA5-Land (OpenMeteo) | 9 km | 1950-present | ⚠️ Not recommended — see note |

**⚠️ ERA5-Land precipitation caveat:**
ERA5-Land is NOT higher-resolution precipitation forcing. It's ERA5 (31km) precipitation fed through a higher-resolution land surface model rerun. The precipitation input is still 31km. For genuinely higher-resolution precip, use IFS-based products.

```python
# For 2022+ data: always use 9km IFS precip
source_map = {
    "default": "google",
    "tp": "openmeteo_ifs",  # Genuine 9km precipitation (2022+)
}

# For pre-2022: stuck with 31km ERA5
source_map = {
    "default": "google",
    # No 9km precip available before 2022
}
```

**9km IFS precip is always preferred when available (2022+):**
- Genuine 9km resolution (not downscaled)
- Better orographic precipitation representation
- Critical for snow/hydro in complex terrain

### Temperature (`t2m`)

ERA5-Land has 9km t2m, but this is just ERA5 t2m interpolated to 9km grid — minimal benefit.
TPS2 lapse rate correction from pressure levels is more physically meaningful.

### Wind (`u10`, `v10`)

Stick with ERA5 (31km). Wind is already heavily parameterized;
higher resolution doesn't necessarily mean more accurate.

### Radiation (`ssrd`, `strd`)

Stick with ERA5 (31km). Cloud effects dominate; resolution less important.

---

## Decision Logic

```
given: start_date, end_date, required_variables, bbox, priority_strategy="speed"

0. Set base priority order from strategy
   if priority_strategy == "speed":
       base_order = [openmeteo, google, cds]
   else:  # reliability
       base_order = [google, s3zarr, cds, openmeteo]

1. Check regional source applicability
   if bbox_within(bbox, s3zarr.bbox) and time_within([start, end], s3zarr.time_range):
       regional_available = True
       # S3 Zarr becomes top priority

2. Determine if pressure levels are needed
   plev_needed = any variable in {t, z, u, v, r, q} is required

3. Determine time range category
   if start_date >= 2022-01-01:
       era = "post2022"
   else:
       era = "pre2022"

4. Check special variable requirements
   strd_needed = "strd" in required_variables
   q_needed = "q" in required_variables and "r" not sufficient
   high_res_precip = config.get("high_res_precip", False)

5. Select backend

   if regional_available:
       try s3zarr
       # Fall through if fails

   # Check for blacklisted variables that force Google/CDS
   blacklisted_vars = {"u", "v", "q", "strd"} & set(required_variables)
   if blacklisted_vars:
       # Must use Google/CDS for these — OpenMeteo has errors or missing
       use google
       fallback: cds
       # Note: Can still fetch tp separately from OM IFS for 9km precip

   elif era == "post2022" and plev_needed:
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

6. Handle computed/missing variables
   if strd_needed and backend is openmeteo:
       compute strd from t2m, d2m, cloud cover

   if q_needed and backend has only r:
       compute q from r, t, p

7. Handle high-res precip (always preferred when available)
   if start_date >= 2022-01-01 and "tp" in required_variables:
       fetch tp from openmeteo_ifs (genuine 9km)
       return as separate dataset for TPS2 to align
   # Note: Do NOT use ERA5-Land for precip — it's fake 9km (31km input)

8. Surface geopotential accuracy
   if z accuracy matters (model terrain vs DEM):
       prefer google/cds over openmeteo
```

---

## Source Mixing Rules

When using multiple sources (e.g., 9km precip + 31km temp):

1. **Fetch independently** - don't try to align grids in downloader
2. **Return as dict** - `{"google": ds_main, "openmeteo_era5land": ds_precip}`
3. **TPS2 handles alignment** - interpolates each to unit centroids

**Constraint:** All sources must cover the same time range.
If one source has gaps, either:
- Fill from fallback source
- Raise error and let user decide

---

## Computed Variables

### strd (longwave radiation)

When `strd` is not available from the source (e.g., OpenMeteo), estimate from:
- `t2m` (2m temperature)
- `d2m` (2m dewpoint) → vapour pressure
- Optional: cloud fraction (improves accuracy)

Common parameterisations:
- **Dilley & O'Brien (1998)**: clear-sky from t2m and vapour pressure
- **Prata (1996)**: clear-sky emissivity from precipitable water
- **Konzelmann et al. (1994)**: includes cloud correction

These are already available in TopoPyScale's `meteo_util.py`.

### q (specific humidity)

When `q` is unavailable but `r` (relative humidity) and `t` (temperature) are present:
```
e_sat = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))
e = r/100 * e_sat
q = 0.622 * e / (p - 0.378 * e)
```

---

## Speed Benchmarks

Measured on Fan Mountains bbox (68.0-69.0E, 39.0-39.8N), 3 days, 3H resolution, 4x5 grid:

| Backend | Surface only | Surface + plev | Notes |
|---------|-------------|---------------|-------|
| S3 Zarr | ~5 s/day | ~5 s/day | Regional, pre-subset |
| OM ERA5 | 0.4 s/day | - | No plev available |
| OM IFS | - | 1.0 s/day | 2022+ only |
| Google | 14 s/day | 14 s/day | Includes plev |
| CDS | 30-60 s/day | 30-60 s/day | Queue dependent |

Open-Meteo pre-fetches the full date range on first call, so per-day cost is amortised.

---

## Fallback Logic

If preferred source fails:

```
1. Try preferred source
2. If unavailable/error → try next in priority
3. If all fail → raise clear error with suggestions
```

Example:
```
Attempting openmeteo... [FAILED: strd not available]
Falling back to google... [SUCCESS]
```

---

## Configuration Examples

### Default (speed priority)
```python
ERA5Loader(
    priority_strategy="speed",  # Default - fastest source first
    ...
)
```

### Reliability priority
```python
ERA5Loader(
    priority_strategy="reliability",  # Google first, complete variables
    ...
)
```

### Explicit backend (bypass priority logic)
```python
ERA5Loader(backend="google", ...)  # Force specific backend
```

### High-res precipitation
```python
ERA5Loader(
    source_map={
        "default": "google",
        "tp": "openmeteo_era5land",
    },
    ...
)
```

### Regional archive
```python
ERA5Loader(
    backend="s3zarr",
    backend_kwargs={"zarr_url": "s3://spi-pamir-c7-sdsc/era5_data/central_asia.zarr/"},
    ...
)
```

### Explicit fallback chain
```python
ERA5Loader(
    backend_priority=["s3zarr", "google", "cds"],
    ...
)
```

### Offline mode (local cache only)
```python
ERA5Loader(
    backend="local",  # Only read from cache, never download
    cache_dir="./cache/",
    ...
)
```

---

## Future Additions

As new sources become available, add entries here:

- [ ] ICON (DWD) - European high-res
- [ ] GFS (NOAA) - Global, free
- [ ] COSMO - Alpine region
- [ ] CMIP6 - Climate projections
- [ ] Custom user sources - local NWP output

---

## Changelog

| Date | Change |
|------|--------|
| 2026-02-06 | Merged source definitions, added regional sources, 9km precip |
| Previous | Initial version with speed-first priority |
