# Backend Priority and Source Selection

> **This document governs backend/source selection logic in nwp-downloader.**
> Code should implement these rules. Update this doc as priorities evolve.

---

## TL;DR — Recommended Strategy

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ DEFAULT: Google ARCO-ERA5                                                   │
│   - Fastest for both small AND large regions                                │
│   - Complete variables (strd, q, pressure levels)                           │
│   - No rate limits, no credentials required                                 │
│   - ~11-15s/day regardless of region size                                   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ REGIONAL OVERRIDE: s3zarr                                                   │
│   IF bbox ⊆ Central Asia (43-90°E, 24-58°N) AND time ≤ 2023:               │
│     → s3zarr FOR EVERYTHING (~5s/day, complete vars)                       │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ USE OpenMeteo S3 ONLY WHEN:                                                 │
│   - Exact 0.25° grid alignment required (API returns off-grid coords)      │
│   - Google unavailable                                                      │
│   - ⚠️ S3 is 3-4x SLOWER than Google even for surface-only fetching       │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ USE OpenMeteo API WHEN:                                                     │
│   - Surface only, small region, speed critical (~0.1s/day)                 │
│   - Want 9km IFS precipitation (2022+)                                      │
│   - Caveat: Missing strd, q; coordinates are off-grid                       │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Benchmark Results (2026-02-06)

| Scenario | Google | OpenMeteo S3 | Winner |
|----------|--------|--------------|--------|
| Small region (25 pts, 3 days) | 44.7s | 13.2s | S3 faster |
| Large region (861 pts, 2 days) | 26.5s | 149.1s | **Google 6x faster** |
| Surface-only, medium region (14 days) | 6.2s/day | 23.9s/day | **Google 4x faster** |

**Key insight:** Google ARCO is faster across the board. S3 has ~4-5s per-variable
connection overhead that dominates even for surface-only fetching. **Use Google as default.**

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

### Recommended: Google-First Strategy

Based on benchmarks (2026-02-06), **Google ARCO is the best default**:

```python
ERA5Loader(
    backend="google",  # Recommended default
    # OR
    backend="auto",    # Auto-selects based on bbox/time
    ...
)
```

### Backend Selection Logic

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ STEP 0: Regional Check                                                      │
│   IF bbox ⊆ Central Asia AND time ≤ 2023:                                  │
│     → s3zarr (fastest + complete for this region)                          │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ STEP 1: Default — Google ARCO                                               │
│   - Fastest for both small AND large regions                                │
│   - Complete variables (strd, q, all pressure levels)                       │
│   - Scales well: ~11-15s/day regardless of region size                      │
│   - No credentials required                                                 │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ STEP 2: Fallbacks                                                           │
│   IF Google unavailable:                                                    │
│     → OpenMeteo API (surface only, fast but incomplete)                    │
│     → OpenMeteo S3 (exact grid, but slow for large regions)                │
│     → CDS (complete but slow, requires credentials)                        │
└─────────────────────────────────────────────────────────────────────────────┘
```

### When to Use Hybrid Mode

**⚠️ Generally NOT recommended.** Google ARCO is faster in almost all scenarios.

Hybrid mode (S3 surface + Google plev) may be beneficial when:

1. **Exact 0.25° grid alignment required** — OpenMeteo API returns off-grid coords
2. **9km IFS precipitation needed** — genuine high-res precip (2022+)
3. **Google unavailable** and API would hit rate limits

**⚠️ Benchmark results show Google is faster:**
- Small region: Hybrid is marginally faster for first request only
- Large region: Google **6x faster** than Hybrid
- Surface-only: Google **3-4x faster** than S3

### Why NOT Hybrid by Default?

OpenMeteo S3 has per-file connection overhead (~4-5s per variable file). Even for
surface-only fetching, this overhead makes S3 ~4x slower than Google's Zarr access.
The overhead does not amortize over more days — S3 stays at ~23-25s/day while
Google drops from 8s to 6s/day with longer time periods.

### Single Mode

For simple use cases or debugging, single mode uses one backend for everything:

| Time range | Need plev? | Priority 1 | Priority 2 | Priority 3 |
|------------|-----------|------------|------------|------------|
| Any | yes | **s3zarr** (if regional) | google | cds |
| Any | no | **s3zarr** (if regional) | openmeteo_s3 | google |

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

### Hybrid Mode (Default)

```python
def select_hybrid_sources(bbox, start_date, end_date, pressure_levels):
    """
    Returns dict mapping variable groups to backends.

    Returns:
        {
            "surface": (backend_name, backend_kwargs),
            "precip": (backend_name, backend_kwargs),
            "plev_strd": (backend_name, backend_kwargs),
        }
    """

    # ─── STEP 0: Regional override ───────────────────────────────────
    # s3zarr is fastest AND complete for Central Asia ≤2023
    if bbox_within(bbox, S3ZARR_BBOX) and end_date <= "2023-12-31":
        return {
            "surface": ("s3zarr", {}),
            "precip": ("s3zarr", {}),
            "plev_strd": ("s3zarr", {}),
        }

    # ─── STEP 1: Surface variables ───────────────────────────────────
    # OpenMeteo S3: exact 0.25° grid, no rate limits, fast
    surface_backend = select_with_fallback(
        ["openmeteo_s3", "openmeteo", "google", "cds"],
        required_vars=["t2m", "d2m", "sp", "ssrd", "u10", "v10"],
    )

    # ─── STEP 2: Precipitation ───────────────────────────────────────
    # 2022+: IFS 9km is genuine high-res (not interpolated like ERA5-Land)
    if start_date >= pd.Timestamp("2022-01-01"):
        precip_backend = select_with_fallback(
            ["openmeteo_ifs", "openmeteo_s3", "google", "cds"],
            required_vars=["tp"],
        )
    else:
        precip_backend = select_with_fallback(
            ["openmeteo_s3", "openmeteo", "google", "cds"],
            required_vars=["tp"],
        )

    # ─── STEP 3: Pressure levels + longwave ──────────────────────────
    # Only Google/CDS have: u, v, q on plev + strd
    # OpenMeteo blacklisted: u/v have 43m/s errors, q/strd missing
    if pressure_levels:
        plev_backend = select_with_fallback(
            ["google", "cds"],  # OpenMeteo blacklisted for plev
            required_vars=["t", "z", "u", "v", "q", "strd"],
        )
    else:
        # Surface-only mode: still need strd from Google
        plev_backend = select_with_fallback(
            ["google", "cds"],
            required_vars=["strd"],
        )

    return {
        "surface": surface_backend,
        "precip": precip_backend,
        "plev_strd": plev_backend,
    }
```

### Single Mode

For debugging or simple use cases:

```python
def select_single_source(bbox, start_date, end_date, pressure_levels):
    """Select one backend for all variables."""

    # Regional always wins
    if bbox_within(bbox, S3ZARR_BBOX) and end_date <= "2023-12-31":
        return "s3zarr"

    # Need pressure levels or strd → must use Google/CDS
    if pressure_levels or needs_strd:
        return select_with_fallback(["google", "cds"])

    # Surface only → OpenMeteo S3 preferred
    return select_with_fallback(["openmeteo_s3", "openmeteo", "google", "cds"])
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

### By Region Size (Switzerland, 3 days, 850 hPa)

| Backend | Small (25 pts) | Large (861 pts) | Notes |
|---------|---------------|-----------------|-------|
| **Google ARCO** | 44.7s (14.9s/day) | **26.5s (13.3s/day)** | **Best for large regions** |
| OpenMeteo S3 | 13.2s (4.4s/day) | 149.1s (74.6s/day) | Per-file overhead hurts |
| OpenMeteo API | 0.2s (0.1s/day) | Rate limited | Surface only, fast |
| Hybrid (S3+Google) | 40.8s (13.6s/day) | 162.3s (81.1s/day) | Two fetches overhead |

### Time Scaling — Fair Surface-Only Comparison

Western Alps region (~5° x 3°), surface variables only, same 6 vars both backends:

| Days | Google Total | Google/day | S3 Total | S3/day | Winner |
|------|--------------|------------|----------|--------|--------|
| 1 | 8.2s | 8.2s | 24.8s | 24.8s | Google 3.0x faster |
| 7 | 46.6s | 6.7s | 163.5s | 23.4s | Google 3.5x faster |
| 14 | 86.6s | 6.2s | 334.4s | 23.9s | Google 3.9x faster |

**Key insights:**
- **Google is 3-4x faster than S3** even for surface-only fetching
- **Google benefits from amortization**: startup overhead spreads over more days (8.2s → 6.2s/day)
- **S3 has constant per-day cost**: ~23-25s/day regardless of period length
- S3's per-file connection overhead (~4-5s per variable) dominates
- **Conclusion: Google ARCO is the best default for all scenarios**

### Regional Archive (Central Asia)

| Backend | Per day | Notes |
|---------|---------|-------|
| s3zarr | ~5s | Pre-subset, complete variables |

Open-Meteo API pre-fetches the full date range on first call, so per-day cost is amortised (~0.1s/day after first).

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

### Default (hybrid mode)
```python
ERA5Loader(
    mode="hybrid",  # Default — multi-backend, variable-optimized
    ...
)
# Surface from OpenMeteo S3, precip from IFS (2022+), plev+strd from Google
```

### Single backend (bypass hybrid)
```python
ERA5Loader(
    mode="single",
    backend="google",  # Force single backend for all variables
    ...
)
```

### Regional archive (auto-detected)
```python
# If bbox is within Central Asia and time ≤ 2023, s3zarr is auto-selected
ERA5Loader(
    bbox=(68.0, 39.0, 72.0, 42.0),  # Tajikistan
    start_date="2020-01-01",
    end_date="2020-12-31",
    ...
)
# → Automatically uses s3zarr (fast + complete)
```

### Force regional archive
```python
ERA5Loader(
    mode="single",
    backend="s3zarr",
    backend_kwargs={"zarr_url": "s3://spi-pamir-c7-sdsc/era5_data/central_asia.zarr/"},
    ...
)
```

### Surface-only (no pressure levels)
```python
ERA5Loader(
    pressure_levels=None,  # No pressure levels
    ...
)
# → Surface from OpenMeteo S3, strd from Google (still needed for energy balance)
```

### Explicit fallback chain
```python
ERA5Loader(
    mode="single",
    backend_priority=["s3zarr", "openmeteo_s3", "google", "cds"],
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
| 2026-02-06 | **Hybrid mode**: Multi-backend fetching (S3 surface, IFS precip, Google plev+strd) |
| 2026-02-06 | Added OpenMeteo S3 backend for direct .om file access |
| 2026-02-06 | Merged source definitions, added regional sources, 9km precip |
| Previous | Initial version with speed-first priority |
