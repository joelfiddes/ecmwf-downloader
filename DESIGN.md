# nwp-downloader Design Document

> **Note:** Package to be renamed from `ecmwf-downloader` to `nwp-downloader`

## Purpose

**nwp-downloader** is a unified interface for downloading NWP (Numerical Weather Prediction) data from multiple sources with consistent output conventions.

### Core Problem It Solves

Different NWP data sources have different:
- APIs (CDS REST, Google Cloud Storage, S3, OpenMeteo REST, ESGF)
- Authentication (CDS API key, anonymous, AWS credentials)
- Variable naming (`2m_temperature` vs `t2m` vs `2t`)
- Dimension naming (`valid_time` vs `time`, `pressure_level` vs `level`)
- Coordinate conventions (-180:180 vs 0:360 longitude)
- Data formats (GRIB, NetCDF, Zarr)
- **Resolutions** (ERA5 31km, Open-Meteo 9km precip, CMIP6 ~100km)

nwp-downloader abstracts these differences and provides:
- Unified API per domain: `nwp_downloader.era5`, `nwp_downloader.forecast`, `nwp_downloader.cmip6`
- Consistent output: standardized variable names, dimensions, units
- Bbox/time subsetting at source (efficient, no full-globe downloads)
- **Variable-dependent source selection** (e.g., 9km precip from Open-Meteo, rest from ERA5)

---

## Package Structure

```
nwp_downloader/
├── era5/           # ERA5 reanalysis (1940-present, 5-day delay)
│   ├── loader.py   # ERA5Loader
│   └── backends/   # google, cds, s3zarr, openmeteo
├── forecast/       # IFS/HRES operational forecasts (10-day ahead)
│   ├── loader.py   # ForecastLoader
│   └── backends/   # ecmwf_opendata, openmeteo_ifs
├── cmip6/          # Climate projections (future)
│   ├── loader.py   # CMIP6Loader
│   └── backends/   # esgf, google_cmip6
├── preprocess/     # Transform raw downloads → simulation-ready
│   ├── __init__.py
│   ├── validate.py # Range checks, missing data detection
│   ├── derived.py  # Compute strd, q, etc.
│   ├── conventions.py  # Lat descending, lon -180:180, units
│   └── chunking.py # Rechunk for access pattern
├── common/
│   ├── bbox.py
│   ├── variables.py
│   └── writer.py
└── __init__.py
```

Each module (era5, forecast, cmip6) is a self-contained domain with its own loader and backends, sharing common utilities.

The `preprocess` module transforms raw downloads into simulation-ready datasets.

---

## Variable-Dependent Source Selection

**Key insight:** Different sources excel at different variables. We should be able to mix sources per-variable when efficiency allows.

### Example: Precipitation Resolution

| Source | Precip Resolution | Other Vars |
|--------|-------------------|------------|
| ERA5 (Google/CDS) | 31 km | 31 km |
| Open-Meteo ERA5-Land | **9 km** | Not available |
| Open-Meteo IFS | **9 km** | 9 km (post-2022 only) |

For snow/hydrology applications, 9km precipitation is a **big win** (4x resolution).

### Proposed API

```python
from nwp_downloader.era5 import ERA5Loader

loader = ERA5Loader(
    bbox=(7.7, 45.95, 7.85, 46.05),
    start_date="2020-01-01",
    end_date="2020-12-31",

    # Variable-specific source override
    source_map={
        "default": "google",           # Most variables from ARCO-ERA5
        "tp": "openmeteo_era5land",    # Precipitation from 9km ERA5-Land
    },
)
```

### Implementation Considerations

1. **Coordinate alignment:** Different sources may have different grids. Need interpolation or nearest-neighbor to common grid.
2. **Time alignment:** Ensure all sources have same time steps.
3. **Efficiency:** Batch requests per-source, don't make separate request per variable.
4. **Fallback:** If preferred source unavailable for a variable, fall back to default.

**[DECIDED]:** Option (c) - Return separate datasets, let consumer handle alignment.

**Rationale:** TPS2 downscales to point locations (unit centroids), not grids. Each variable is independently interpolated to the same target points. Grids never need to match - they just get sampled at the same locations. This preserves the 9km precip resolution benefit without requiring grid manipulation in the downloader.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        ERA5Loader                               │
│  (orchestrator: date iteration, parallelism, derived vars)      │
└─────────────────────────┬───────────────────────────────────────┘
                          │
        ┌─────────────────┼─────────────────┐
        │                 │                 │
        ▼                 ▼                 ▼
┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│ GoogleBackend│  │  CDSBackend  │  │OpenMeteoBackend│ ...
│  (ARCO Zarr) │  │  (CDS API)   │  │  (REST API)  │
└──────────────┘  └──────────────┘  └──────────────┘
        │                 │                 │
        └─────────────────┼─────────────────┘
                          │
                          ▼
                   ┌─────────────┐
                   │   Writer    │
                   │ NetCDF/Zarr │
                   └─────────────┘
```

### Components

| Component | Responsibility |
|-----------|---------------|
| `ERA5Loader` | Orchestration: build date list, parallel fetch, derived variable computation, write coordination |
| `ERA5Backend` (abstract) | Fetch one day of data from a specific source, return standardized `(ds_surf, ds_plev)` |
| `Writer` | Persist data to disk (NetCDF daily/yearly or Zarr append) |
| `IFSForecastLoader` | Separate loader for IFS open data forecasts (different data structure, GRIB input) |

---

## Backends

| Backend | Source | Auth | Latency | Historical Coverage |
|---------|--------|------|---------|---------------------|
| `google` | Google ARCO-ERA5 (gs://gcp-public-data-arco-era5) | Anonymous | ~20s/day | 1940-present (5-day delay) |
| `cds` | Copernicus Climate Data Store | API key | ~2-5min/day | 1940-present |
| `s3zarr` | AWS S3 ERA5 Zarr mirror | AWS creds or anonymous | ~20s/day | 1940-present |
| `openmeteo` | Open-Meteo API | Anonymous | ~5s/day | 2022-present (limited) |

### Backend Selection Logic (`select_backend`)

```
If pressure levels needed AND all dates >= 2022:
    openmeteo(ifs) → google → cds
If pressure levels needed AND any date < 2022:
    google → cds
If no pressure levels needed:
    openmeteo(era5) → google → cds
```

**[QUESTION 1]:** Is this priority order correct? Should `google` always be preferred over `openmeteo` for reliability, even if slower?

---

## Output Conventions

### Variable Names (ERA5 short names)

**Surface:**
| Variable | Name | Units |
|----------|------|-------|
| 2m temperature | `t2m` | K |
| 2m dewpoint | `d2m` | K |
| Surface pressure | `sp` | Pa |
| Surface geopotential | `z` | m²/s² |
| Downward SW radiation | `ssrd` | J/m² (accumulated per timestep) |
| Downward LW radiation | `strd` | J/m² (accumulated per timestep) |
| Total precipitation | `tp` | m (accumulated per timestep) |
| 10m U wind | `u10` | m/s |
| 10m V wind | `v10` | m/s |

**Pressure levels:**
| Variable | Name | Units |
|----------|------|-------|
| Geopotential | `z` | m²/s² |
| Temperature | `t` | K |
| U wind | `u` | m/s |
| V wind | `v` | m/s |
| Specific humidity | `q` | kg/kg |
| Relative humidity | `r` | % (0-100) |

### Dimensions

```
time: datetime64[ns]
latitude: float (descending, i.e. north to south)
longitude: float (-180 to 180)
level: int (hPa, ascending, i.e. 300, 500, 700, 850, 1000)
```

**[DECIDED]:** Latitude always descending (north-first). This is the ERA5/ECMWF convention and what TPS2 expects.

---

## Data Flow & Formats

### Input Formats (per backend)

| Backend | Native Format | Notes |
|---------|---------------|-------|
| Google ARCO-ERA5 | Zarr | Cloud-native, efficient |
| CDS | NetCDF | Downloaded per-day |
| S3 Zarr | Zarr | Cloud-native |
| ECMWF OpenData (IFS) | GRIB | Requires cfgrib |
| OpenMeteo | REST/JSON | Converted to xarray in-memory |

### Output Format: Always Zarr

**TPS2 expects Zarr.** nwp-downloader converts all sources to a unified local Zarr cache:

```
{cache_dir}/
  era5.zarr/           # ERA5 reanalysis
    time/
    latitude/
    longitude/
    level/
    t2m/
    tp/
    ...
  forecast.zarr/       # IFS forecasts (if used)
    ...
```

**Benefits of Zarr-only output:**
- Single format for TPS2 to handle
- Efficient append (add days without rewriting)
- Lazy loading (don't load 50GB into memory)
- Chunked access (read only needed bbox/time)
- Cloud-native (could push to cloud storage)

### Conversion Pipeline

```
Google Zarr ──────────────────────────────→ Local Zarr (chunk copy)
CDS NetCDF ───→ xarray ───→ rechunk ──────→ Local Zarr
ECMWF GRIB ───→ cfgrib ───→ xarray ───────→ Local Zarr
OpenMeteo ────→ REST ─────→ xarray ───────→ Local Zarr
```

**[DECIDED]:** No more NetCDF output. No daily+yearly duplication. Single Zarr store per data type.

### Incremental Downloads

```python
# First run: downloads 2020
loader.download(time_range=["2020-01-01", "2020-12-31"])

# Later: extend to 2021 (only fetches new data)
loader.download(time_range=["2020-01-01", "2021-12-31"])
# Detects existing 2020 data, only downloads 2021
```

Essential for large (50-100GB) operational datasets.

---

## Derived Variables

Computed if missing from source:

| Variable | Computed From | When |
|----------|--------------|------|
| Relative humidity (`r`) | `q`, `t`, `level` | If `compute_rh=True` and `r` missing |
| Surface geopotential (`z`) | `sp`, `msl`, `t2m` | IFS forecasts only (no native z) |

**[DECIDED]:** Derived variables computed in downloader. The downloader should deliver data ready-to-use. A preprocessor module handles all transformations needed to get raw downloads into simulation-ready state.

---

## Preprocessor Module

The `preprocess` module ensures all downloaded data is simulation-ready before being written to the final Zarr cache. This runs automatically after download.

### Pipeline

```
Raw Download → Validate → Derive → Conventions → Chunk → Ready Zarr
```

### Components

| Module | Responsibility |
|--------|---------------|
| `validate.py` | Range checks (T: 180-350K, P≥0, RH: 0-100%), missing data detection, temporal continuity |
| `derived.py` | Compute missing variables: strd from (t2m, d2m, cloud), q from (r, t, p) |
| `conventions.py` | Ensure lat descending, lon -180:180, standard units (K, Pa, J/m², m) |
| `chunking.py` | Rechunk for TPS2 access pattern (time-major for point extraction) |

### Derived Variable Computation

| Variable | Computed From | Method |
|----------|---------------|--------|
| `strd` | t2m, d2m, (cloud) | Dilley & O'Brien (1998) or Konzelmann et al. (1994) |
| `q` | r, t, p | Clausius-Clapeyron: q = 0.622 * e / (p - 0.378*e) |
| `r` | q, t, p | Inverse of above |

### When Preprocessing Runs

- **On download:** Each day is preprocessed before appending to Zarr
- **On open:** Quick validation check (time range, variables present)
- **Manual:** `nwp_downloader preprocess --validate-only` CLI command

---

## IFSForecastLoader

Separate class for ECMWF IFS open data forecasts:

- Downloads via `ecmwf.opendata` client
- GRIB format → requires cfgrib
- 10-day forecast (steps 0-240h)
- Deaccumulates/reaccumulates radiation and precip
- Interpolates from 3h/6h native resolution to 1h

**[QUESTION 6]:** Should `IFSForecastLoader` follow the same Backend/Writer pattern as `ERA5Loader`, or is the separate class approach correct given the different data structure?

---

## What ecmwf-downloader Does NOT Do

1. **Spatial interpolation** - returns data on native ERA5 grid (0.25°), consumer handles interpolation to points/clusters
2. **Downscaling** - no lapse rate correction, terrain effects, etc.
3. **Quality control** - no outlier detection, gap filling
4. **Unit conversion** - outputs in native ERA5 units (K, Pa, J/m², m)

**[QUESTION 7]:** Is this boundary correct? Are there any operations that should move into or out of ecmwf-downloader?

---

## Proposed API

### Size-Aware Loading

```python
from nwp_downloader.era5 import ERA5Loader

loader = ERA5Loader(
    bbox=(7.7, 45.95, 7.85, 46.05),
    time_range=["2020-01-01", "2020-12-31"],
    pressure_levels=[700, 850, 1000],
    cache_dir="./cache/",
    source_map={
        "default": "google",
        "tp": "openmeteo_era5land",  # 9km precip
    },
)

# Smart open: downloads if missing, returns lazy xarray handle
ds = loader.open()  # Returns lazy Dataset pointing to local Zarr

# For small jobs (< ~1GB): can request in-memory
ds = loader.fetch()  # Returns loaded Dataset (fits in memory)
```

### Behavior by Job Size

| Estimated Size | `open()` behavior | `fetch()` behavior |
|----------------|-------------------|-------------------|
| < 100 MB | Cache to Zarr, return lazy | Fetch to memory directly |
| 100 MB - 1 GB | Cache to Zarr, return lazy | Cache to Zarr, load to memory |
| > 1 GB | Cache to Zarr, return lazy | Raises error (use `open()`) |

### Integration with TPS2

```python
# topopyscale2/inputs/era5.py
class ERA5Source:
    def fetch(self, bbox, time_range) -> xr.Dataset:
        loader = ERA5Loader(
            bbox=bbox,
            time_range=time_range,
            cache_dir=self.config.cache_dir,
            source_map=self.config.source_map,
        )
        # Returns lazy handle to local Zarr cache
        # Incremental: only downloads missing data
        return loader.open()
```

TPS2 then:
1. Iterates over spatial units
2. For each unit, selects nearest grid cell(s) from lazy dataset
3. Loads only that small subset into memory
4. Processes and writes output

**Result:** 50GB ERA5 dataset never fully loaded into memory.

---

## Future Considerations

1. **Additional backends:** ICON, GFS, COSMO, HRES deterministic
2. **Ensemble support:** ERA5 EDA (10 members)
3. **Caching layer:** Avoid re-downloading same bbox/time
4. **Async/streaming:** For operational forecasting pipelines

**[QUESTION 9]:** Should the package be renamed to something more general (e.g., `meteo-downloader`, `nwp-downloader`) if we add non-ECMWF sources?

---

## Summary of Open Questions

| # | Question | Status |
|---|----------|--------|
| 1 | Backend priority order: google vs openmeteo? | **Decided: See `priorities-nwp.md`** |
| 2 | Latitude order: always descending or configurable? | **Decided: always descending** |
| 3 | Daily files: delete after yearly merge? | **Decided: No daily files. Single Zarr store.** |
| 4 | Zarr pass-through mode (no local cache)? | **Decided: Local Zarr cache is default. Direct cloud read possible but not primary mode.** |
| 5 | Derived variables: here or in consumer? | **Decided: in downloader (preprocess module)** |
| 6 | IFSForecastLoader: separate class or unified pattern? | **Decided: separate module** (`nwp_downloader.forecast`) |
| 7 | Scope boundary: any operations to add/remove? | Open |
| 8 | Streaming API (return dataset, no disk write)? | **Decided: Hybrid. `fetch()` for small jobs, `open()` returns lazy handle to local Zarr for large jobs.** |
| 9 | Rename package for non-ECMWF sources? | **Decided: yes, `nwp-downloader`** |
| 10 | Variable-dependent sources: how to handle coordinate alignment? | **Decided: (c) separate datasets** |

## Decisions Made

1. **Package name:** `nwp-downloader` (not ECMWF-specific)
2. **Module structure:** One module per domain (`era5`, `forecast`, `cmip6`)
3. **Variable-dependent sources:** Supported (e.g., 9km precip from Open-Meteo)
4. **Coordinate alignment:** None - return separate datasets per source, consumer (TPS2) handles point-based interpolation
5. **Output format:** Always Zarr (convert NetCDF/GRIB sources to Zarr)
6. **No duplication:** Single Zarr store, no daily+yearly files
7. **Size-aware API:** `fetch()` for small jobs (memory), `open()` for large jobs (lazy Zarr handle)
8. **Incremental downloads:** Only fetch missing time chunks, support resume
9. **Source priorities:** Documented in `priorities-nwp.md` (code should implement these rules)
