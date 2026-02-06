# TopoPyScale Flexible Data Acquisition — Design Sketch

## Overview

A refactored data acquisition and downscaling system that dynamically selects
optimal data sources based on variable availability, resolution, access speed,
and request limits. Replaces the current rigid CDS-only pipeline.

---

## 1. Configuration

The user provides a config (YAML or dict) that defines **what** they need.
The system figures out **where** to get it and **how** to downscale it.

```yaml
project:
  name: "northern_kazakhstan"
  time:
    start: "2020-10-01"
    end: "2021-06-30"
    freq: "1h"

  # Domain is DERIVED from points/clusters, not specified as bbox
  # Points come from TopoSUB clusters or user-supplied coords
  points_file: "cluster_centroids.csv"  # lon, lat, elevation
  # OR
  points:
    - {lon: 68.5, lat: 51.2, elevation: 450}
    - {lon: 69.1, lat: 52.0, elevation: 820}

  # Variables requested — user just says what they need
  variables:
    - temperature
    - precipitation
    - shortwave_radiation
    - longwave_radiation
    - wind
    - humidity
    - surface_pressure

  # Pressure levels: "auto" derives from max domain elevation
  # or user can specify explicitly
  pressure_levels: auto   # or [1000, 925, 850, 700, 600, 500]

  # Source priority (user can reorder or exclude)
  source_priority:
    - open_meteo_api
    - open_meteo_s3
    - arco_era5
    - ecmwf_ifs_api
    - cds

  # Downscaling mode
  downscaling: auto  # or "full_physics", "lapse_rate_only", "none"
```

---

## 2. Domain Analysis (from points)

Before any data retrieval, the system analyses the point set to determine:

```
┌─────────────────────────────────────────────────────────┐
│                    DOMAIN ANALYSER                       │
│                                                         │
│  Input: list of (lon, lat, elevation)                   │
│                                                         │
│  Derives:                                               │
│  ├── bbox (with padding for 3×3 interpolation stencil)  │
│  ├── max_elevation → required pressure levels           │
│  │   e.g. 3000m → need ≥700 hPa                        │
│  │        5000m → need ≥500 hPa                         │
│  │        1000m → need ≥850 hPa                         │
│  ├── point_density → retrieval strategy                 │
│  │   (sparse points / clustered / dense)                │
│  ├── unique_grid_cells (per source resolution)          │
│  └── spatial_extent_km²                                 │
│                                                         │
│  Output: DomainSpec                                     │
│  ├── .bbox                                              │
│  ├── .pressure_levels   [1000, 925, 850, 700, 600, 500] │
│  ├── .n_points          2000                            │
│  ├── .density_class     "dense"                         │
│  ├── .extent_km2        450000                          │
│  └── .grid_cells        {era5_025: 1800, ifs_01: 12000} │
└─────────────────────────────────────────────────────────┘
```

### Pressure level selection from elevation

```
max_elev (m)  →  min_pressure (hPa)  →  levels needed
   500            950                    [1000, 975, 950]
  1000            900                    [1000, 975, 950, 925, 900]
  2000            800                    [1000, 925, 850, 800]
  3000            700                    [1000, 925, 850, 800, 700]
  4000            600                    [1000, 925, 850, 700, 600]
  5000            500                    [1000, 925, 850, 700, 600, 500]

Rule of thumb: min_pressure ≈ 1013 * exp(-max_elev / 7400) - margin
Always include at least 2 levels above max elevation for interpolation
```

---

## 3. Variable Registry

Central registry that maps abstract variable names to source-specific
parameter names and tracks what each source can provide.

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        VARIABLE REGISTRY                                │
│                                                                         │
│  Variable         │ Level     │ Downscaling method (by data available)  │
│  ─────────────────┼───────────┼─────────────────────────────────────────│
│  temperature      │ pressure  │ vertical profile interpolation          │
│                   │ surface   │ lapse rate correction (fallback)        │
│  ─────────────────┼───────────┼─────────────────────────────────────────│
│  precipitation    │ surface   │ IDW + optional elevation scaling        │
│  ─────────────────┼───────────┼─────────────────────────────────────────│
│  sw_radiation     │ surface   │ slope/aspect/horizon (if direct+diffuse)│
│                   │           │ simple aspect correction (if GHI only)  │
│  ─────────────────┼───────────┼─────────────────────────────────────────│
│  lw_radiation     │ surface   │ elevation + humidity correction         │
│                   │           │ parameterised from T+RH (fallback)      │
│  ─────────────────┼───────────┼─────────────────────────────────────────│
│  wind_speed       │ pressure  │ vertical profile + terrain exposure     │
│                   │ surface   │ terrain exposure correction (fallback)  │
│  ─────────────────┼───────────┼─────────────────────────────────────────│
│  wind_direction   │ pressure  │ vertical interpolation                  │
│                   │ surface   │ passthrough                             │
│  ─────────────────┼───────────┼─────────────────────────────────────────│
│  humidity         │ pressure  │ specific humidity profile interpolation │
│  (spec. hum / RH) │ surface   │ convert RH↔Q with T,P (fallback)       │
│  ─────────────────┼───────────┼─────────────────────────────────────────│
│  surface_pressure │ surface   │ barometric correction to elevation      │
│  ─────────────────┼───────────┼─────────────────────────────────────────│
│  geopotential     │ pressure  │ (needed for vertical interpolation)     │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 4. Source Capabilities Matrix

Each data source declares what it can provide. This is a static registry
(updated when sources change) that the planner queries.

```
┌──────────────────────────────────────────────────────────────────────────────────────┐
│                           SOURCE CAPABILITIES                                        │
│                                                                                      │
│  Source            │ Res    │ Plev │ LW  │ Q   │ Access   │ Limits      │ Period      │
│  ──────────────────┼────────┼──────┼─────┼─────┼──────────┼─────────────┼─────────────│
│  open_meteo_api    │ 9km*   │  ❌  │ ❌  │ ❌  │ REST     │ 10k/day free│ 1940+       │
│                    │        │      │     │     │ (points) │ €30/mo paid │             │
│  ──────────────────┼────────┼──────┼─────┼─────┼──────────┼─────────────┼─────────────│
│  open_meteo_s3     │ 9km*   │  ❌  │ ❌  │ ❌  │ S3 chunk │ unlimited   │ 1940+       │
│  (OM-files)        │        │      │     │     │ (spatial)│             │             │
│  ──────────────────┼────────┼──────┼─────┼─────┼──────────┼─────────────┼─────────────│
│  arco_era5         │ 0.25°  │  ✅  │ ✅  │ ✅  │ Zarr/GCS │ unlimited   │ 1940+       │
│  (Google Cloud)    │        │ 37lv │     │     │ (chunks) │ (egress $$) │             │
│  ──────────────────┼────────┼──────┼─────┼─────┼──────────┼─────────────┼─────────────│
│  ecmwf_ifs_api     │ 9km    │  ?   │ ?   │ ?   │ REST/API │ varies      │ 2017+       │
│  (new ECMWF API)   │        │      │     │     │          │             │             │
│  ──────────────────┼────────┼──────┼─────┼─────┼──────────┼─────────────┼─────────────│
│  cds               │ 0.25°  │  ✅  │ ✅  │ ✅  │ Queue    │ slow but    │ 1940+       │
│  (Copernicus)      │        │ 37lv │     │     │ (bbox)   │ unlimited   │             │
│  ──────────────────┼────────┼──────┼─────┼─────┼──────────┼─────────────┼─────────────│
│                                                                                      │
│  * "9km" = blended IFS (2017+) / ERA5-Land (0.1°) / ERA5 (0.25°)                    │
│                                                                                      │
│  Key: Plev = pressure levels, LW = longwave radiation, Q = specific humidity         │
└──────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 5. Source Selection Algorithm

The planner walks the priority list and tries to satisfy each requested
variable from the highest-priority source. Variables that can't be
fulfilled cascade to the next source.

```
┌──────────────────────────────────────────────────────────────────────┐
│                       SOURCE PLANNER                                 │
│                                                                      │
│  Input:                                                              │
│    - requested_variables: [T, P, SW, LW, wind, humidity, Ps]        │
│    - domain_spec: (bbox, pressure_levels, n_points, density...)     │
│    - source_priority: [open_meteo_api, open_meteo_s3, arco, cds]    │
│                                                                      │
│  Algorithm:                                                          │
│                                                                      │
│    remaining = set(requested_variables)                               │
│    plan = {}                                                         │
│                                                                      │
│    for source in source_priority:                                    │
│        capabilities = source.can_provide(remaining, domain_spec)     │
│        # capabilities checks:                                        │
│        #   - does source have this variable?                         │
│        #   - at required pressure levels?                            │
│        #   - for requested time period?                              │
│        #   - within rate limits for this domain size?                │
│                                                                      │
│        satisfiable = remaining ∩ capabilities                        │
│        if satisfiable:                                               │
│            plan[source] = satisfiable                                │
│            remaining -= satisfiable                                  │
│                                                                      │
│    if remaining:                                                     │
│        # Some variables unfulfilled                                  │
│        for var in remaining:                                         │
│            if has_parameterisation(var):                              │
│                # e.g. LW can be parameterised from T + RH            │
│                plan['parameterised'].add(var)                        │
│            else:                                                     │
│                warn(f"Cannot obtain {var} from any source")          │
│                                                                      │
│    return RetrievalPlan(plan)                                        │
│                                                                      │
│  Example output for full-physics run:                                │
│    open_meteo_api  →  [T_surface, P_precip, SW_rad, wind_10m,       │
│                        humidity_2m, surface_pressure]                │
│    arco_era5       →  [T_plev, wind_plev, Q_plev, geopotential,     │
│                        LW_radiation]                                 │
│                                                                      │
│  Example output for T+P only run:                                    │
│    open_meteo_api  →  [T_surface, P_precip]                         │
│    (nothing else needed)                                             │
│                                                                      │
└──────────────────────────────────────────────────────────────────────┘
```

---

## 6. Retrieval Strategy (per source)

Once the planner knows WHICH source provides WHICH variables, the
retriever decides HOW to fetch from each source efficiently.

```
┌──────────────────────────────────────────────────────────────────────┐
│                    RETRIEVAL STRATEGY SELECTOR                       │
│                                                                      │
│  Input: source, domain_spec                                          │
│                                                                      │
│  Step 1: Compute unique grid cells for this source's resolution      │
│                                                                      │
│    For each point, identify the 3×3 stencil of grid cells            │
│    needed for IDW interpolation. Merge overlapping stencils.         │
│    → unique_cells                                                    │
│                                                                      │
│  Step 2: Compare to bbox                                             │
│                                                                      │
│    bbox_cells = cells covering padded bounding box                   │
│    fill_ratio = len(unique_cells) / bbox_cells                       │
│                                                                      │
│  Step 3: Select strategy                                             │
│                                                                      │
│    ┌─────────────────────────────────────────────────────┐           │
│    │  fill_ratio > 0.5                                   │           │
│    │  OR source prefers bbox (CDS, Zarr)                 │           │
│    │  → BBOX retrieval                                   │           │
│    ├─────────────────────────────────────────────────────┤           │
│    │  fill_ratio < 0.1                                   │           │
│    │  AND source supports point queries (Open-Meteo API) │           │
│    │  → POINT retrieval (unique cells only)              │           │
│    ├─────────────────────────────────────────────────────┤           │
│    │  fill_ratio 0.1–0.5                                 │           │
│    │  Points form spatial clusters                       │           │
│    │  → SUB-BBOX retrieval (DBSCAN cluster → mini boxes) │           │
│    └─────────────────────────────────────────────────────┘           │
│                                                                      │
│  Step 4: Estimate cost                                               │
│                                                                      │
│    n_cells × n_variables × n_timesteps × source_cost_factor          │
│    → estimated API calls / download size / time                      │
│                                                                      │
└──────────────────────────────────────────────────────────────────────┘

Visual example — same 10 points, different strategies:

  BBOX (wasteful for sparse points)    POINT-BASED (efficient)
  ┌──────────────────────────┐         Only fetch what's needed:
  │ · · · · · · · · · · · · │
  │ · · · · · · · · · ■ · · │              ■ ■ ■
  │ · · · · · · · · · · · · │              ■ ■ ■    ← 3×3 patch
  │ · ■ · · · · · · · · · · │              ■ ■ ■
  │ · · · · · · · · · · · · │
  │ · · · · · · · · · · · · │         ■ ■ ■
  │ · · · · · ■ · · · · · · │         ■ ■ ■         (etc.)
  │ · · · · · · · · · · · · │         ■ ■ ■
  └──────────────────────────┘
  Downloads ~500 cells                 Downloads ~90 cells
```

---

## 7. Downscaling Method Selection

The downscaler adapts based on what data was actually retrieved.

```
┌──────────────────────────────────────────────────────────────────────┐
│                    DOWNSCALING DISPATCHER                            │
│                                                                      │
│  For each variable, select method based on available data:           │
│                                                                      │
│  TEMPERATURE                                                         │
│  ├── Have pressure level profiles + geopotential?                    │
│  │   YES → vertical_profile_interpolation()                          │
│  │         Interpolate T through atmospheric column to point elev    │
│  │         ✓ Best accuracy, captures inversions                      │
│  └── Only surface T?                                                 │
│      YES → lapse_rate_correction()                                   │
│            T_local = T_grid + lapse_rate × (elev_local - elev_grid)  │
│            Default: -6.5°C/km, or month-varying, or from config     │
│            ✓ Good enough for most hydrology applications             │
│                                                                      │
│  PRECIPITATION                                                       │
│  ├── Always surface variable                                         │
│  └── IDW from 3×3 stencil + optional elevation scaling               │
│                                                                      │
│  SHORTWAVE RADIATION                                                 │
│  ├── Have direct + diffuse components?                               │
│  │   YES → full_radiation_partitioning()                             │
│  │         Separate direct beam (slope/aspect/horizon dependent)     │
│  │         from diffuse (sky view factor dependent)                  │
│  └── Only GHI (global horizontal)?                                   │
│      YES → simple_terrain_correction()                               │
│            Apply cos(incidence_angle) scaling                        │
│            ✓ Approximate but reasonable                              │
│                                                                      │
│  LONGWAVE RADIATION                                                  │
│  ├── Have LW from source?                                            │
│  │   YES → elevation_correction()                                    │
│  │         Adjust for T and humidity at local elevation              │
│  └── No LW available?                                                │
│      YES → parameterise_from_T_RH()                                  │
│            Use Stefan-Boltzmann + emissivity model                   │
│            e.g., Dilley & O'Brien (1998), Prata (1996)              │
│            ✓ ~10-20 W/m² RMSE typical, adequate for many uses       │
│                                                                      │
│  WIND                                                                │
│  ├── Have pressure level wind profiles?                              │
│  │   YES → vertical_interpolation + terrain_exposure()               │
│  └── Only 10m wind?                                                  │
│      YES → terrain_exposure_correction()                             │
│            Wind exposure index from DEM                              │
│                                                                      │
│  HUMIDITY                                                            │
│  ├── Have specific humidity on pressure levels?                      │
│  │   YES → vertical_interpolation()                                  │
│  └── Only surface RH or dewpoint?                                    │
│      YES → convert_and_adjust()                                      │
│            RH → Q using T and P at local elevation                   │
│                                                                      │
└──────────────────────────────────────────────────────────────────────┘
```

---

## 8. Overall Pipeline Flow

```
USER CONFIG
    │
    ▼
┌─────────────┐     ┌──────────────────┐
│ Parse config │────▶│  Domain Analyser │
│ (YAML/dict)  │     │                  │
└─────────────┘     │  - bbox from pts │
                    │  - max elevation │
                    │  - pressure lvls │
                    │  - point density │
                    └────────┬─────────┘
                             │
                             ▼
                    ┌──────────────────┐
                    │  Source Planner   │
                    │                  │
                    │  Walk priority   │
                    │  list, match     │
                    │  variables to    │
                    │  best sources    │
                    └────────┬─────────┘
                             │
                             ▼
                    ┌──────────────────┐
                    │  Retrieval Plan   │
                    │                  │
                    │  source A → vars │     ┌─────────────────┐
                    │  source B → vars │────▶│ Retrieval Engine │
                    │  parameterised   │     │                  │
                    └──────────────────┘     │  Per source:     │
                                            │  - select bbox/  │
                                            │    point/sub-bbox│
                                            │  - fetch data    │
                                            │  - cache locally │
                                            └────────┬─────────┘
                                                     │
                                                     ▼
                                            ┌──────────────────┐
                                            │ Raw Data Store    │
                                            │ (unified xarray)  │
                                            │                   │
                                            │ Surface vars from │
                                            │ source A (9km)    │
                                            │ Plev vars from    │
                                            │ source B (0.25°)  │
                                            └────────┬──────────┘
                                                     │
                                                     ▼
                                            ┌──────────────────┐
                                            │ Spatial Interp    │
                                            │                   │
                                            │ IDW from 3×3      │
                                            │ stencil to each   │
                                            │ point/cluster     │
                                            └────────┬──────────┘
                                                     │
                                                     ▼
                                            ┌──────────────────┐
                                            │ Downscale         │
                                            │ Dispatcher        │
                                            │                   │
                                            │ Per variable:     │
                                            │ select method     │
                                            │ based on what     │
                                            │ data is available │
                                            └────────┬──────────┘
                                                     │
                                                     ▼
                                            ┌──────────────────┐
                                            │ Output            │
                                            │                   │
                                            │ Per point/cluster │
                                            │ timeseries:       │
                                            │ - FSM format      │
                                            │ - CSV             │
                                            │ - NetCDF          │
                                            │ - CryoGrid        │
                                            │ - SNOWPACK        │
                                            └──────────────────┘
```

---

## 9. Key Design Principles

1. **Variables are first-class citizens** — each variable independently
   selects its source and downscaling method. No more "all or nothing."

2. **Graceful degradation** — if pressure levels aren't available, fall
   back to lapse rates. If LW isn't available, parameterise it. The user
   gets the best possible result from whatever data is accessible.

3. **Source-aware retrieval** — the system knows each source's access
   pattern (point queries, bbox, chunked) and optimises accordingly.

4. **Separation of concerns**:
   - Config: what do I need?
   - Planner: where do I get it?
   - Retriever: how do I fetch it efficiently?
   - Downscaler: how do I adjust it for local terrain?

5. **Transparent** — the system should log/report what it did:
   "Temperature from Open-Meteo (lapse rate correction),
    LW parameterised from T+RH (Prata 1996),
    Geopotential from ARCO-ERA5 (vertical interpolation)"

6. **Points are the primitive** — everything works on (lon, lat, elev)
   tuples. Whether those come from TopoSUB clusters, a user list, or
   a regular grid doesn't matter to the retrieval/downscaling system.
