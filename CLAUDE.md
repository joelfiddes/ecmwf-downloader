# CLAUDE.md

## Project Overview

`ecmwf-downloader` is a standalone pip-installable package for downloading ECMWF climate data (ERA5 reanalysis + IFS forecasts) from multiple sources with flexible output formats. No CDO dependency — pure Python with xarray/cfgrib.

## Development Commands

```bash
# Install in development mode
pip install -e ".[all,dev]"

# Run tests
pytest tests/

# Lint
flake8 src/ --select=E9,F63,F7,F82 --show-source
```

## Architecture

### ERA5 Pipeline
- **Backends** (`backends/`): Each implements `fetch_day(date) -> (ds_surf, ds_plev)` with standardized names/dims
  - `cds.py`: CDS API backend (requires ~/.cdsapirc)
  - `google.py`: Google ARCO-ERA5 netCDF (anonymous, no auth needed)
  - `s3zarr.py`: S3-compatible zarr store
- **Writers** (`writer.py`): NetCDFWriter (daily files + yearly merge) and ZarrWriter (single store)
- **Orchestrator** (`loader.py`): ERA5Loader wires backend + writer + RH computation

### IFS Forecast Pipeline
- **IFSForecastLoader** (`ifs_forecast.py`): Downloads ECMWF open data IFS forecasts, processes with cfgrib/xarray (no CDO)

### Output Conventions
All output uses ERA5 native conventions:
- Radiation (ssrd, strd): accumulated per timestep (J/m²)
- Precipitation (tp): accumulated per timestep (m)
- Standard dims: time, latitude, longitude, level (hPa ascending)
