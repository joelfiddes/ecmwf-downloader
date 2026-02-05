"""Backend registry for ERA5 data sources."""

from __future__ import annotations

import importlib
import logging
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from ecmwf_downloader.base import ERA5Backend

logger = logging.getLogger(__name__)

BACKEND_REGISTRY: dict[str, str] = {
    "cds": "ecmwf_downloader.backends.cds:CDSBackend",
    "google": "ecmwf_downloader.backends.google:GoogleCloudBackend",
    "s3zarr": "ecmwf_downloader.backends.s3zarr:S3ZarrBackend",
    "openmeteo": "ecmwf_downloader.backends.openmeteo:OpenMeteoBackend",
}

# Key dependency for each backend used by the availability check.
_BACKEND_DEPS: dict[str, str] = {
    "openmeteo": "requests",
    "google": "gcsfs",
    "cds": "cdsapi",
}


def _check_available(backend_name: str) -> bool:
    """Return True if the key dependency for *backend_name* can be imported."""
    dep = _BACKEND_DEPS.get(backend_name)
    if dep is None:
        return False
    try:
        importlib.import_module(dep)
        return True
    except ImportError:
        return False


def select_backend(
    start_date: pd.Timestamp,
    end_date: pd.Timestamp,
    pressure_levels: list[int] | None,
) -> tuple[str, dict]:
    """Choose the fastest available backend for the given parameters.

    Priority logic:
        plev needed + all post-2022 → openmeteo(ifs) → google → cds
        plev needed + pre-2022      → google → cds
        no plev needed              → openmeteo(era5) → google → cds

    Returns:
        (backend_name, backend_kwargs) ready to pass to ``get_backend``.

    Raises:
        RuntimeError: If no backend dependency is installed.
    """
    plev_needed = bool(pressure_levels)
    all_post_2022 = start_date >= pd.Timestamp("2022-01-01")

    if plev_needed and all_post_2022:
        candidates = [
            ("openmeteo", {"model": "ifs", "start_date": str(start_date.date()), "end_date": str(end_date.date())}),
            ("google", {}),
            ("cds", {}),
        ]
    elif plev_needed:
        candidates = [
            ("google", {}),
            ("cds", {}),
        ]
    else:
        candidates = [
            ("openmeteo", {"model": "era5", "start_date": str(start_date.date()), "end_date": str(end_date.date())}),
            ("google", {}),
            ("cds", {}),
        ]

    for name, kwargs in candidates:
        if _check_available(name):
            logger.info("Auto-selected backend: %s", name)
            return name, kwargs

    tried = [name for name, _ in candidates]
    raise RuntimeError(
        f"No backend available. Tried: {tried}. "
        "Install one of: requests (openmeteo), gcsfs (google), cdsapi (cds)."
    )


def get_backend(name: str, **kwargs) -> "ERA5Backend":
    """Instantiate a backend by name.

    Args:
        name: One of 'cds', 'google', 's3zarr', 'openmeteo'.
        **kwargs: Passed to the backend constructor.

    Returns:
        An ERA5Backend instance.
    """
    if name not in BACKEND_REGISTRY:
        raise ValueError(
            f"Unknown backend '{name}'. Available: {list(BACKEND_REGISTRY.keys())}"
        )

    module_path, class_name = BACKEND_REGISTRY[name].rsplit(":", 1)
    module = importlib.import_module(module_path)
    cls = getattr(module, class_name)
    return cls(**kwargs)
