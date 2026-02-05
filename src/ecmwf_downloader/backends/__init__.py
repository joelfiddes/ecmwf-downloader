"""Backend registry for ERA5 data sources."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ecmwf_downloader.base import ERA5Backend

BACKEND_REGISTRY: dict[str, str] = {
    "cds": "ecmwf_downloader.backends.cds:CDSBackend",
    "google": "ecmwf_downloader.backends.google:GoogleCloudBackend",
    "s3zarr": "ecmwf_downloader.backends.s3zarr:S3ZarrBackend",
}


def get_backend(name: str, **kwargs) -> "ERA5Backend":
    """Instantiate a backend by name.

    Args:
        name: One of 'cds', 'google', 's3zarr'.
        **kwargs: Passed to the backend constructor.

    Returns:
        An ERA5Backend instance.
    """
    if name not in BACKEND_REGISTRY:
        raise ValueError(
            f"Unknown backend '{name}'. Available: {list(BACKEND_REGISTRY.keys())}"
        )

    module_path, class_name = BACKEND_REGISTRY[name].rsplit(":", 1)
    import importlib

    module = importlib.import_module(module_path)
    cls = getattr(module, class_name)
    return cls(**kwargs)
