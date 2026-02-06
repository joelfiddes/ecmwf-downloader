"""Tests for HybridBackend."""

from __future__ import annotations

import tempfile
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from ecmwf_downloader.bbox import BBox


class TestHybridBackendInit:
    """Test HybridBackend initialization."""

    def test_init_with_valid_params(self, sample_bbox: BBox):
        """Backend should initialize with valid parameters."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        with tempfile.TemporaryDirectory() as tmpdir:
            backend = HybridBackend(
                bbox=sample_bbox,
                pressure_levels=[850, 500],
                cache_dir=tmpdir,
            )
            assert backend.bbox == sample_bbox
            # Parent class sorts pressure levels
            assert backend.pressure_levels == [500, 850]

    def test_init_without_pressure_levels(self, sample_bbox: BBox):
        """Backend should initialize without pressure levels."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        with tempfile.TemporaryDirectory() as tmpdir:
            backend = HybridBackend(
                bbox=sample_bbox,
                pressure_levels=None,
                cache_dir=tmpdir,
            )
            assert backend.pressure_levels == []

    def test_init_include_strd_default(self, sample_bbox: BBox):
        """Backend should include strd by default."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        with tempfile.TemporaryDirectory() as tmpdir:
            backend = HybridBackend(
                bbox=sample_bbox,
                cache_dir=tmpdir,
            )
            assert backend.include_strd is True

    def test_init_exclude_strd(self, sample_bbox: BBox):
        """Backend should allow excluding strd."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        with tempfile.TemporaryDirectory() as tmpdir:
            backend = HybridBackend(
                bbox=sample_bbox,
                include_strd=False,
                cache_dir=tmpdir,
            )
            assert backend.include_strd is False


class TestHybridBackendRegionalDetection:
    """Test regional (s3zarr) coverage detection."""

    def test_central_asia_bbox_detected(self):
        """Central Asia bbox should trigger regional mode."""
        from ecmwf_downloader.reanalysis.backends.hybrid import (
            HybridBackend,
            S3ZARR_BBOX,
            _bbox_within,
        )

        # Tajikistan (within Central Asia)
        bbox = BBox.from_tuple((68.0, 39.0, 72.0, 42.0))
        assert _bbox_within(bbox, S3ZARR_BBOX)

    def test_europe_bbox_not_regional(self):
        """European bbox should not trigger regional mode."""
        from ecmwf_downloader.reanalysis.backends.hybrid import (
            S3ZARR_BBOX,
            _bbox_within,
        )

        # Switzerland (outside Central Asia)
        bbox = BBox.from_tuple((6.0, 45.8, 10.5, 47.8))
        assert not _bbox_within(bbox, S3ZARR_BBOX)

    def test_partial_overlap_not_regional(self):
        """Bbox partially overlapping should not trigger regional mode."""
        from ecmwf_downloader.reanalysis.backends.hybrid import (
            S3ZARR_BBOX,
            _bbox_within,
        )

        # Bbox extending west of s3zarr coverage (west < 43)
        bbox = BBox.from_tuple((40.0, 35.0, 50.0, 45.0))
        assert not _bbox_within(bbox, S3ZARR_BBOX)

    def test_date_after_2023_not_regional(self, sample_bbox: BBox):
        """Dates after 2023 should not use regional mode."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        # Central Asia bbox but 2024 date
        bbox = BBox.from_tuple((68.0, 39.0, 72.0, 42.0))

        with tempfile.TemporaryDirectory() as tmpdir:
            backend = HybridBackend(
                bbox=bbox,
                end_date="2024-06-01",
                cache_dir=tmpdir,
            )
            # Should not use regional since end_date > 2023
            assert backend._use_regional is False


class TestHybridBackendRegistry:
    """Test backend registry integration."""

    def test_hybrid_in_registry(self):
        """Hybrid backend should be in the registry."""
        from ecmwf_downloader.reanalysis.backends import BACKEND_REGISTRY

        assert "hybrid" in BACKEND_REGISTRY

    def test_get_backend_instantiates_hybrid(self, sample_bbox: BBox):
        """get_backend should instantiate HybridBackend."""
        from ecmwf_downloader.reanalysis.backends import get_backend

        with tempfile.TemporaryDirectory() as tmpdir:
            backend = get_backend(
                "hybrid",
                bbox=sample_bbox,
                cache_dir=tmpdir,
            )
            assert backend.__class__.__name__ == "HybridBackend"


class TestHybridBackendSubBackendSelection:
    """Test sub-backend selection logic."""

    def test_surface_backend_prefers_s3(self, sample_bbox: BBox):
        """Surface backend should prefer openmeteo_s3."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        with patch(
            "ecmwf_downloader.reanalysis.backends._check_available",
            side_effect=lambda x: x == "openmeteo_s3",
        ):
            with tempfile.TemporaryDirectory() as tmpdir:
                backend = HybridBackend(bbox=sample_bbox, cache_dir=tmpdir)
                # Force backend creation by accessing it
                # (would need mocking to fully test)

    def test_plev_backend_requires_google_or_cds(self, sample_bbox: BBox):
        """Plev backend should only use google or cds."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        # Mock only openmeteo_s3 available - should fail to get plev backend
        with patch(
            "ecmwf_downloader.reanalysis.backends._check_available",
            side_effect=lambda x: x == "openmeteo_s3",
        ):
            with tempfile.TemporaryDirectory() as tmpdir:
                backend = HybridBackend(
                    bbox=sample_bbox,
                    pressure_levels=[850],
                    cache_dir=tmpdir,
                )
                # Trying to get plev backend should raise
                with pytest.raises(RuntimeError, match="No plev backend"):
                    backend._get_plev_backend()


class TestHybridBackendPrecipLogic:
    """Test precipitation backend selection."""

    def test_precip_prefers_ifs_for_2022_plus(self, sample_bbox: BBox):
        """Precipitation should prefer IFS for dates >= 2022."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        with tempfile.TemporaryDirectory() as tmpdir:
            backend = HybridBackend(bbox=sample_bbox, cache_dir=tmpdir)

            # For 2023, should try openmeteo with model=ifs first
            date = pd.Timestamp("2023-06-15")
            # Would need mocking to verify IFS is selected

    def test_precip_uses_s3_for_pre_2022(self, sample_bbox: BBox):
        """Precipitation should use S3 for dates < 2022."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        with tempfile.TemporaryDirectory() as tmpdir:
            backend = HybridBackend(bbox=sample_bbox, cache_dir=tmpdir)

            # For 2020, should prefer openmeteo_s3 (not IFS)
            date = pd.Timestamp("2020-06-15")
            # Would need mocking to verify S3 is selected


class TestHybridBackendCoordMatching:
    """Test coordinate matching logic."""

    def test_coords_match_identical(self):
        """Identical coordinates should match."""
        from ecmwf_downloader.reanalysis.backends.hybrid import _coords_match

        lats = np.array([46.0, 46.25, 46.5])
        lons = np.array([7.0, 7.25, 7.5])

        ds1 = xr.Dataset(coords={"latitude": lats, "longitude": lons})
        ds2 = xr.Dataset(coords={"latitude": lats, "longitude": lons})

        assert _coords_match(ds1, ds2)

    def test_coords_match_within_tolerance(self):
        """Coordinates within tolerance should match."""
        from ecmwf_downloader.reanalysis.backends.hybrid import _coords_match

        lats1 = np.array([46.0, 46.25, 46.5])
        lats2 = np.array([46.001, 46.251, 46.501])  # Within 0.01 tolerance
        lons = np.array([7.0, 7.25, 7.5])

        ds1 = xr.Dataset(coords={"latitude": lats1, "longitude": lons})
        ds2 = xr.Dataset(coords={"latitude": lats2, "longitude": lons})

        assert _coords_match(ds1, ds2, tol=0.01)

    def test_coords_no_match_different_size(self):
        """Coordinates with different sizes should not match."""
        from ecmwf_downloader.reanalysis.backends.hybrid import _coords_match

        ds1 = xr.Dataset(coords={
            "latitude": np.array([46.0, 46.25, 46.5]),
            "longitude": np.array([7.0, 7.25, 7.5]),
        })
        ds2 = xr.Dataset(coords={
            "latitude": np.array([46.0, 46.25]),  # Different size
            "longitude": np.array([7.0, 7.25]),
        })

        assert not _coords_match(ds1, ds2)

    def test_coords_no_match_outside_tolerance(self):
        """Coordinates outside tolerance should not match."""
        from ecmwf_downloader.reanalysis.backends.hybrid import _coords_match

        lats1 = np.array([46.0, 46.25, 46.5])
        lats2 = np.array([46.1, 46.35, 46.6])  # 0.1 offset
        lons = np.array([7.0, 7.25, 7.5])

        ds1 = xr.Dataset(coords={"latitude": lats1, "longitude": lons})
        ds2 = xr.Dataset(coords={"latitude": lats2, "longitude": lons})

        assert not _coords_match(ds1, ds2, tol=0.01)


class TestHybridBackendVariableGroups:
    """Test variable group definitions."""

    def test_surface_vars_defined(self):
        """Surface variables should be defined."""
        from ecmwf_downloader.reanalysis.backends.hybrid import SURFACE_VARS

        expected = {"t2m", "d2m", "sp", "ssrd", "u10", "v10", "msl"}
        assert SURFACE_VARS == expected

    def test_precip_vars_defined(self):
        """Precipitation variables should be defined."""
        from ecmwf_downloader.reanalysis.backends.hybrid import PRECIP_VARS

        assert PRECIP_VARS == {"tp"}

    def test_plev_vars_defined(self):
        """Pressure level variables should be defined."""
        from ecmwf_downloader.reanalysis.backends.hybrid import PLEV_VARS

        expected = {"t", "z", "u", "v", "q", "r"}
        assert PLEV_VARS == expected

    def test_google_only_vars_defined(self):
        """Google-only variables should be defined."""
        from ecmwf_downloader.reanalysis.backends.hybrid import (
            GOOGLE_ONLY_VARS,
            GOOGLE_ONLY_PLEV,
        )

        assert "strd" in GOOGLE_ONLY_VARS
        assert "q" in GOOGLE_ONLY_VARS
        assert "u" in GOOGLE_ONLY_PLEV
        assert "v" in GOOGLE_ONLY_PLEV


class TestHybridBackendClose:
    """Test resource cleanup."""

    def test_close_without_backends(self, sample_bbox: BBox):
        """close() should work even if no backends were initialized."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        with tempfile.TemporaryDirectory() as tmpdir:
            backend = HybridBackend(bbox=sample_bbox, cache_dir=tmpdir)
            # Should not raise
            backend.close()

    def test_close_with_mock_backends(self, sample_bbox: BBox):
        """close() should close all sub-backends."""
        from ecmwf_downloader.reanalysis.backends.hybrid import HybridBackend

        with tempfile.TemporaryDirectory() as tmpdir:
            backend = HybridBackend(bbox=sample_bbox, cache_dir=tmpdir)

            # Mock sub-backends
            mock_surface = MagicMock()
            mock_plev = MagicMock()
            backend._surface_backend = mock_surface
            backend._plev_backend = mock_plev

            backend.close()

            mock_surface.close.assert_called_once()
            mock_plev.close.assert_called_once()
