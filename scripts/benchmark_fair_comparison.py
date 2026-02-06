#!/usr/bin/env python
"""Fair comparison: Google-only vs Hybrid (S3 + Google).

Both fetch the same variables for the same region and time period.
"""

from __future__ import annotations

import sys
import tempfile
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ecmwf_downloader.bbox import BBox
from ecmwf_downloader.reanalysis.backends import get_backend


def benchmark_google_only(bbox, start_date, end_date, pressure_levels, cache_dir):
    """Fetch everything from Google ARCO-ERA5."""
    print("\n" + "=" * 70)
    print("BENCHMARK: Google-only (all variables from single source)")
    print("=" * 70)

    bbox_obj = BBox.from_tuple(bbox)
    dates = pd.date_range(start_date, end_date, freq="D")

    # Use explicit plev_vars to avoid relative_humidity issue
    backend = get_backend(
        "google",
        bbox=bbox_obj,
        pressure_levels=pressure_levels,
        cache_dir=cache_dir,
        plev_vars=[
            "temperature",
            "geopotential",
            "u_component_of_wind",
            "v_component_of_wind",
            "specific_humidity",
        ],
    )

    total_start = time.perf_counter()
    results = []

    for i, date in enumerate(dates):
        print(f"  [{i+1}/{len(dates)}] Fetching {date.strftime('%Y-%m-%d')}...", end=" ", flush=True)
        day_start = time.perf_counter()
        ds_surf, ds_plev = backend.fetch_day(date)
        elapsed = time.perf_counter() - day_start
        print(f"{elapsed:.1f}s")
        results.append({
            "date": date,
            "elapsed": elapsed,
            "surf_vars": list(ds_surf.data_vars),
            "plev_vars": list(ds_plev.data_vars) if ds_plev else [],
        })

    total_elapsed = time.perf_counter() - total_start
    backend.close()

    print(f"\nGoogle-only TOTAL: {total_elapsed:.1f}s ({total_elapsed/len(dates):.1f}s/day)")
    print(f"  Surface vars: {results[0]['surf_vars']}")
    print(f"  Plev vars:    {results[0]['plev_vars']}")

    return total_elapsed, results


def benchmark_hybrid(bbox, start_date, end_date, pressure_levels, cache_dir):
    """Fetch using Hybrid (S3 for surface, Google for plev+strd)."""
    print("\n" + "=" * 70)
    print("BENCHMARK: Hybrid (S3 surface + Google plev+strd)")
    print("=" * 70)

    bbox_obj = BBox.from_tuple(bbox)
    dates = pd.date_range(start_date, end_date, freq="D")

    backend = get_backend(
        "hybrid",
        bbox=bbox_obj,
        pressure_levels=pressure_levels,
        cache_dir=cache_dir,
        start_date=start_date,
        end_date=end_date,
    )

    total_start = time.perf_counter()
    results = []

    for i, date in enumerate(dates):
        print(f"  [{i+1}/{len(dates)}] Fetching {date.strftime('%Y-%m-%d')}...", end=" ", flush=True)
        day_start = time.perf_counter()
        ds_surf, ds_plev = backend.fetch_day(date)
        elapsed = time.perf_counter() - day_start
        print(f"{elapsed:.1f}s")
        results.append({
            "date": date,
            "elapsed": elapsed,
            "surf_vars": list(ds_surf.data_vars),
            "plev_vars": list(ds_plev.data_vars) if ds_plev else [],
        })

    total_elapsed = time.perf_counter() - total_start
    backend.close()

    print(f"\nHybrid TOTAL: {total_elapsed:.1f}s ({total_elapsed/len(dates):.1f}s/day)")
    print(f"  Surface vars: {results[0]['surf_vars']}")
    print(f"  Plev vars:    {results[0]['plev_vars']}")

    return total_elapsed, results


def main():
    # Test configuration
    bbox = (7.5, 46.0, 8.5, 47.0)  # Switzerland
    start_date = "2023-06-01"
    end_date = "2023-06-03"  # 3 days
    pressure_levels = [500, 700, 850]

    print("=" * 70)
    print("FAIR COMPARISON: Google-only vs Hybrid")
    print("=" * 70)
    print(f"BBox:     {bbox}")
    print(f"Dates:    {start_date} to {end_date}")
    print(f"Plev:     {pressure_levels}")
    print("=" * 70)

    with tempfile.TemporaryDirectory() as cache_dir:
        # Run Google-only benchmark
        google_time, google_results = benchmark_google_only(
            bbox, start_date, end_date, pressure_levels, cache_dir
        )

    with tempfile.TemporaryDirectory() as cache_dir:
        # Run Hybrid benchmark
        hybrid_time, hybrid_results = benchmark_hybrid(
            bbox, start_date, end_date, pressure_levels, cache_dir
        )

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Backend':<20} {'Total':<12} {'Per Day':<12} {'Speedup':<12}")
    print("-" * 70)

    n_days = len(google_results)
    print(f"{'Google-only':<20} {google_time:>8.1f}s {google_time/n_days:>9.1f}s {'(baseline)':<12}")

    if hybrid_time > 0:
        speedup = google_time / hybrid_time
        speedup_str = f"{speedup:.2f}x" if speedup >= 1 else f"{1/speedup:.2f}x slower"
        print(f"{'Hybrid (S3+Google)':<20} {hybrid_time:>8.1f}s {hybrid_time/n_days:>9.1f}s {speedup_str:<12}")

    print("-" * 70)

    # Variable comparison
    print("\nVariables retrieved:")
    print(f"  Google surface: {google_results[0]['surf_vars']}")
    print(f"  Hybrid surface: {hybrid_results[0]['surf_vars']}")
    print(f"  Google plev:    {google_results[0]['plev_vars']}")
    print(f"  Hybrid plev:    {hybrid_results[0]['plev_vars']}")


if __name__ == "__main__":
    main()
