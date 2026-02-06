#!/usr/bin/env python
"""Fair time scaling comparison: Google surface-only vs S3 surface-only.

Both fetch the same variables (surface only) for increasing time periods.
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


def benchmark_scaling(backend_name: str, bbox: tuple, days_list: list[int], cache_dir: str):
    """Benchmark time scaling for a backend (surface only)."""
    print(f"\n{'='*60}")
    print(f"TIME SCALING: {backend_name} (surface only)")
    print(f"{'='*60}")

    bbox_obj = BBox.from_tuple(bbox)
    results = []

    for n_days in days_list:
        start_date = "2023-06-01"
        end_date = pd.Timestamp(start_date) + pd.Timedelta(days=n_days-1)
        end_date_str = end_date.strftime("%Y-%m-%d")
        dates = pd.date_range(start_date, end_date_str, freq="D")

        print(f"\n  {n_days} days ({start_date} to {end_date_str}):")

        # Configure backend for surface-only
        kwargs = {
            "bbox": bbox_obj,
            "cache_dir": cache_dir,
        }

        if backend_name == "google":
            kwargs["pressure_levels"] = []  # Surface only
            # Use explicit surface vars to skip pressure level fetch
            # Must pass a non-empty list to avoid default which includes relative_humidity
            kwargs["surf_vars"] = [
                "2m_temperature",
                "2m_dewpoint_temperature",
                "surface_pressure",
                "mean_sea_level_pressure",
                "total_precipitation",
                "surface_solar_radiation_downwards",
            ]
            kwargs["plev_vars"] = ["temperature"]  # Dummy - won't be used with empty pressure_levels
        # openmeteo_s3 doesn't take pressure_levels arg

        try:
            backend = get_backend(backend_name, **kwargs)

            total_start = time.perf_counter()
            for i, date in enumerate(dates):
                day_start = time.perf_counter()
                ds_surf, ds_plev = backend.fetch_day(date)
                day_elapsed = time.perf_counter() - day_start
                print(f"    Day {i+1}: {day_elapsed:.1f}s — vars: {list(ds_surf.data_vars)}")

            total_elapsed = time.perf_counter() - total_start
            per_day = total_elapsed / n_days

            backend.close()

            results.append({
                "days": n_days,
                "total_s": total_elapsed,
                "per_day_s": per_day,
            })

            print(f"    TOTAL: {total_elapsed:.1f}s ({per_day:.1f}s/day)")

        except Exception as e:
            print(f"    ERROR: {e}")
            results.append({
                "days": n_days,
                "total_s": 0,
                "per_day_s": 0,
            })

    return results


def main():
    # Medium bbox: Western Alps
    bbox = (5.0, 44.0, 10.0, 47.0)  # ~5° x 3° = ~15 square degrees
    days_list = [1, 7, 14]  # Practical range

    print("="*60)
    print("FAIR TIME SCALING: Google vs S3 (surface only)")
    print("="*60)
    print(f"BBox: {bbox} (Western Alps, ~5° x 3°)")
    print(f"Days: {days_list}")
    print("Both backends fetch surface variables only")
    print("="*60)

    all_results = {}

    # Google surface-only
    with tempfile.TemporaryDirectory() as cache_dir:
        all_results["google"] = benchmark_scaling("google", bbox, days_list, cache_dir)

    # S3 surface-only
    with tempfile.TemporaryDirectory() as cache_dir:
        all_results["openmeteo_s3"] = benchmark_scaling("openmeteo_s3", bbox, days_list, cache_dir)

    # Summary table
    print("\n" + "="*70)
    print("SUMMARY: Time Scaling (Surface Only)")
    print("="*70)
    print(f"{'Days':<6} {'Google Total':<14} {'Google/day':<12} {'S3 Total':<14} {'S3/day':<12}")
    print("-"*70)

    for i, n_days in enumerate(days_list):
        g = all_results["google"][i]
        s = all_results["openmeteo_s3"][i]
        print(f"{n_days:<6} {g['total_s']:>10.1f}s   {g['per_day_s']:>8.1f}s/day {s['total_s']:>10.1f}s   {s['per_day_s']:>8.1f}s/day")

    # Speedup comparison
    print("\n" + "-"*70)
    print("SPEEDUP (S3 vs Google):")
    print("-"*70)
    for i, n_days in enumerate(days_list):
        g = all_results["google"][i]
        s = all_results["openmeteo_s3"][i]
        if s["per_day_s"] > 0 and g["per_day_s"] > 0:
            speedup = g["per_day_s"] / s["per_day_s"]
            print(f"  {n_days} days: S3 is {speedup:.1f}x {'faster' if speedup > 1 else 'slower'} than Google")


if __name__ == "__main__":
    main()
