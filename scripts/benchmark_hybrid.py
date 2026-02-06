#!/usr/bin/env python
"""Benchmark script for hybrid backend vs single backends.

Usage:
    python scripts/benchmark_hybrid.py [--scenario N] [--all]

Scenarios:
    1. Regional (s3zarr) - Central Asia, pre-2023
    2. Hybrid surface+plev - Switzerland, 2023
    3. Hybrid pre-2022 - Switzerland, 2019 (no IFS precip)
    4. Large region - Scandinavia surface only
    5. Backend comparison - Same bbox, different backends
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ecmwf_downloader.bbox import BBox


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""
    scenario: str
    backend: str
    bbox: str
    dates: str
    plev: str
    total_seconds: float
    per_day_seconds: float
    n_days: int
    surface_vars: list[str]
    plev_vars: list[str]
    error: str | None = None


@contextmanager
def timer():
    """Context manager for timing."""
    start = time.perf_counter()
    result = {"elapsed": 0.0}
    try:
        yield result
    finally:
        result["elapsed"] = time.perf_counter() - start


def run_benchmark(
    name: str,
    backend_name: str,
    bbox: tuple[float, float, float, float],
    start_date: str,
    end_date: str,
    pressure_levels: list[int] | None = None,
    cache_dir: str | None = None,
) -> BenchmarkResult:
    """Run a single benchmark."""
    from ecmwf_downloader.reanalysis.backends import get_backend, _check_available

    # Check if backend is available
    if backend_name != "hybrid" and not _check_available(backend_name):
        print(f"\n  Skipping {backend_name}: dependencies not installed")
        return BenchmarkResult(
            scenario=name,
            backend=backend_name,
            bbox=f"{bbox[0]:.1f},{bbox[1]:.1f},{bbox[2]:.1f},{bbox[3]:.1f}",
            dates=f"{start_date} to {end_date}",
            plev=str(pressure_levels) if pressure_levels else "None",
            total_seconds=0,
            per_day_seconds=0,
            n_days=0,
            surface_vars=[],
            plev_vars=[],
            error="Dependencies not installed",
        )

    bbox_obj = BBox.from_tuple(bbox)
    start = pd.Timestamp(start_date)
    end = pd.Timestamp(end_date)
    dates = pd.date_range(start, end, freq="D")
    n_days = len(dates)

    cache = cache_dir or tempfile.mkdtemp(prefix=f"bench_{backend_name}_")

    print(f"\n{'='*60}")
    print(f"Scenario: {name}")
    print(f"Backend:  {backend_name}")
    print(f"BBox:     {bbox}")
    print(f"Dates:    {start_date} to {end_date} ({n_days} days)")
    print(f"Plev:     {pressure_levels or 'None'}")
    print(f"{'='*60}")

    surface_vars = []
    plev_vars = []
    error = None

    try:
        # Initialize backend
        backend_kwargs = {
            "bbox": bbox_obj,
            "cache_dir": cache,
        }
        # Most backends require pressure_levels (even if empty list)
        if backend_name in ["google", "cds", "s3zarr"]:
            backend_kwargs["pressure_levels"] = pressure_levels or []
            # Use variables that exist in Google ARCO (no relative_humidity)
            backend_kwargs["plev_vars"] = ["temperature", "geopotential", "u_component_of_wind",
                                           "v_component_of_wind", "specific_humidity"]
        elif backend_name == "openmeteo":
            backend_kwargs["pressure_levels"] = pressure_levels or []
            backend_kwargs["start_date"] = start_date
            backend_kwargs["end_date"] = end_date
        elif backend_name == "hybrid":
            if pressure_levels:
                backend_kwargs["pressure_levels"] = pressure_levels
            backend_kwargs["start_date"] = start_date
            backend_kwargs["end_date"] = end_date
        # openmeteo_s3 doesn't take pressure_levels

        backend = get_backend(backend_name, **backend_kwargs)

        # Fetch all days
        with timer() as t:
            for i, date in enumerate(dates):
                print(f"  Fetching {date.strftime('%Y-%m-%d')} ({i+1}/{n_days})...", end=" ")
                day_start = time.perf_counter()
                ds_surf, ds_plev = backend.fetch_day(date)
                day_elapsed = time.perf_counter() - day_start
                print(f"{day_elapsed:.1f}s")

                if i == 0:
                    surface_vars = list(ds_surf.data_vars)
                    plev_vars = list(ds_plev.data_vars) if ds_plev else []

        backend.close()

    except Exception as e:
        error = str(e)
        print(f"  ERROR: {error}")
        t = {"elapsed": 0.0}

    total = t["elapsed"]
    per_day = total / n_days if n_days > 0 else 0

    result = BenchmarkResult(
        scenario=name,
        backend=backend_name,
        bbox=f"{bbox[0]:.1f},{bbox[1]:.1f},{bbox[2]:.1f},{bbox[3]:.1f}",
        dates=f"{start_date} to {end_date}",
        plev=str(pressure_levels) if pressure_levels else "None",
        total_seconds=total,
        per_day_seconds=per_day,
        n_days=n_days,
        surface_vars=surface_vars,
        plev_vars=plev_vars,
        error=error,
    )

    print(f"\nResult: {total:.1f}s total, {per_day:.1f}s/day")
    print(f"Surface vars: {surface_vars}")
    print(f"Plev vars:    {plev_vars}")

    return result


def scenario_1_regional() -> list[BenchmarkResult]:
    """Scenario 1: Regional (s3zarr) - Central Asia, pre-2023."""
    results = []

    # Central Asia bbox (within s3zarr coverage)
    bbox = (68.0, 39.0, 72.0, 42.0)  # Tajikistan
    start = "2020-01-01"
    end = "2020-01-03"
    plev = [500, 700, 850]

    # Test hybrid (should auto-select s3zarr)
    results.append(run_benchmark(
        "1. Regional (hybrid→s3zarr)",
        "hybrid",
        bbox, start, end, plev,
    ))

    # Test s3zarr directly
    results.append(run_benchmark(
        "1. Regional (s3zarr direct)",
        "s3zarr",
        bbox, start, end, plev,
    ))

    return results


def scenario_2_hybrid_plev() -> list[BenchmarkResult]:
    """Scenario 2: Hybrid surface+plev - Switzerland, 2023."""
    results = []

    # Switzerland (outside s3zarr coverage)
    bbox = (7.5, 46.0, 8.5, 47.0)
    start = "2023-06-01"
    end = "2023-06-03"
    plev = [500, 700, 850]

    # Test hybrid
    results.append(run_benchmark(
        "2. Hybrid plev (hybrid)",
        "hybrid",
        bbox, start, end, plev,
    ))

    # Test google for comparison
    results.append(run_benchmark(
        "2. Hybrid plev (google)",
        "google",
        bbox, start, end, plev,
    ))

    return results


def scenario_3_hybrid_pre2022() -> list[BenchmarkResult]:
    """Scenario 3: Hybrid pre-2022 - Switzerland, 2019 (no IFS precip)."""
    results = []

    bbox = (7.5, 46.0, 8.5, 47.0)
    start = "2019-01-01"
    end = "2019-01-03"
    plev = [850]

    # Test hybrid (should use S3 for surface, Google for plev)
    results.append(run_benchmark(
        "3. Pre-2022 (hybrid)",
        "hybrid",
        bbox, start, end, plev,
    ))

    # Test google for comparison
    results.append(run_benchmark(
        "3. Pre-2022 (google)",
        "google",
        bbox, start, end, plev,
    ))

    return results


def scenario_4_large_region() -> list[BenchmarkResult]:
    """Scenario 4: Large region - Scandinavia surface only."""
    results = []

    # Scandinavia (~6000 grid points)
    bbox = (5.0, 55.0, 30.0, 71.0)
    start = "2023-01-01"
    end = "2023-01-02"  # Just 2 days for large region
    plev = None

    # Test hybrid (should use S3 for surface)
    results.append(run_benchmark(
        "4. Large region (hybrid)",
        "hybrid",
        bbox, start, end, plev,
    ))

    # Note: API would likely hit rate limits, skip by default
    # results.append(run_benchmark(
    #     "4. Large region (openmeteo API)",
    #     "openmeteo",
    #     bbox, start, end, plev,
    # ))

    return results


def scenario_5_backend_comparison() -> list[BenchmarkResult]:
    """Scenario 5: Backend comparison - Same bbox, different backends."""
    results = []

    # Small region for fair comparison
    bbox = (7.7, 45.9, 8.0, 46.1)  # Tiny Swiss bbox
    start = "2023-03-01"
    end = "2023-03-03"
    plev = [850]  # Include one pressure level

    backends = ["openmeteo", "google"]  # Skip hybrid for now (needs S3)

    for backend in backends:
        try:
            results.append(run_benchmark(
                f"5. Comparison ({backend})",
                backend,
                bbox, start, end, plev,
            ))
        except Exception as e:
            print(f"  Skipping {backend}: {e}")

    return results


def print_summary(results: list[BenchmarkResult]):
    """Print summary table of results."""
    print("\n" + "=" * 80)
    print("BENCHMARK SUMMARY")
    print("=" * 80)

    # Group by scenario
    scenarios = {}
    for r in results:
        if r.scenario not in scenarios:
            scenarios[r.scenario] = []
        scenarios[r.scenario].append(r)

    print(f"\n{'Scenario':<40} {'Backend':<15} {'Total':<10} {'Per Day':<10} {'Status'}")
    print("-" * 80)

    for scenario, runs in scenarios.items():
        for r in runs:
            status = "OK" if not r.error else f"ERR: {r.error[:20]}"
            print(f"{r.scenario:<40} {r.backend:<15} {r.total_seconds:>7.1f}s {r.per_day_seconds:>7.1f}s {status}")

    # Print speedups
    print("\n" + "-" * 80)
    print("SPEEDUPS (vs Google baseline)")
    print("-" * 80)

    for scenario, runs in scenarios.items():
        google_time = next((r.total_seconds for r in runs if "google" in r.backend.lower()), None)
        if google_time and google_time > 0:
            for r in runs:
                if r.backend != "google" and r.total_seconds > 0:
                    speedup = google_time / r.total_seconds
                    print(f"{r.scenario}: {r.backend} is {speedup:.1f}x {'faster' if speedup > 1 else 'slower'}")


def main():
    parser = argparse.ArgumentParser(description="Benchmark hybrid backend")
    parser.add_argument("--scenario", type=int, choices=[1, 2, 3, 4, 5],
                        help="Run specific scenario (1-5)")
    parser.add_argument("--all", action="store_true",
                        help="Run all scenarios")
    parser.add_argument("--quick", action="store_true",
                        help="Run quick comparison only (scenario 5)")
    args = parser.parse_args()

    results = []

    if args.quick or (not args.scenario and not args.all):
        print("Running quick comparison (scenario 5)...")
        results.extend(scenario_5_backend_comparison())

    elif args.all:
        print("Running all scenarios...")
        results.extend(scenario_1_regional())
        results.extend(scenario_2_hybrid_plev())
        results.extend(scenario_3_hybrid_pre2022())
        results.extend(scenario_4_large_region())
        results.extend(scenario_5_backend_comparison())

    elif args.scenario == 1:
        results.extend(scenario_1_regional())
    elif args.scenario == 2:
        results.extend(scenario_2_hybrid_plev())
    elif args.scenario == 3:
        results.extend(scenario_3_hybrid_pre2022())
    elif args.scenario == 4:
        results.extend(scenario_4_large_region())
    elif args.scenario == 5:
        results.extend(scenario_5_backend_comparison())

    print_summary(results)

    return results


if __name__ == "__main__":
    main()
