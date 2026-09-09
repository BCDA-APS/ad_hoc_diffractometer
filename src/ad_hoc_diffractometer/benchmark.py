# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
benchmark.py — forward/inverse performance measurement.

Measures the throughput (operations per second) and round-trip accuracy
of :meth:`~geometry.AdHocDiffractometer.forward` and
:meth:`~geometry.AdHocDiffractometer.inverse` across all preset
geometries and their declared diffraction modes.

Usage from the command line::

    python -m ad_hoc_diffractometer.benchmark

Usage from Python::

    from ad_hoc_diffractometer.benchmark import benchmark_all, benchmark_geometry

    results = benchmark_all()
    results = benchmark_geometry("fourcv")

Each result is a dict with keys:

- ``geometry`` (str): geometry name
- ``mode`` (str): mode name
- ``status`` (str): one of ``"ok"``, ``"no_solutions"``,
  ``"not_implemented"``, ``"error"``
- ``forward_ops_per_sec`` (float or None)
- ``inverse_ops_per_sec`` (float or None)
- ``forward_inverse_ratio`` (float or None): dimensionless ratio
  ``forward_ops/sec / inverse_ops/sec``.  This workstation-independent
  metric measures how fast ``forward()`` is relative to ``inverse()``
  on the same machine.  **Higher is better.**  A value of 1.0 means
  they are equally fast; lower values indicate ``forward()`` is slower.
  Since ``inverse()`` is a single direct computation (no iteration),
  its speed characterizes the machine; the ratio captures pure
  algorithmic efficiency of the forward solver.
- ``round_trip_max_error`` (float or None)
- ``n_reflections`` (int): number of reflections attempted
- ``n_solutions`` (int): total solutions returned
- ``error_message`` (str or None): exception text when status is ``"error"``
"""

from __future__ import annotations

import logging
import time

from .factories import list_geometries
from .factories import make_geometry
from .lattice import Lattice
from .mode import REQUIRED
from .orientation import ub_identity

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default benchmark parameters
# ---------------------------------------------------------------------------

DEFAULT_WAVELENGTH: float = 1.5406
"""Cu Kα wavelength in Å."""

DEFAULT_LATTICE_A: float = 4.0
"""Cubic lattice parameter in Å."""

DEFAULT_REFLECTIONS: list[tuple[float, float, float]] = [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 0),
    (1, 1, 1),
]
"""Reflections to benchmark.  All reachable with a=4 Å at Cu Kα."""

DEFAULT_N_ITER: int = 100
"""Number of timing iterations per operation."""


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _setup_geometry(name: str):
    """Create a geometry with a standard cubic sample and UB = B."""
    g = make_geometry(name)
    g.wavelength = DEFAULT_WAVELENGTH
    g.sample.lattice = Lattice(a=DEFAULT_LATTICE_A)
    ub_identity(g.sample)
    return g


def _prepare_mode(geometry, mode_name: str) -> None:
    """Activate a mode and configure any special prerequisites.

    Sets ``mode_name`` on the geometry and provides the minimum setup
    required for modes that need extra configuration:

    - fixed_psi modes: sets ``azimuth``
    - double_diffraction modes: sets h2/k2/l2 extras
    - zone modes: sets z0/z1 extras to a generic (h,k,0) plane
    - surface/reference modes (incidence, emergence, incidence_equals_emergence):
      sets ``surface_normal``
    """
    geometry.mode_name = mode_name
    cs = geometry.modes[mode_name]

    # Double-diffraction modes: replace REQUIRED sentinels with values
    for key in ("h2", "k2", "l2"):
        if key in cs.extras and cs.extras[key] is REQUIRED:
            defaults = {"h2": 0.0, "k2": 1.0, "l2": 0.0}
            cs.extras[key] = defaults[key]

    # Zone modes: replace REQUIRED z0/z1 sentinels with a generic plane
    for key, default in (("z0", (1, 0, 0)), ("z1", (0, 1, 0))):
        if key in cs.extras and cs.extras[key] is REQUIRED:
            cs.extras[key] = default

    # Reference-vector modes: set azimuth if needed
    for c in cs._constraints:
        cname = getattr(c, "_name", getattr(c, "name", ""))
        if cname == "psi" and geometry.azimuth is None:
            geometry.azimuth = (0, 0, 1)

    # Surface modes: set surface_normal if needed.
    for c in cs._constraints:
        cname = getattr(c, "_name", getattr(c, "name", ""))
        if cname in ("incidence", "emergence", "incidence_equals_emergence"):
            if geometry.surface_normal is None:
                geometry.surface_normal = (0, 0, 1)


def _time_forward(geometry, reflections, n_iter: int) -> tuple[float, list]:
    """Time forward() calls and collect solutions.

    Returns (ops_per_sec, all_solutions) where all_solutions is a list
    of (hkl, solutions) pairs from the final iteration.
    """
    all_solutions = []

    start = time.perf_counter()
    for _ in range(n_iter):
        all_solutions.clear()
        for hkl in reflections:
            solutions = geometry.forward(*hkl)
            all_solutions.append((hkl, solutions))
    elapsed = time.perf_counter() - start

    total_ops = n_iter * len(reflections)
    ops_per_sec = total_ops / elapsed if elapsed > 0 else float("inf")
    return ops_per_sec, all_solutions


def _time_inverse(geometry, solutions_by_hkl, n_iter: int) -> tuple[float, float]:
    """Time inverse() calls and measure round-trip error.

    Parameters
    ----------
    solutions_by_hkl : list of (hkl, solutions)
        Output from _time_forward.

    Returns (ops_per_sec, max_error).
    """
    # Flatten to (hkl, angle_dict) pairs
    pairs = []
    for hkl, solutions in solutions_by_hkl:
        for sol in solutions:
            pairs.append((hkl, sol))

    if not pairs:
        return 0.0, 0.0

    max_error = 0.0

    start = time.perf_counter()
    for _ in range(n_iter):
        for hkl, angles in pairs:
            hkl_back = geometry.inverse(angles)
            err = max(abs(a - b) for a, b in zip(hkl, hkl_back, strict=False))
            if err > max_error:
                max_error = err
    elapsed = time.perf_counter() - start

    total_ops = n_iter * len(pairs)
    ops_per_sec = total_ops / elapsed if elapsed > 0 else float("inf")
    return ops_per_sec, max_error


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def benchmark_mode(
    geometry_name: str,
    mode_name: str,
    *,
    reflections: list[tuple[float, float, float]] | None = None,
    n_iter: int = DEFAULT_N_ITER,
) -> dict:
    """Benchmark a single geometry + mode combination.

    Parameters
    ----------
    geometry_name : str
        Name of a registered geometry (e.g. ``"fourcv"``).
    mode_name : str
        Name of a mode declared on the geometry.
    reflections : list of (h, k, l) or None
        Reflections to use.  Defaults to :data:`DEFAULT_REFLECTIONS`.
    n_iter : int
        Number of timing iterations.  Defaults to :data:`DEFAULT_N_ITER`.

    Returns
    -------
    dict
        Result dict with keys: ``geometry``, ``mode``, ``status``,
        ``forward_ops_per_sec``, ``inverse_ops_per_sec``,
        ``round_trip_max_error``, ``n_reflections``, ``n_solutions``,
        ``error_message``.
    """
    if reflections is None:
        reflections = DEFAULT_REFLECTIONS

    result = {
        "geometry": geometry_name,
        "mode": mode_name,
        "status": "ok",
        "forward_ops_per_sec": None,
        "inverse_ops_per_sec": None,
        "forward_inverse_ratio": None,
        "round_trip_max_error": None,
        "n_reflections": len(reflections),
        "n_solutions": 0,
        "error_message": None,
    }

    try:
        g = _setup_geometry(geometry_name)
        _prepare_mode(g, mode_name)

        fwd_ops, solutions_by_hkl = _time_forward(g, reflections, n_iter)
        n_solutions = sum(len(sols) for _, sols in solutions_by_hkl)

        result["forward_ops_per_sec"] = fwd_ops
        result["n_solutions"] = n_solutions

        if n_solutions == 0:
            result["status"] = "no_solutions"
            return result

        inv_ops, max_err = _time_inverse(g, solutions_by_hkl, n_iter)
        result["inverse_ops_per_sec"] = inv_ops
        result["round_trip_max_error"] = max_err

        # Workstation-independent ratio: forward speed as a fraction of
        # inverse speed.  Higher is better.  Since inverse() is a single
        # direct computation (no iteration), its speed characterizes the
        # machine; the ratio captures pure algorithmic efficiency of the
        # forward solver.
        if fwd_ops > 0 and inv_ops > 0:  # pragma: no branch
            result["forward_inverse_ratio"] = fwd_ops / inv_ops

    except NotImplementedError as exc:
        result["status"] = "not_implemented"
        result["error_message"] = str(exc)
    except Exception as exc:  # noqa: BLE001
        result["status"] = "error"
        result["error_message"] = str(exc)

    return result


def benchmark_geometry(
    name: str,
    *,
    reflections: list[tuple[float, float, float]] | None = None,
    n_iter: int = DEFAULT_N_ITER,
    verbose: bool = True,
) -> list[dict]:
    """Benchmark all modes of a single geometry.

    Parameters
    ----------
    name : str
        Name of a registered geometry (e.g. ``"fourcv"``).
    reflections : list of (h, k, l) or None
        Reflections to use.  Defaults to :data:`DEFAULT_REFLECTIONS`.
    n_iter : int
        Number of timing iterations.
    verbose : bool
        If True, print a results table to stdout.

    Returns
    -------
    list of dict
        One result dict per mode.
    """
    g = make_geometry(name)
    mode_names = sorted(g.modes.keys())

    results = []
    for mode_name in mode_names:
        r = benchmark_mode(name, mode_name, reflections=reflections, n_iter=n_iter)
        results.append(r)

    if verbose:
        _print_results(results)

    return results


def benchmark_all(
    *,
    reflections: list[tuple[float, float, float]] | None = None,
    n_iter: int = DEFAULT_N_ITER,
    verbose: bool = True,
) -> list[dict]:
    """Benchmark all modes of all registered geometries.

    Parameters
    ----------
    reflections : list of (h, k, l) or None
        Reflections to use.  Defaults to :data:`DEFAULT_REFLECTIONS`.
    n_iter : int
        Number of timing iterations.
    verbose : bool
        If True, print a results table to stdout.

    Returns
    -------
    list of dict
        One result dict per geometry+mode combination.
    """
    all_results = []
    for name in sorted(list_geometries()):
        results = benchmark_geometry(
            name, reflections=reflections, n_iter=n_iter, verbose=False
        )
        all_results.extend(results)

    if verbose:
        _print_results(all_results)

    return all_results


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

_STATUS_WIDTH = 15
_OPS_WIDTH = 10


def _print_results(results: list[dict]) -> None:
    """Print a formatted table of benchmark results to stdout."""
    import builtins

    header = (
        f"{'geometry':<12s}  {'mode':<32s}  {'status':<{_STATUS_WIDTH}s}"
        f"  {'fwd ops/s':>{_OPS_WIDTH}s}  {'inv ops/s':>{_OPS_WIDTH}s}"
        f"  {'fwd/inv':>7s}"
        f"  {'round-trip err':>14s}  {'solns':>5s}"
    )
    separator = "-" * len(header)

    builtins.print(separator)
    builtins.print(header)
    builtins.print(separator)

    for r in results:
        fwd = (
            f"{r['forward_ops_per_sec']:>{_OPS_WIDTH},.0f}"
            if r["forward_ops_per_sec"] is not None
            else f"{'-':>{_OPS_WIDTH}s}"
        )
        inv = (
            f"{r['inverse_ops_per_sec']:>{_OPS_WIDTH},.0f}"
            if r["inverse_ops_per_sec"] is not None
            else f"{'-':>{_OPS_WIDTH}s}"
        )
        ratio = (
            f"{r['forward_inverse_ratio']:>7.4f}"
            if r.get("forward_inverse_ratio") is not None
            else f"{'-':>7s}"
        )
        err = (
            f"{r['round_trip_max_error']:>14.2e}"
            if r["round_trip_max_error"] is not None
            else f"{'-':>14s}"
        )
        builtins.print(
            f"{r['geometry']:<12s}  {r['mode']:<32s}  {r['status']:<{_STATUS_WIDTH}s}"
            f"  {fwd}  {inv}  {ratio}  {err}  {r['n_solutions']:>5d}"
        )

    builtins.print(separator)


# ---------------------------------------------------------------------------
# CLI entry point: python -m ad_hoc_diffractometer.benchmark
# ---------------------------------------------------------------------------

if __name__ == "__main__":  # pragma: no cover
    benchmark_all()
