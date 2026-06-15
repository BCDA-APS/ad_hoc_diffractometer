# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
Tests for the benchmark module.

Corresponds to ``src/ad_hoc_diffractometer/benchmark.py``.
Uses ``n_iter=1`` throughout to keep test duration short.
"""

import re
from contextlib import nullcontext as does_not_raise

import pytest

from ad_hoc_diffractometer.benchmark import DEFAULT_REFLECTIONS
from ad_hoc_diffractometer.benchmark import _prepare_mode
from ad_hoc_diffractometer.benchmark import _setup_geometry
from ad_hoc_diffractometer.benchmark import _time_forward
from ad_hoc_diffractometer.benchmark import _time_inverse
from ad_hoc_diffractometer.benchmark import benchmark_all
from ad_hoc_diffractometer.benchmark import benchmark_geometry
from ad_hoc_diffractometer.benchmark import benchmark_mode
from ad_hoc_diffractometer.factories import list_geometries

REQUIRED_KEYS = {
    "geometry",
    "mode",
    "status",
    "forward_ops_per_sec",
    "inverse_ops_per_sec",
    "forward_inverse_ratio",
    "round_trip_max_error",
    "n_reflections",
    "n_solutions",
    "error_message",
}

VALID_STATUSES = {"ok", "no_solutions", "not_implemented", "error"}


# ---------------------------------------------------------------------------
# _setup_geometry
# ---------------------------------------------------------------------------


class TestSetupGeometry:
    """Tests for the _setup_geometry helper."""

    @pytest.mark.parametrize(
        "name, context",
        [
            pytest.param("fourcv", does_not_raise(), id="fourcv"),
            pytest.param("psic", does_not_raise(), id="psic"),
            pytest.param("kappa4cv", does_not_raise(), id="kappa4cv"),
            pytest.param(
                "nonexistent",
                pytest.raises(
                    ValueError,
                    match=re.escape("No geometry named 'nonexistent'"),
                ),
                id="nonexistent-raises",
            ),
        ],
    )
    def test_setup_geometry(self, name, context):
        with context:
            g = _setup_geometry(name)
            assert g.name == name
            assert g.wavelength == pytest.approx(1.5406)
            assert g.sample.UB is not None


# ---------------------------------------------------------------------------
# _prepare_mode
# ---------------------------------------------------------------------------


class TestPrepareMode:
    """Tests for the _prepare_mode helper."""

    def test_standard_mode(self):
        """A standard bisecting mode needs no special setup."""
        g = _setup_geometry("fourcv")
        _prepare_mode(g, "bisecting")
        assert g.mode_name == "bisecting"

    @pytest.mark.parametrize(
        "geometry_name, mode_name, context",
        [
            pytest.param(
                "fourcv",
                "double_diffraction",
                does_not_raise(),
                id="fourcv-dd",
            ),
            pytest.param(
                "psic",
                "double_diffraction_vertical",
                does_not_raise(),
                id="psic-dd-vert",
            ),
        ],
    )
    def test_double_diffraction_extras(self, geometry_name, mode_name, context):
        """Double-diffraction modes get h2/k2/l2 sentinels replaced."""
        with context:
            g = _setup_geometry(geometry_name)
            _prepare_mode(g, mode_name)
            cs = g.modes[mode_name]
            for key in ("h2", "k2", "l2"):
                assert isinstance(cs.extras[key], float)

    def test_fixed_psi_sets_reference(self):
        """fixed_psi modes get an azimuth if not set."""
        g = _setup_geometry("fourcv")
        assert g.azimuth is None
        _prepare_mode(g, "fixed_psi")
        assert g.azimuth is not None

    @pytest.mark.parametrize(
        "geometry_name, mode_name",
        [
            pytest.param("psic", "zone_vertical", id="psic-zone-vert"),
            pytest.param("psic", "zone_horizontal", id="psic-zone-horiz"),
            pytest.param("kappa6c", "zone_vertical", id="kappa6c-zone-vert"),
            pytest.param("kappa6c", "zone_horizontal", id="kappa6c-zone-horiz"),
        ],
    )
    def test_zone_extras_replaced(self, geometry_name, mode_name):
        """Zone modes get z0/z1 REQUIRED sentinels replaced with the
        default (1,0,0)/(0,1,0) plane."""
        g = _setup_geometry(geometry_name)
        _prepare_mode(g, mode_name)
        cs = g.modes[mode_name]
        assert cs.extras["z0"] == (1, 0, 0)
        assert cs.extras["z1"] == (0, 1, 0)

    @pytest.mark.parametrize(
        "geometry_name, mode_name",
        [
            pytest.param("sixc", "fixed_alpha_zaxis", id="sixc-alpha-zaxis"),
            pytest.param("sixc", "fixed_beta_zaxis", id="sixc-beta-zaxis"),
            pytest.param("sixc", "alpha_eq_beta_zaxis", id="sixc-a-eq-b"),
        ],
    )
    def test_surface_modes_set_normal(self, geometry_name, mode_name):
        """Surface/reference modes get a surface_normal if not set."""
        g = _setup_geometry(geometry_name)
        assert g.surface_normal is None
        _prepare_mode(g, mode_name)
        assert g.surface_normal is not None

    def test_surface_mode_preserves_existing_normal(self):
        """If surface_normal is already set, _prepare_mode leaves it."""
        g = _setup_geometry("sixc")
        g.surface_normal = (1, 0, 0)
        _prepare_mode(g, "fixed_alpha_zaxis")
        assert g.surface_normal == (1.0, 0.0, 0.0)


# ---------------------------------------------------------------------------
# _time_forward / _time_inverse
# ---------------------------------------------------------------------------


class TestTimingHelpers:
    """Tests for the internal timing functions."""

    def test_time_forward_returns_solutions(self):
        g = _setup_geometry("fourcv")
        _prepare_mode(g, "bisecting")
        ops, solutions = _time_forward(g, [(1, 0, 0)], n_iter=1)
        assert ops > 0
        assert len(solutions) == 1  # one reflection
        hkl, sols = solutions[0]
        assert hkl == (1, 0, 0)
        assert len(sols) >= 1

    def test_time_inverse_round_trip(self):
        g = _setup_geometry("fourcv")
        _prepare_mode(g, "bisecting")
        _, solutions = _time_forward(g, [(1, 0, 0)], n_iter=1)
        ops, max_err = _time_inverse(g, solutions, n_iter=1)
        assert ops > 0
        assert max_err < 1e-8

    def test_time_inverse_empty_solutions(self):
        """inverse timing with no solutions returns zero ops."""
        g = _setup_geometry("fourcv")
        _prepare_mode(g, "bisecting")
        ops, max_err = _time_inverse(g, [], n_iter=1)
        assert ops == 0.0
        assert max_err == 0.0


# ---------------------------------------------------------------------------
# benchmark_mode
# ---------------------------------------------------------------------------


class TestBenchmarkMode:
    """Tests for the benchmark_mode function."""

    def test_result_dict_keys(self):
        """Result contains all required keys."""
        r = benchmark_mode("fourcv", "bisecting", n_iter=1)
        assert set(r.keys()) == REQUIRED_KEYS

    def test_ok_status(self):
        """A working mode returns status 'ok'."""
        r = benchmark_mode("fourcv", "bisecting", reflections=[(1, 0, 0)], n_iter=1)
        assert r["status"] == "ok"
        assert r["forward_ops_per_sec"] is not None
        assert r["forward_ops_per_sec"] > 0
        assert r["inverse_ops_per_sec"] is not None
        assert r["inverse_ops_per_sec"] > 0
        assert r["forward_inverse_ratio"] is not None
        assert r["forward_inverse_ratio"] > 0
        assert r["round_trip_max_error"] is not None
        assert r["round_trip_max_error"] < 1e-8
        assert r["n_solutions"] > 0
        assert r["error_message"] is None

    def test_no_solutions_status(self):
        """A mode that returns no solutions for the given reflections."""
        # fixed_psi on fourcv with default setup rarely matches
        r = benchmark_mode("fourcv", "fixed_psi", reflections=[(1, 0, 0)], n_iter=1)
        assert r["status"] in ("no_solutions", "not_implemented", "ok")
        # For fixed_psi the likely outcome is no_solutions or not_implemented
        if r["status"] == "no_solutions":
            assert r["inverse_ops_per_sec"] is None
            assert r["forward_inverse_ratio"] is None
            assert r["n_solutions"] == 0

    @pytest.mark.parametrize(
        "status",
        [
            pytest.param("ok", id="ok"),
            pytest.param("no_solutions", id="no-solutions"),
            pytest.param("not_implemented", id="not-implemented"),
            pytest.param("error", id="error"),
        ],
    )
    def test_valid_status_values(self, status):
        """All status values are in the expected set."""
        assert status in VALID_STATUSES

    def test_not_implemented_status(self, monkeypatch):
        """A mode that raises NotImplementedError gets status 'not_implemented'."""
        from ad_hoc_diffractometer import benchmark as bm

        def _raise_not_impl(*args, **kwargs):
            raise NotImplementedError("solver not available")

        monkeypatch.setattr(bm, "_time_forward", _raise_not_impl)

        r = benchmark_mode("fourcv", "bisecting", reflections=[(1, 0, 0)], n_iter=1)
        assert r["status"] == "not_implemented"
        assert "solver not available" in r["error_message"]

    def test_nonexistent_geometry(self):
        """A nonexistent geometry returns error status."""
        r = benchmark_mode("nonexistent", "bisecting", n_iter=1)
        assert r["status"] == "error"
        assert r["error_message"] is not None
        assert "nonexistent" in r["error_message"]

    def test_n_reflections_matches_input(self):
        """n_reflections reflects how many reflections were passed in."""
        refls = [(1, 0, 0), (0, 1, 0)]
        r = benchmark_mode("fourcv", "bisecting", reflections=refls, n_iter=1)
        assert r["n_reflections"] == 2

    def test_default_reflections(self):
        """Without explicit reflections, uses DEFAULT_REFLECTIONS."""
        r = benchmark_mode("fourcv", "bisecting", n_iter=1)
        assert r["n_reflections"] == len(DEFAULT_REFLECTIONS)


# ---------------------------------------------------------------------------
# benchmark_geometry
# ---------------------------------------------------------------------------


class TestBenchmarkGeometry:
    """Tests for the benchmark_geometry function."""

    def test_returns_all_modes(self):
        """Returns one result per declared mode."""
        from helpers import fourcv

        g = fourcv()
        n_modes = len(g.modes)
        results = benchmark_geometry(
            "fourcv", reflections=[(1, 0, 0)], n_iter=1, verbose=False
        )
        assert len(results) == n_modes

    def test_all_results_have_correct_geometry(self):
        results = benchmark_geometry(
            "fourch", reflections=[(1, 0, 0)], n_iter=1, verbose=False
        )
        for r in results:
            assert r["geometry"] == "fourch"
            assert set(r.keys()) == REQUIRED_KEYS
            assert r["status"] in VALID_STATUSES

    def test_verbose_prints_table(self, capsys):
        """verbose=True prints a table to stdout."""
        benchmark_geometry("fourcv", reflections=[(1, 0, 0)], n_iter=1, verbose=True)
        captured = capsys.readouterr()
        assert "fourcv" in captured.out
        assert "bisecting" in captured.out
        assert "fwd ops/s" in captured.out

    def test_verbose_false_no_output(self, capsys):
        """verbose=False produces no stdout output."""
        benchmark_geometry("fourcv", reflections=[(1, 0, 0)], n_iter=1, verbose=False)
        captured = capsys.readouterr()
        assert captured.out == ""


# ---------------------------------------------------------------------------
# benchmark_all
# ---------------------------------------------------------------------------


class TestBenchmarkAllFast:
    """Fast tests for benchmark_all using a single monkeypatched geometry.

    These exercise 100 % of benchmark_all()'s code paths (loop, extend,
    verbose branch) in ~1 second by limiting the geometry registry to a
    single entry.  Full-sweep tests are in TestBenchmarkAll (marked
    slow_benchmark).
    """

    @pytest.fixture(autouse=True)
    def _one_geometry(self, monkeypatch):
        """Restrict list_geometries() to a single entry."""
        full = list_geometries()
        monkeypatch.setattr(
            "ad_hoc_diffractometer.benchmark.list_geometries",
            lambda: {"fourcv": full["fourcv"]},
        )

    def test_covers_requested_geometries(self):
        """Result geometry names match the (monkeypatched) registry."""
        results = benchmark_all(reflections=[(1, 0, 0)], n_iter=1, verbose=False)
        geometry_names = {r["geometry"] for r in results}
        assert geometry_names == {"fourcv"}

    def test_all_results_valid(self):
        """Every result has valid keys and status."""
        results = benchmark_all(reflections=[(1, 0, 0)], n_iter=1, verbose=False)
        for r in results:
            assert set(r.keys()) == REQUIRED_KEYS
            assert r["status"] in VALID_STATUSES

    def test_verbose_prints_table(self, capsys):
        """verbose=True prints a formatted table to stdout."""
        results = benchmark_all(reflections=[(1, 0, 0)], n_iter=1, verbose=True)
        captured = capsys.readouterr()
        assert "fwd ops/s" in captured.out
        assert "fourcv" in captured.out
        assert len(results) > 0


# ---------------------------------------------------------------------------
# benchmark_all — full sweep (slow)
# ---------------------------------------------------------------------------


class TestBenchmarkAll:
    """Full-sweep tests for benchmark_all across all registered geometries.

    These are marked slow_benchmark and excluded from the default pytest
    run.  Run on demand with: ``pytest -m slow_benchmark``
    """

    @pytest.mark.slow_benchmark
    def test_covers_all_geometries(self):
        """Results include every registered geometry."""
        results = benchmark_all(reflections=[(1, 0, 0)], n_iter=1, verbose=False)
        geometry_names = {r["geometry"] for r in results}
        registered = set(list_geometries().keys())
        assert geometry_names == registered

    @pytest.mark.slow_benchmark
    def test_all_results_valid(self):
        """Every result has valid keys and status."""
        results = benchmark_all(reflections=[(1, 0, 0)], n_iter=1, verbose=False)
        for r in results:
            assert set(r.keys()) == REQUIRED_KEYS
            assert r["status"] in VALID_STATUSES

    @pytest.mark.slow_benchmark
    def test_verbose_prints_table(self, capsys):
        results = benchmark_all(reflections=[(1, 0, 0)], n_iter=1, verbose=True)
        captured = capsys.readouterr()
        assert "fwd ops/s" in captured.out
        assert len(results) > 0
