# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
Regression tests for issue #307.

Issue #307: for geometry ``psic``, mode ``fixed_incidence_vertical``,
``forward(h, k, l)`` returned solutions whose ``chi`` and ``phi`` were
pinned at ``0`` and whose ``inverse()`` did **not** reproduce the
requested reflection.  The returned motor angles satisfied the incidence
constraint but were not a valid forward solution for the given hkl.

Root cause: ``fixed_incidence_vertical`` leaves three sample stages free
(``eta, chi, phi``), but ``_solve_surface`` did a 1-D Newton over only the
first free sample stage (``eta``), leaving ``chi`` and ``phi`` at their
current values and enforcing only the incidence residual — never the
Bragg/Q condition.

Fix: ``_solve_surface`` routes the multi-free-sample case to
``_solve_reference_three_sample``, which scans the outer sample stage,
solves the inner ``(chi, phi)`` pair against the Bragg condition for each
value, and roots the scan on the surface residual.  Every returned
solution reproduces the requested hkl by construction.
"""

from __future__ import annotations

import warnings
from contextlib import nullcontext as does_not_raise

import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.reference import incidence_angle


def _silicon_psic():
    """psic geometry with a Si UB from two reflections (issue text setup)."""
    g = ahd.make_geometry("psic")
    g.sample.lattice = ahd.Lattice(a=5.431)
    g.wavelength = 1.0
    r1 = g.add_reflection(
        "r1",
        (0, 0, 1),
        {"mu": 0, "eta": 22, "chi": 93, "phi": 0, "nu": 0, "delta": 40},
        wavelength=1.0,
    )
    r2 = g.add_reflection(
        "r2",
        (1, 0, 0),
        {"mu": 0, "eta": 22, "chi": 180, "phi": 0, "nu": 0, "delta": 40},
        wavelength=1.0,
    )
    ahd.ub_from_two_reflections_bl1967(g.sample, r1, r2)
    return g


def _set_incidence(g, value):
    g.mode_name = "fixed_incidence_vertical"
    g.surface_normal = (0, 0, 1)
    g._modes[g.mode_name] = g.mode.with_constraint_values(incidence=value)
    g.mode_name = g.mode_name


@pytest.mark.parametrize(
    "hkl, incidence, context",
    [
        pytest.param((0, 1, 1), 2.0, does_not_raise(), id="011-inc2"),
        pytest.param((1, 0, 1), 1.0, does_not_raise(), id="101-inc1"),
    ],
)
def test_fixed_incidence_solutions_round_trip(hkl, incidence, context):
    """Each fixed_incidence_vertical solution reproduces hkl and honors incidence."""
    with context:
        g = _silicon_psic()
        _set_incidence(g, incidence)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sols = g.forward(*hkl)

        assert sols, "expected at least one accessible solution"
        for sol in sols:
            rt = g.inverse(sol)
            assert rt == pytest.approx(tuple(float(v) for v in hkl), abs=1e-3)
            assert incidence_angle(g, angles=sol) == pytest.approx(incidence, abs=1e-4)


def test_fixed_incidence_inaccessible_returns_empty():
    """An inaccessible incidence target returns [] (not wrong-hkl settings)."""
    g = _silicon_psic()
    _set_incidence(g, 89.0)  # too steep to be reachable for (0, 1, 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sols = g.forward(0, 1, 1)
    assert sols == []
