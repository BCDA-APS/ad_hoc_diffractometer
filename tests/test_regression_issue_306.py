# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
Regression tests for issue #306.

Issue #306: for geometry ``psic``, mode ``fixed_psi_vertical``,
``forward(h, k, l)`` returned an empty solution list even when the ``psi``
constraint target equaled the reflection's natural ψ (so the ψ validation
filter passed and no ``UserWarning`` was emitted).  The same reflection
solved normally in ``bisecting_vertical``.

Root cause: once ψ is validated it imposes no further restriction (every
Bragg solution shares the natural ψ for a given hkl + UB), so the psic
``fixed_psi_*`` mode is degenerate in its third free sample angle.  The
solver delegated to ``_solve_fixed_sample``, which froze the outer sample
stage (``eta``) at its current motor value and grid-searched only the
inner ``(chi, phi)`` pair — no ``(chi, phi)`` satisfied Bragg with ``eta``
frozen, so the result was empty.

Fix: ``_solve_psi_mode`` now routes the psic ``fixed_psi_*`` family to a
synthetic bisecting ConstraintSet (the lone non-chi/non-phi free sample
stage takes ttheta/2, paired with the active detector), so
``fixed_psi_vertical`` reproduces ``bisecting_vertical``.
"""

from __future__ import annotations

import warnings
from contextlib import nullcontext as does_not_raise

import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.reference import natural_psi


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


def _set_psi(g, hkl):
    """Select fixed_psi_vertical with the constraint at the natural ψ."""
    g.mode_name = "fixed_psi_vertical"
    g.azimuth = (0, -1, 0)
    nat = natural_psi(g, *hkl)
    g._modes[g.mode_name] = g.mode.with_constraint_values(psi=nat)
    g.mode_name = g.mode_name
    return nat


@pytest.mark.parametrize(
    "hkl, context",
    [
        pytest.param((0, 0, 1), does_not_raise(), id="001"),
        pytest.param((1, 0, 0), does_not_raise(), id="100"),
    ],
)
def test_fixed_psi_vertical_matches_bisecting(hkl, context):
    """fixed_psi_vertical returns the same solution set as bisecting_vertical."""
    with context:
        g = _silicon_psic()
        _set_psi(g, hkl)
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # any ψ filter warning would raise
            psi_sols = g.forward(*hkl)

        g.mode_name = "bisecting_vertical"
        bis_sols = g.forward(*hkl)

        assert len(psi_sols) == len(bis_sols)
        assert len(psi_sols) > 0
        # Each psi-mode solution matches a bisecting solution motor-for-motor.
        for ps in psi_sols:
            assert any(all(abs(ps[k] - bs[k]) < 1e-4 for k in ps) for bs in bis_sols)


def test_fixed_psi_vertical_solutions_round_trip():
    """Every fixed_psi_vertical solution reproduces the requested hkl."""
    g = _silicon_psic()
    _set_psi(g, (0, 0, 1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sols = g.forward(0, 0, 1)
    assert sols
    for sol in sols:
        rt = g.inverse(sol)
        assert rt == pytest.approx((0.0, 0.0, 1.0), abs=1e-3)
