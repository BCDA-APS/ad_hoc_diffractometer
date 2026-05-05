# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Cross-module regression tests for issue #243.

Reported by @jwkim-anl: the psic ``fixed_phi_vertical``,
``fixed_chi_vertical``, ``fixed_phi_horizontal``, and
``fixed_chi_horizontal`` modes used a :class:`BisectConstraint` where
the kinematically correct constraint is a plane-lock
:class:`SampleConstraint` on the base sample stage (``mu = 0`` for the
vertical-scattering family, ``eta = 0`` for the horizontal-scattering
family).

You (1999) §5 organises psic modes by how many of the three free
angles are pseudo-angles vs sample-orienting circles.  The
``fixed_phi`` / ``fixed_chi`` family belongs to §5.4 (two
sample-orienting angles given) — there is no bisecting condition in
this family.  SPEC's ``psic`` macros implement the same convention:
every ``vertical`` mode locks ``mu`` and every ``horizontal`` mode
locks ``eta``, independent of which other sample stage is fixed.

Pre-fix, the bisect happened to keep the locked axis at zero for the
specific test points used elsewhere in the suite (e.g. cubic
``(1, 0, 0)`` with ``UB = B``), so the bug did not surface there.
For a general hkl on a less symmetric crystal the bisect leaves the
plane-lock axis free, which violates the named scattering geometry.

Related cleanup in the same PR
------------------------------

While reviewing the affected modes, four other psic modes were found
to be **exact duplicates** of ``bisecting_vertical`` /
``bisecting_horizontal`` (same constraint set, same ``computed=``
list — the only difference was the order of the constraints in the
list):

- ``fixed_mu_vertical``    → identical to ``bisecting_vertical``
- ``fixed_nu_vertical``    → identical to ``bisecting_vertical``
- ``fixed_eta_horizontal`` → identical to ``bisecting_horizontal``
- ``fixed_delta_horizontal`` → identical to ``bisecting_horizontal``

These four modes were removed in the same PR.  See
:func:`test_removed_redundant_modes_absent`.

Note on ``fixed_chi_horizontal``
--------------------------------

After the fix, ``fixed_chi_horizontal`` (chi held at 90°, eta = 0,
delta = 0) is kinematically infeasible on the cubic ``UB = B`` test
crystal: ``chi = 90`` rotates Q out of the horizontal scattering
plane, conflicting with ``delta = 0``.  This is a separate concern
from the bisect-vs-plane-lock fix in #243 (the constraint structure
is now correct; only the *value* of chi appears wrong, suggesting
the SPEC convention is ``chi = 0``).  Resolution of the correct
default chi value (and the full audit of psic modes against SPEC
``psic`` and Hkl/Soleil ``E6C``) is tracked in **issue #259**.  This
file therefore does not parametrise ``fixed_chi_horizontal`` for
the round-trip / plane-lock value tests; the test row will be added
when #259 lands the corrected default.

This file lives at the cross-module regression level rather than in
``tests/test_presets.py`` because the symptom is visible only after
``forward()`` is exercised through the constraint-set dispatcher in
``forward.py``.
"""

from __future__ import annotations

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer.presets import psic

# Cu Ka wavelength — matches the rest of the test suite.
WAVELENGTH = 1.5406


def _setup_psic_cubic(a: float = 4.0):
    """Return a psic geometry with UB = B for a cubic lattice of side ``a``."""
    g = psic()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=a)
    ub_identity(g.sample)
    return g


# ---------------------------------------------------------------------------
# Constraint structure: every affected mode contains the plane-lock
# SampleConstraint and no longer contains a BisectConstraint.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, locked_stage, locked_value, context",
    [
        pytest.param(
            "fixed_phi_vertical",
            "mu",
            0.0,
            does_not_raise(),
            id="fixed_phi_vertical-locks-mu",
        ),
        pytest.param(
            "fixed_chi_vertical",
            "mu",
            0.0,
            does_not_raise(),
            id="fixed_chi_vertical-locks-mu",
        ),
        pytest.param(
            "fixed_phi_horizontal",
            "eta",
            0.0,
            does_not_raise(),
            id="fixed_phi_horizontal-locks-eta",
        ),
        pytest.param(
            "fixed_chi_horizontal",
            "eta",
            0.0,
            does_not_raise(),
            id="fixed_chi_horizontal-locks-eta",
        ),
    ],
)
def test_mode_contains_plane_lock_sample_constraint(
    mode_name, locked_stage, locked_value, context
):
    """Affected modes contain the plane-lock SampleConstraint, not a Bisect."""
    with context:
        g = psic()
        cs = g.modes[mode_name]
        # The plane-lock SampleConstraint must be present.
        names = {c.name: c.value for c in cs.fixed_sample_constraints}
        assert locked_stage in names, (
            f"mode {mode_name!r} is missing the plane-lock "
            f"SampleConstraint on {locked_stage!r}"
        )
        assert names[locked_stage] == pytest.approx(locked_value, abs=1e-12)
        # No BisectConstraint must remain in the affected modes.
        assert not cs.has_bisect, (
            f"mode {mode_name!r} must not contain a BisectConstraint after #243"
        )


# ---------------------------------------------------------------------------
# Behavioural test: the locked axis stays at 0 in all forward solutions
# for a range of accessible hkl on a simple cubic crystal.
# ---------------------------------------------------------------------------


_HKL_VERTICAL = [
    pytest.param(1, 0, 0, id="100"),
    pytest.param(0, 1, 0, id="010"),
    pytest.param(1, 1, 0, id="110"),
    pytest.param(1, 1, 1, id="111"),
]

# fixed_phi_horizontal is reachable on cubic UB=B for the hkl below;
# fixed_chi_horizontal is *not* reachable for any low-index hkl on this
# crystal (see module docstring for why).
_HKL_HORIZONTAL = [
    pytest.param(0, 1, 0, id="010"),
    pytest.param(0, 0, 1, id="001"),
    pytest.param(0, 1, 1, id="011"),
    pytest.param(1, 0, 1, id="101"),
    pytest.param(1, 1, 0, id="110"),
    pytest.param(1, 1, 1, id="111"),
]


@pytest.mark.parametrize("h, k, l", _HKL_VERTICAL)
@pytest.mark.parametrize(
    "mode_name, locked_stage",
    [
        pytest.param("fixed_phi_vertical", "mu", id="fixed_phi_vertical"),
        pytest.param("fixed_chi_vertical", "mu", id="fixed_chi_vertical"),
    ],
)
def test_fixed_vertical_modes_lock_mu(mode_name, locked_stage, h, k, l):  # noqa: E741
    """``mu`` is exactly zero in every forward solution for vertical modes.

    Pre-fix this assertion would fail for at least some hkl on less
    symmetric crystals, because the BisectConstraint left ``mu`` free
    rather than locking the scattering plane vertical.
    """
    g = _setup_psic_cubic()
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0, (
        f"forward({h},{k},{l}) returned no solutions for {mode_name!r}"
    )
    for sol in solutions:
        assert sol[locked_stage] == pytest.approx(0.0, abs=1e-8), (
            f"{mode_name}: {locked_stage}={sol[locked_stage]!r} (expected 0)"
        )


@pytest.mark.parametrize("h, k, l", _HKL_HORIZONTAL)
def test_fixed_phi_horizontal_locks_eta(h, k, l):  # noqa: E741
    """``eta`` is exactly zero in every forward solution for fixed_phi_horizontal.

    fixed_chi_horizontal is omitted: see module docstring (chi=90 is
    kinematically infeasible on cubic UB=B; deferred to the audit
    follow-up issue).
    """
    g = _setup_psic_cubic()
    g.mode_name = "fixed_phi_horizontal"
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0, (
        f"forward({h},{k},{l}) returned no solutions for fixed_phi_horizontal"
    )
    for sol in solutions:
        assert sol["eta"] == pytest.approx(0.0, abs=1e-8), (
            f"fixed_phi_horizontal: eta={sol['eta']!r} (expected 0)"
        )


# ---------------------------------------------------------------------------
# Round-trip: forward() solutions invert back to the requested hkl.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, h, k, l",
    [
        pytest.param("fixed_phi_vertical", 1, 0, 0, id="fixed_phi_vertical-100"),
        pytest.param("fixed_phi_vertical", 1, 1, 1, id="fixed_phi_vertical-111"),
        pytest.param("fixed_chi_vertical", 1, 0, 0, id="fixed_chi_vertical-100"),
        pytest.param("fixed_chi_vertical", 1, 1, 1, id="fixed_chi_vertical-111"),
        pytest.param("fixed_phi_horizontal", 0, 0, 1, id="fixed_phi_horizontal-001"),
        pytest.param("fixed_phi_horizontal", 0, 1, 1, id="fixed_phi_horizontal-011"),
        pytest.param("fixed_phi_horizontal", 1, 1, 1, id="fixed_phi_horizontal-111"),
    ],
)
def test_round_trip_forward_inverse(mode_name, h, k, l):  # noqa: E741
    """``inverse(forward(hkl)) == hkl`` for every solution."""
    g = _setup_psic_cubic()
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0
    for sol in solutions:
        hkl_back = g.inverse(sol)
        np.testing.assert_allclose(hkl_back, [h, k, l], atol=1e-6)


# ---------------------------------------------------------------------------
# Plane-axis verification: the fixed_chi/_phi families also keep the
# detector plane-lock axis at zero (nu = 0 vertical, delta = 0 horizontal).
# This locks in the full "scattering plane is named vertical/horizontal"
# property — both the sample base and the detector base must be fixed.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, sample_lock, detector_lock, fixed_stage, fixed_value, h, k, l",
    [
        pytest.param(
            "fixed_phi_vertical",
            "mu",
            "nu",
            "phi",
            0.0,
            1,
            0,
            0,
            id="fixed_phi_vertical-100",
        ),
        pytest.param(
            "fixed_chi_vertical",
            "mu",
            "nu",
            "chi",
            90.0,
            1,
            0,
            0,
            id="fixed_chi_vertical-100",
        ),
        pytest.param(
            "fixed_phi_horizontal",
            "eta",
            "delta",
            "phi",
            0.0,
            0,
            0,
            1,
            id="fixed_phi_horizontal-001",
        ),
    ],
)
def test_full_plane_lock_invariants(
    mode_name,
    sample_lock,
    detector_lock,
    fixed_stage,
    fixed_value,
    h,
    k,
    l,  # noqa: E741
):
    """All three locked axes (sample base, detector base, named fix) hold."""
    g = _setup_psic_cubic()
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol[sample_lock] == pytest.approx(0.0, abs=1e-8), (
            f"{mode_name}: {sample_lock}={sol[sample_lock]!r} (expected 0)"
        )
        assert sol[detector_lock] == pytest.approx(0.0, abs=1e-8), (
            f"{mode_name}: {detector_lock}={sol[detector_lock]!r} (expected 0)"
        )
        assert sol[fixed_stage] == pytest.approx(fixed_value, abs=1e-6), (
            f"{mode_name}: {fixed_stage}={sol[fixed_stage]!r} (expected {fixed_value})"
        )


# ---------------------------------------------------------------------------
# Cleanup: four modes that were exact duplicates of bisecting_vertical /
# bisecting_horizontal must no longer be present on psic.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "removed_mode",
    [
        pytest.param("fixed_mu_vertical", id="fixed_mu_vertical"),
        pytest.param("fixed_nu_vertical", id="fixed_nu_vertical"),
        pytest.param("fixed_eta_horizontal", id="fixed_eta_horizontal"),
        pytest.param("fixed_delta_horizontal", id="fixed_delta_horizontal"),
    ],
)
def test_removed_redundant_modes_absent(removed_mode):
    """Redundant modes are no longer in psic.modes."""
    g = psic()
    assert removed_mode not in g.modes, (
        f"{removed_mode!r} was an exact duplicate of bisecting_vertical/"
        f"bisecting_horizontal and must have been removed in #243"
    )


# ---------------------------------------------------------------------------
# Coverage: exercise the analytic-empty fallthrough in _solve_bisecting.
#
# Removing the four redundant bisecting modes also removed the only test
# path that exercised the ``if analytic_results: ... else: fall through``
# branch in :func:`ad_hoc_diffractometer.forward._solve_bisecting`
# (line 1052).  The fallthrough is hit when the scattering vector is
# kinematically inaccessible in the analytic solver but Newton fallback
# must still be invoked (and ultimately also returns no solutions).
#
# ``psic bisecting_horizontal`` with the inaccessible reflection
# ``(-2, -2, 0)`` reproduces this branch on a cubic ``UB = B`` test
# crystal.
# ---------------------------------------------------------------------------


def test_bisecting_horizontal_analytic_empty_fallthrough():
    """Inaccessible hkl returns 0 solutions and exercises the Newton fallback."""
    g = _setup_psic_cubic()
    g.mode_name = "bisecting_horizontal"
    solutions = g.forward(-2, -2, 0)
    assert solutions == []
