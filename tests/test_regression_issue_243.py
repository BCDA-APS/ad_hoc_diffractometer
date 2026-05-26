# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
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

After the #243 fix the constraint *structure* of
``fixed_chi_horizontal`` was correct (``eta = 0``, ``delta = 0``,
``chi`` fixed) but the *value* of chi was wrong: it carried over the
``chi = 90`` default from ``fixed_chi_vertical`` and from the
four-circle ``fixed_chi`` modes, where 90° is the canonical
"spinning Q" symmetry value.  In the horizontal-scattering psic
sub-geometry ``chi = 90`` instead tilts the phi axis out of the
horizontal plane, conflicting with ``delta = 0`` and rendering every
low-index hkl unreachable on a cubic ``UB = B`` crystal.

Issue #259 audited the psic modes against SPEC ``psic`` and Hkl/Soleil
``E6C`` and resolved the default to ``chi = 0`` — the value that keeps
the phi axis in the horizontal scattering plane (the chi-circle axis
in the residual sub-geometry lies along the longitudinal/beam
direction with ``eta = 0``).  This file now parametrises
``fixed_chi_horizontal`` for the eta-lock, round-trip and full
plane-lock invariant tests using hkl that are reachable with the new
default.

This file lives at the cross-module regression level rather than in
``tests/test_presets.py`` because the symptom is visible only after
``forward()`` is exercised through the constraint-set dispatcher in
``forward.py``.
"""

from __future__ import annotations

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from helpers import psic

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import ub_identity

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

# fixed_phi_horizontal is reachable on cubic UB=B for the hkl below.
_HKL_HORIZONTAL = [
    pytest.param(0, 1, 0, id="010"),
    pytest.param(0, 0, 1, id="001"),
    pytest.param(0, 1, 1, id="011"),
    pytest.param(1, 0, 1, id="101"),
    pytest.param(1, 1, 0, id="110"),
    pytest.param(1, 1, 1, id="111"),
]

# fixed_chi_horizontal at the new (#259) default chi = 0 is reachable
# only for hkl with h = 0 on cubic UB = B, because chi = 0 leaves the
# phi-circle axis along the longitudinal (beam) direction so phi
# rotations cannot move the h-component of Q out of the beam plane.
_HKL_FIXED_CHI_HORIZONTAL = [
    pytest.param(0, 1, 0, id="010"),
    pytest.param(0, 0, 1, id="001"),
    pytest.param(0, 1, 1, id="011"),
    pytest.param(0, 0, 2, id="002"),
    pytest.param(0, 2, 1, id="021"),
    pytest.param(0, 1, 2, id="012"),
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
    """``eta`` is exactly zero in every forward solution for fixed_phi_horizontal."""
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


@pytest.mark.parametrize("h, k, l", _HKL_FIXED_CHI_HORIZONTAL)
def test_fixed_chi_horizontal_locks_eta(h, k, l):  # noqa: E741
    """``eta`` is exactly zero in every forward solution for fixed_chi_horizontal.

    Resolved in issue #259: the default chi value was corrected from
    90° (the four-circle "spinning Q" value, kinematically infeasible
    here) to 0° (which keeps the phi axis in the horizontal scattering
    plane).  Reachable hkl on cubic UB = B all have h = 0; see the
    ``_HKL_FIXED_CHI_HORIZONTAL`` list.
    """
    g = _setup_psic_cubic()
    g.mode_name = "fixed_chi_horizontal"
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0, (
        f"forward({h},{k},{l}) returned no solutions for fixed_chi_horizontal"
    )
    for sol in solutions:
        assert sol["eta"] == pytest.approx(0.0, abs=1e-8), (
            f"fixed_chi_horizontal: eta={sol['eta']!r} (expected 0)"
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
        pytest.param("fixed_chi_horizontal", 0, 1, 0, id="fixed_chi_horizontal-010"),
        pytest.param("fixed_chi_horizontal", 0, 0, 1, id="fixed_chi_horizontal-001"),
        pytest.param("fixed_chi_horizontal", 0, 1, 1, id="fixed_chi_horizontal-011"),
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
        pytest.param(
            "fixed_chi_horizontal",
            "eta",
            "delta",
            "chi",
            0.0,
            0,
            0,
            1,
            id="fixed_chi_horizontal-001",
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
# (line 1052).  The fallthrough is hit when the analytic solver returns
# ``[]`` and the Newton fallback must still be invoked.
#
# Under the issue #280 corrected outermost-leftmost composition the
# analytic solver returns ``[]`` whenever ``Q_phi`` is parallel to the
# phi-axis (phi-rotation degenerate; analytic decomposition cannot
# enumerate the multiple phi-equivalent solutions).  The Newton solver
# then generates the phi-representatives.  Reflection ``(0, 0, 1)`` on
# psic bisecting_horizontal puts ``Q_phi`` along ``+z`` while the
# phi-axis is ``-z``, satisfying this condition.
# ---------------------------------------------------------------------------


def test_bisecting_horizontal_analytic_empty_fallthrough():
    """Analytic returns empty for phi-degenerate target; Newton fills in.

    The historic assertion (``solutions == []``) used reflection
    ``(-2, -2, 0)``, which under the pre-#280 buggy composition was
    geometrically inaccessible.  Under the corrected composition that
    reflection is reachable; no kinematically reachable hkl yields the
    "analytic empty AND Newton empty" combination on this geometry.

    The replacement assertion still exercises the fallthrough branch
    (verified via a monkey-patched spy in
    ``test_bisecting_horizontal_analytic_returns_empty``) and confirms
    the Newton fallback finishes with correct, round-tripping
    solutions.
    """
    g = _setup_psic_cubic()
    g.mode_name = "bisecting_horizontal"
    solutions = g.forward(0, 0, 1)
    assert len(solutions) >= 1
    for sol in solutions:
        hkl_back = g.inverse(sol)
        assert abs(hkl_back[0]) < 1e-6
        assert abs(hkl_back[1]) < 1e-6
        assert abs(hkl_back[2] - 1.0) < 1e-6


def test_bisecting_horizontal_analytic_returns_empty():
    """Verify the analytic solver returns ``[]`` on the phi-degenerate target."""
    from ad_hoc_diffractometer import forward as fmod

    g = _setup_psic_cubic()
    g.mode_name = "bisecting_horizontal"

    captured: list[list[tuple[float, float]]] = []
    orig = fmod._solve_bisecting_analytic

    def spy(ctx, chi_stage, phi_stage, angles, Q_phi_target):
        result = orig(ctx, chi_stage, phi_stage, angles, Q_phi_target)
        captured.append(list(result))
        return result

    fmod._solve_bisecting_analytic = spy
    try:
        g.forward(0, 0, 1)
    finally:
        fmod._solve_bisecting_analytic = orig

    # The analytic solver must have been invoked at least once and
    # returned ``[]`` (phi-degenerate target → defer to Newton fallback).
    assert captured, "analytic solver was never called"
    assert captured[0] == []
