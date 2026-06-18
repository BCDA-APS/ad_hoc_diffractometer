# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
Regression tests for issue #279.

Before the fix, the surface-mode solver
:func:`ad_hoc_diffractometer.forward._solve_surface` wrote the
``2θ`` angle to the **last** stage in ``geometry.detector_stages``
unconditionally.  For psic that last stage is ``delta``, and the
horizontal surface modes
(``fixed_incidence_horizontal``, ``fixed_emergence_horizontal``,
``specular_horizontal``) declare a
:class:`~ad_hoc_diffractometer.mode.DetectorConstraint` that pins
``delta = 0``.  The solver therefore overwrote the pinned value
moments after applying it and left the truly active detector stage
``nu`` at zero — producing "solutions" with ``delta = 2θ`` and
``nu = 0`` that violated the mode's own constraint.

The mirror failure occurred in the vertical surface modes
(``fixed_incidence_vertical``, ``fixed_emergence_vertical``,
``specular_vertical``), where the mode pins ``nu = 0`` and
the active detector stage is ``delta`` — but the dispatch picked
``delta`` for both roles regardless, so the constraint happened to
agree with the active stage by accident (the resulting ``nu = 0``
matched the pinned value), masking the underlying dispatch bug on
that axis.

The fix selects the active detector stage by *excluding* the
detector stage named by the mode's non-qaz
:class:`~ad_hoc_diffractometer.mode.DetectorConstraint`.  These
tests pin the contract:

* Every horizontal surface mode on psic produces solutions that
  satisfy ``delta == 0`` and place ``2θ`` on ``nu``.
* Every vertical surface mode on psic produces solutions that
  satisfy ``nu == 0`` and place ``2θ`` on ``delta``.
* The legacy single-detector-pin-free surface geometries (zaxis,
  sixc) continue to dispatch ``2θ`` to the last detector stage as
  before — this guards against an over-broad fix.

These checks intentionally do **not** assert ``forward → inverse``
round-trip on every psic surface-mode hkl: ``_solve_surface``
rocks a single sample stage at a time, and psic has four sample
stages, so the legacy 1-D solver cannot in general reach an
arbitrary hkl from this mode family.  That orthogonal limitation
(routing psic surface modes to a multi-sample-stage solver) is a
separate concern from the detector-dispatch bug fixed here.
"""

from __future__ import annotations

from contextlib import nullcontext as does_not_raise

import pytest
from helpers import psic
from helpers import sixc
from helpers import zaxis

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import ub_identity

WAVELENGTH = 1.5406  # Cu Kα

# Surface-reference modes on psic whose DetectorConstraint pins the
# inner (delta) stage so the active detector for 2θ must be ``nu``.
_PSIC_HORIZONTAL_SURFACE_MODES = (
    "fixed_incidence_horizontal",
    "fixed_emergence_horizontal",
    "specular_horizontal",
)

# Surface-reference modes on psic whose DetectorConstraint pins the
# outer (nu) stage so the active detector for 2θ must be ``delta``.
_PSIC_VERTICAL_SURFACE_MODES = (
    "fixed_incidence_vertical",
    "fixed_emergence_vertical",
    "specular_vertical",
)


def _setup_psic_cubic(a: float = 4.0):
    """psic with UB = B for a cubic lattice and surface normal +z."""
    g = psic()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=a)
    ub_identity(g.sample)
    g.surface_normal = (0, 0, 1)
    return g


# ---------------------------------------------------------------------------
# Horizontal surface modes: delta must stay at 0, nu must carry 2θ.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, context",
    [pytest.param(m, does_not_raise(), id=m) for m in _PSIC_HORIZONTAL_SURFACE_MODES],
)
def test_psic_horizontal_surface_honors_delta_pin(mode_name, context):
    """Horizontal surface modes must leave ``delta = 0`` and put 2θ on ``nu``.

    The mode's ``DetectorConstraint('delta', 0.0)`` is applied first;
    the active detector stage that receives the 2θ magnitude is
    therefore the *other* detector stage (``nu``).
    """
    with context:
        g = _setup_psic_cubic()
        g.mode_name = mode_name
        sols = g.forward(1, 0, 0)
        assert sols, f"{mode_name}: expected at least one candidate solution"
        for sol in sols:
            assert sol["delta"] == pytest.approx(0.0, abs=1e-8), (
                f"{mode_name}: delta pin violated, delta={sol['delta']}, nu={sol['nu']}"
            )
            # The 2θ magnitude for (1,0,0) on a=4 Å with λ=1.5406 Å is
            # 2·arcsin(λ / 2a) ≈ 22.206°.  nu must be at ±2θ (the sign
            # depends on the sample-stage seeding chosen by the Newton
            # search).
            assert abs(sol["nu"]) == pytest.approx(22.2062, abs=1e-3), (
                f"{mode_name}: nu should carry the 2θ magnitude, "
                f"nu={sol['nu']}, delta={sol['delta']}"
            )


# ---------------------------------------------------------------------------
# Vertical surface modes: nu must stay at 0, delta must carry 2θ.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, context",
    [pytest.param(m, does_not_raise(), id=m) for m in _PSIC_VERTICAL_SURFACE_MODES],
)
def test_psic_vertical_surface_honors_nu_pin(mode_name, context):
    """Vertical surface modes must leave ``nu = 0`` and put 2θ on ``delta``.

    The mode's ``DetectorConstraint('nu', 0.0)`` is applied first;
    the active detector stage is therefore ``delta``.  Before the
    fix the dispatch picked ``delta`` regardless of the pin, so the
    vertical case appeared to "work" because the pin happened to
    name the other stage — but the dispatch was still hard-coded
    rather than constraint-aware.
    """
    with context:
        g = _setup_psic_cubic()
        g.mode_name = mode_name
        sols = g.forward(1, 0, 0)
        assert sols, f"{mode_name}: expected at least one candidate solution"
        for sol in sols:
            assert sol["nu"] == pytest.approx(0.0, abs=1e-8), (
                f"{mode_name}: nu pin violated, nu={sol['nu']}, delta={sol['delta']}"
            )
            assert abs(sol["delta"]) == pytest.approx(22.2062, abs=1e-3), (
                f"{mode_name}: delta should carry the 2θ magnitude, "
                f"delta={sol['delta']}, nu={sol['nu']}"
            )


# ---------------------------------------------------------------------------
# Legacy surface geometries (no chi): the fix must not disturb them.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "geometry_factory, mode_name, h, k, l, context",
    [
        # sixc surface-zaxis modes: detector pin freezes ``gamma``;
        # remaining detector stage ``delta`` carries 2θ.  The fix
        # must not alter the existing healthy round-trip count.
        pytest.param(
            sixc,
            "fixed_incidence_zaxis",
            1,
            0,
            0,
            does_not_raise(),
            id="sixc-fixed_incidence_zaxis-100",
        ),
        pytest.param(
            sixc,
            "specular_zaxis",
            1,
            0,
            0,
            does_not_raise(),
            id="sixc-specular_zaxis-100",
        ),
        # zaxis reflectivity mode: confirm at least one returned
        # solution still round-trips, i.e. the legacy path is intact.
        pytest.param(
            zaxis,
            "reflectivity",
            1,
            0,
            0,
            does_not_raise(),
            id="zaxis-reflectivity-100",
        ),
    ],
)
def test_legacy_surface_geometries_unchanged(
    geometry_factory,
    mode_name,
    h,
    k,
    l,  # noqa: E741
    context,
):
    """Non-psic surface geometries must continue to round-trip.

    The fix narrows the active-detector dispatch only when the
    mode's DetectorConstraint pins a real detector stage by name.
    These geometries either have a single detector stage or pin
    the outer stage in a way that already coincides with the legacy
    behavior; the regression catch is that *at least one* returned
    solution round-trips to the requested hkl, which is the
    pre-existing baseline.
    """
    with context:
        g = geometry_factory()
        g.wavelength = WAVELENGTH
        g.sample.lattice = ahd.Lattice(a=4.0)
        ub_identity(g.sample)
        g.surface_normal = (0, 0, 1)
        g.mode_name = mode_name

        sols = g.forward(h, k, l)
        assert sols, f"{mode_name}: expected at least one solution"
        ok = 0
        for sol in sols:
            rt = g.inverse(sol)
            if all(abs(a - b) < 1e-3 for a, b in zip(rt, (h, k, l), strict=False)):
                ok += 1
        assert ok >= 1, (
            f"{mode_name}: expected at least one round-trip-correct "
            f"solution, got 0 of {len(sols)}"
        )


# ---------------------------------------------------------------------------
# Direct dispatch check: surface solver returns [] when no detector free.
# ---------------------------------------------------------------------------


def test_solve_surface_returns_empty_when_only_detector_stage_pinned():
    """``_solve_surface`` must return ``[]`` on a single-detector-stage
    geometry whose only detector stage is pinned by the mode.

    Before the fix the solver would have happily overwritten the only
    detector stage with ``2θ``, silently violating the user-declared
    constraint.  After the fix it correctly identifies that there is
    no active stage and returns an empty solution list, so the caller
    can decide what to do.

    The :class:`~ad_hoc_diffractometer.mode.ConstraintSet` validator
    forbids two :class:`~ad_hoc_diffractometer.mode.DetectorConstraint`
    instances in the same mode, so this regression uses ``fourcv``
    (which has a single ``ttheta`` detector stage) and a synthetic
    mode that pins that stage.
    """
    import numpy as np

    from ad_hoc_diffractometer.forward import _solve_surface
    from ad_hoc_diffractometer.mode import ConstraintSet
    from ad_hoc_diffractometer.mode import DetectorConstraint
    from ad_hoc_diffractometer.mode import ReferenceConstraint
    from ad_hoc_diffractometer.mode import SampleConstraint

    g = ahd.make_geometry("fourcv")
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=4.0)
    ub_identity(g.sample)
    g.surface_normal = (0, 0, 1)
    det_name = g.detector_stages[-1].name  # the only detector stage

    over_constrained = ConstraintSet(
        constraints=[
            SampleConstraint("omega", 0.0),
            DetectorConstraint(det_name, 0.0),
            ReferenceConstraint("alpha_i", 0.0),
        ],
    )

    result = _solve_surface(g, np.array([0.0, 0.0, 1.0]), 30.0, over_constrained)
    assert result == []
