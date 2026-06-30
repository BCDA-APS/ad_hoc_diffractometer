# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
Regression tests for issue #304.

Issue #304 reported that ``ReferenceConstraint.evaluate()`` /
``.is_satisfied()`` were stubs raising ``NotImplementedError``, with code
comments pointing at the long-closed "Issue J" (#157).  In fact the
forward solvers for every reference-constraint mode
(``incidence`` / ``emergence`` / ``specular`` / ``psi`` / ``omega``)
were already implemented in :mod:`ad_hoc_diffractometer.forward`; the
dispatcher routes those modes to dedicated ``_solve_*`` functions and
never called the stubbed methods.

The fix makes ``ReferenceConstraint.evaluate()`` /
``.is_satisfied()`` compute real pseudo-angle residuals using
:mod:`ad_hoc_diffractometer.reference`, and wires the reference
constraints (except ``psi``, which is enforced upstream as a validation
filter) into the post-solve verification loop in
:func:`ad_hoc_diffractometer.forward._validate_solutions`.

These tests span the ``mode``, ``forward``, and ``reference`` modules,
so they live in a cross-module regression file per the project's testing
conventions.

The tests below confirm:

* every affected reference-constraint mode solves end-to-end through
  ``forward()`` once its reference vector is set (no spurious
  ``ConstraintViolation`` from the verification loop), and
* the verification loop rejects a solution whose reference residual
  exceeds tolerance.
"""

from __future__ import annotations

from contextlib import nullcontext as does_not_raise

import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import ConstraintViolation
from ad_hoc_diffractometer import ReferenceConstraint
from ad_hoc_diffractometer.forward import _validate_solutions

WAVELENGTH = 1.5406  # Cu Kα


def _setup_cubic(name, a=4.0, **kwargs):
    """Return a cubic geometry ready for forward()."""
    g = ahd.make_geometry(name)
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=a)
    ahd.ub_identity(g.sample)
    for key, val in kwargs.items():
        setattr(g, key, val)
    return g


# ---------------------------------------------------------------------------
# Every affected reference-constraint mode solves end-to-end.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "geometry, mode_name, ref_attr, ref_val, context",
    [
        pytest.param(
            "psic",
            "fixed_incidence_vertical",
            "surface_normal",
            (0, 0, 1),
            does_not_raise(),
            id="psic-incidence-vertical",
        ),
        pytest.param(
            "psic",
            "fixed_emergence_vertical",
            "surface_normal",
            (0, 0, 1),
            does_not_raise(),
            id="psic-emergence-vertical",
        ),
        pytest.param(
            "psic",
            "specular_vertical",
            "surface_normal",
            (0, 0, 1),
            does_not_raise(),
            id="psic-specular-vertical",
        ),
        pytest.param(
            "psic",
            "specular_horizontal",
            "surface_normal",
            (0, 0, 1),
            does_not_raise(),
            id="psic-specular-horizontal",
        ),
        pytest.param(
            "psic",
            "fixed_omega_vertical",
            "surface_normal",
            (0, 0, 1),
            does_not_raise(),
            id="psic-omega-vertical",
        ),
        pytest.param(
            "sixc",
            "fixed_incidence_zaxis",
            "surface_normal",
            (0, 0, 1),
            does_not_raise(),
            id="sixc-incidence-zaxis",
        ),
        pytest.param(
            "sixc",
            "specular_zaxis",
            "surface_normal",
            (0, 0, 1),
            does_not_raise(),
            id="sixc-specular-zaxis",
        ),
        pytest.param(
            "zaxis",
            "zaxis",
            "surface_normal",
            (0, 0, 1),
            does_not_raise(),
            id="zaxis-zaxis",
        ),
        pytest.param(
            "zaxis",
            "reflectivity",
            "surface_normal",
            (0, 0, 1),
            does_not_raise(),
            id="zaxis-reflectivity",
        ),
        pytest.param(
            "s2d2",
            "reflectivity",
            "surface_normal",
            (0, 0, 1),
            does_not_raise(),
            id="s2d2-reflectivity",
        ),
    ],
)
def test_reference_mode_solves_without_spurious_violation(
    geometry, mode_name, ref_attr, ref_val, context
):
    """Reference-constraint modes solve end-to-end with verification active.

    The post-solve verification loop now evaluates the reference-constraint
    residual; a correct solver must produce solutions that pass it without
    raising :class:`ConstraintViolation`.
    """
    with context:
        g = _setup_cubic(geometry, **{ref_attr: ref_val})
        g.mode_name = mode_name
        assert g.mode.is_implemented(g)
        # The mode default emergence target (0°) is physically inaccessible
        # for cubic (1, 0, 1) with surface_normal (0, 0, 1): the achievable
        # emergence range is ~0.7°–56°, so emergence=0 yields no Bragg-valid
        # solution.  Use an accessible target so the test exercises the
        # verification loop on real solutions rather than (formerly) on
        # the buggy wrong-hkl output (issues #304, #307).
        if mode_name == "fixed_emergence_vertical":
            g._modes[mode_name] = g.mode.with_constraint_values(  # noqa: SLF001
                emergence=5.0
            )
            g.mode_name = mode_name
        sols = g.forward(1, 0, 1)
        assert len(sols) > 0
        # Re-running the verification loop explicitly must not raise.
        _validate_solutions(sols, g.mode, g)
        # Every returned solution satisfies the reference constraint.
        rc = next(c for c in g.mode.constraints if c.category == "reference")
        if rc.name != "psi" and rc.has_reference_vector(g):
            for sol in sols:
                assert abs(rc.evaluate(sol, g)) < 1e-6


# ---------------------------------------------------------------------------
# The verification loop rejects a bad reference solution.
# ---------------------------------------------------------------------------


def test_verification_rejects_bad_reference_solution():
    """A solution violating a reference constraint raises ConstraintViolation."""
    g = _setup_cubic("psic", surface_normal=(0, 0, 1))
    g.mode_name = "fixed_incidence_vertical"
    good = g.forward(1, 0, 1)[0]

    # A mode demanding incidence=45° is not met by the incidence=0 solution.
    bad_mode = ahd.ConstraintSet(
        [ReferenceConstraint("incidence", 45.0)],
        extras={"n_hat": ahd.REQUIRED, "incidence": None},
    )
    with pytest.raises(ConstraintViolation):
        _validate_solutions([good], bad_mode, g)


# ---------------------------------------------------------------------------
# psi is skipped by the verification loop (enforced upstream as a filter).
# ---------------------------------------------------------------------------


def test_verification_skips_reference_without_vector():
    """A reference constraint whose vector is unset is skipped, not raised."""
    g = _setup_cubic("psic")  # no surface_normal set
    assert g.surface_normal is None
    mode = ahd.ConstraintSet(
        [ReferenceConstraint("incidence", 0.0)],
        extras={"n_hat": ahd.REQUIRED, "incidence": None},
    )
    # Any motor-angle dict works; the guard short-circuits before evaluate().
    angles = {s.name: s.angle for s in g._stages.values()}  # noqa: SLF001
    _validate_solutions([angles], mode, g)  # must not raise


def test_verification_swallows_reference_evaluate_error(monkeypatch):
    """If evaluate() raises despite a set vector, the loop skips that solution."""
    g = _setup_cubic("psic", surface_normal=(0, 0, 1))
    mode = ahd.ConstraintSet(
        [ReferenceConstraint("incidence", 0.0)],
        extras={"n_hat": ahd.REQUIRED, "incidence": None},
    )
    angles = {s.name: s.angle for s in g._stages.values()}  # noqa: SLF001

    def _boom(*_args, **_kwargs):
        raise ValueError("synthetic incidence failure for issue #304 coverage")

    monkeypatch.setattr("ad_hoc_diffractometer.reference.incidence_angle", _boom)
    _validate_solutions([angles], mode, g)  # must not raise


def test_psi_constraint_skipped_by_verification_loop():
    """psi residuals (subject to ±360° wrap) do not trip the verification loop."""
    g = _setup_cubic("fourcv", azimuth=(0, 0, 1))
    import numpy as np

    from ad_hoc_diffractometer.forward import _compute_natural_psi

    Q_phi = g.sample.UB @ np.array([1.0, 1.0, 1.0])
    natural = _compute_natural_psi(g, Q_phi)
    assert natural is not None

    g.modes["fixed_psi"] = ahd.ConstraintSet(
        [ReferenceConstraint("psi", natural)],
        computed=g.modes["fixed_psi"].computed,
        extras={"n_hat": ahd.REQUIRED, "psi": None},
    )
    g.mode_name = "fixed_psi"
    sols = g.forward(1, 1, 1)
    assert len(sols) > 0
    # Must not raise even though a naive psi residual could be ~360°.
    _validate_solutions(sols, g.mode, g)
