# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
Regression tests for issue #278.

Issue #278 surfaced confusion around the ψ-validation-filter model
(issue #176): when ``forward()`` is called on a ``fixed_psi_*`` mode
with a target ψ that does not match the natural ψ for the requested
(h, k, l), the solver returned a silent empty list — which looked
like a bug to callers expecting a "set ψ and rock the sample to
achieve it" semantics.

In fact ψ is uniquely determined by ``Q_phi = UB @ (h, k, l)`` and
the azimuthal reference vector: every Bragg solution of a fixed
reflection has the same ψ.  No motor configuration can change it.
The empty return is therefore mathematically correct, but the lack
of feedback made the cause non-obvious.

The fix in this PR:

1. Emits a :class:`UserWarning` from
   :func:`ad_hoc_diffractometer.forward._solve_psi_mode` whenever the
   target ψ does not match the natural ψ (or ψ is undefined).  The
   message names the natural value and points the user at
   :func:`~ad_hoc_diffractometer.reference.natural_psi`.

2. Adds the public helper
   :func:`~ad_hoc_diffractometer.reference.natural_psi(g, h, k, l)`
   so callers can discover the natural value programmatically before
   calling :meth:`~diffractometer.AdHocDiffractometer.forward`.

3. Clarifies the validation-filter semantics in the YAML comments
   above every ``fixed_psi*`` mode in the demo geometries.

These tests verify the cross-module contracts:

- ``forward()`` on a ψ-mode with the *natural* ψ returns solutions
  (this is the path the user wants to take after the warning tells
  them the correct value).
- ``forward()`` on a ψ-mode with a *mismatched* ψ emits the warning
  and returns ``[]``.
- ``forward()`` on a reflection with undefined ψ (Q ∥ reference)
  emits the "undefined" warning and returns ``[]``.
- The :func:`reference.natural_psi` helper agrees with the warning
  message and with the psi observed in a real Bragg solution.
"""

from __future__ import annotations

import re
import warnings
from contextlib import nullcontext as does_not_raise

import pytest
from helpers import psic

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer.reference import natural_psi
from ad_hoc_diffractometer.reference import psi_angle

WAVELENGTH = 1.5406  # Cu Kα


def _setup_psic_cubic():
    g = psic()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=4.0)
    ub_identity(g.sample)
    g.azimuth = (0, 0, 1)
    g.surface_normal = (0, 0, 1)
    return g


def _set_psi_constraint(mode, value: float) -> None:
    """Mutate the mode's psi ReferenceConstraint value in place."""
    for c in mode.constraints:
        if c.__class__.__name__ == "ReferenceConstraint" and c.name == "psi":
            c._value = float(value)  # noqa: SLF001
            return
    raise AssertionError("mode has no psi ReferenceConstraint")


# ---------------------------------------------------------------------------
# Mismatched-target path emits warning and returns [].
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, h, k, l, context",
    [
        pytest.param(
            "fixed_psi_horizontal",
            1,
            1,
            0,
            does_not_raise(),
            id="horizontal-110",
        ),
        pytest.param(
            "fixed_psi_vertical",
            1,
            1,
            0,
            does_not_raise(),
            id="vertical-110",
        ),
    ],
)
def test_fixed_psi_mismatch_warns_and_returns_empty(
    mode_name,
    h,
    k,
    l,  # noqa: E741
    context,
):
    """When the constraint ψ does not match the natural ψ for the
    reflection, ``forward()`` warns the user (naming the natural value
    and the recommended helper) and returns ``[]``.
    """
    with context:
        g = _setup_psic_cubic()
        g.mode_name = mode_name
        # Default constraint value is 0.0; the natural ψ for (1,1,0) on
        # cubic / ub_identity / azimuthal=(0,0,1) is 90°, so the target
        # 0° is guaranteed to mismatch.
        nat = natural_psi(g, h, k, l)
        assert nat is not None
        with pytest.warns(UserWarning, match=re.escape(f"{nat:.4f}")):
            solutions = g.forward(h, k, l)
        assert solutions == []


# ---------------------------------------------------------------------------
# Undefined-ψ path emits warning and returns [].
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, h, k, l, context",
    [
        # (0,0,1) makes Q_phi parallel to the azimuthal reference (0,0,1)
        # so ψ is undefined regardless of the mode.
        pytest.param(
            "fixed_psi_horizontal",
            0,
            0,
            1,
            does_not_raise(),
            id="horizontal-undef-001",
        ),
        pytest.param(
            "fixed_psi_vertical",
            0,
            0,
            1,
            does_not_raise(),
            id="vertical-undef-001",
        ),
    ],
)
def test_fixed_psi_undefined_warns_and_returns_empty(
    mode_name,
    h,
    k,
    l,  # noqa: E741
    context,
):
    """When ψ is undefined for the reflection (Q ∥ reference), ``forward()``
    warns and returns ``[]``.

    The warning message must mention the azimuth so the
    user can choose a different reference direction if desired.
    """
    with context:
        g = _setup_psic_cubic()
        g.mode_name = mode_name
        assert natural_psi(g, h, k, l) is None
        with pytest.warns(UserWarning, match=re.escape("undefined")):
            solutions = g.forward(h, k, l)
        assert solutions == []


# ---------------------------------------------------------------------------
# Matching-target path returns solutions and does NOT warn.
# ---------------------------------------------------------------------------


def test_fixed_psi_horizontal_with_natural_target_returns_solutions():
    """Setting the constraint ψ to the value returned by ``natural_psi``
    makes the reflection reachable and produces real Bragg solutions.

    This is the recommended user workflow: call ``natural_psi`` to
    discover the achievable ψ, set the constraint, then call
    ``forward()``.  The forward solver must not warn in this case.
    """
    g = _setup_psic_cubic()
    g.mode_name = "fixed_psi_horizontal"
    nat = natural_psi(g, 1, 1, 0)
    assert nat is not None
    _set_psi_constraint(g.modes["fixed_psi_horizontal"], nat)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)  # any UserWarning fails
        solutions = g.forward(1, 1, 0)
    assert len(solutions) >= 1
    # Every solution must report the same ψ as natural_psi (the central
    # claim of issue #176).
    for sol in solutions:
        assert psi_angle(g, angles=sol) == pytest.approx(nat, abs=1e-6)


def test_fixed_psi_vertical_with_natural_target_returns_solutions():
    """Mirror of the horizontal test for the vertical mode."""
    g = _setup_psic_cubic()
    g.mode_name = "fixed_psi_vertical"
    nat = natural_psi(g, 1, 1, 0)
    assert nat is not None
    _set_psi_constraint(g.modes["fixed_psi_vertical"], nat)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        solutions = g.forward(1, 1, 0)
    assert len(solutions) >= 1
    for sol in solutions:
        assert psi_angle(g, angles=sol) == pytest.approx(nat, abs=1e-6)


# ---------------------------------------------------------------------------
# Issue #278's exact reproducer: sapphire + horizontal mode.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "h, k, l, context",
    [
        # Issue #278 lists these as "forward() returns 0 solutions" —
        # the user's own Newton check confirmed they are not reachable
        # at ψ = 0.  These tests pin that the empty return is now
        # accompanied by an informative warning.
        pytest.param(0, 0, 3, does_not_raise(), id="003"),
        pytest.param(0, 0, 6, does_not_raise(), id="006"),
        pytest.param(1, 1, 0, does_not_raise(), id="110"),
        pytest.param(0, 1, 2, does_not_raise(), id="012"),
        pytest.param(1, 0, 4, does_not_raise(), id="104"),
    ],
)
def test_issue_278_sapphire_reproducer_warns(h, k, l, context):  # noqa: E741
    """The exact sapphire reproducer from issue #278: every hkl yields a
    warning + empty list (correct behavior under the validation-filter
    model).
    """
    with context:
        g = psic()
        g.wavelength = 1.0
        g.sample.lattice = ahd.Lattice(
            a=4.758, b=4.758, c=12.991, alpha=90, beta=90, gamma=120
        )
        ub_identity(g.sample)
        g.azimuth = (0, 0, 1)
        g.surface_normal = (0, 0, 1)
        g.mode_name = "fixed_psi_horizontal"

        with pytest.warns(UserWarning):
            solutions = g.forward(h, k, l)
        assert solutions == []
