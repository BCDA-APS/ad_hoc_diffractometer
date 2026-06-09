# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
Regression tests for issue #292.

Issue #292 noted that the ``fixed_alpha_i_*``, ``fixed_beta_out_*``,
``alpha_eq_beta_*`` modes (psic, sixc, zaxis, s2d2) — and by extension
every mode declaring an ``alpha_i``, ``beta_out``, ``psi``, or ``omega``
output slot in ``mode.extras`` — left those slots at their YAML default
of ``None`` even after a successful ``forward()`` call.  This was a
latent contract bug: the slots advertised an output the solver never
delivered.

The fix adds ``_populate_output_extras`` to ``compute_forward``: after
the solver returns a non-empty list of solutions, each declared output
slot is replaced by a Python list of per-solution float values computed
via the corresponding helper in :mod:`ad_hoc_diffractometer.reference`:

* ``alpha_i``  → :func:`~ad_hoc_diffractometer.reference.incidence_angle`
* ``beta_out`` → :func:`~ad_hoc_diffractometer.reference.exit_angle`
* ``psi``      → :func:`~ad_hoc_diffractometer.reference.psi_angle`
* ``omega``    → :func:`~ad_hoc_diffractometer.reference.omega_pseudo`

Empty solution lists reset the slot to ``[]``; a missing reference
vector leaves the slot at ``None`` (with a debug-level log message).
Repeated ``forward()`` calls overwrite the slot.

The tests below pin the populated values against an independent
recomputation of the same helper per solution.
"""

from __future__ import annotations

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.reference import exit_angle
from ad_hoc_diffractometer.reference import incidence_angle
from ad_hoc_diffractometer.reference import omega_pseudo
from ad_hoc_diffractometer.reference import psi_angle

WAVELENGTH = 1.5406  # Cu Kα


def _setup_cubic(name, a=4.0):
    """Return a cubic geometry ready for forward()."""
    g = ahd.make_geometry(name)
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=a)
    ahd.ub_identity(g.sample)
    return g


# ---------------------------------------------------------------------------
# B3 surface mode: psic / fixed_alpha_i_fixed_chi_fixed_phi
# (the only fixed_alpha_i_* mode that is currently implemented for psic).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "h, k, l, alpha_target, context",
    [
        pytest.param(0, 1, 1, 0.0, does_not_raise(), id="011-ai0"),
        pytest.param(0, 1, 1, 5.0, does_not_raise(), id="011-ai5"),
        pytest.param(1, 1, 1, 3.0, does_not_raise(), id="111-ai3"),
    ],
)
def test_psic_b3_populates_alpha_i_and_beta_out_extras(
    h,
    k,
    l,  # noqa: E741
    alpha_target,
    context,
):
    """B3 mode populates ``extras['alpha_i']`` and ``extras['beta_out']``."""
    with context:
        g = _setup_cubic("psic")
        g.surface_normal = (0, 0, 1)
        cs = g.modes["fixed_alpha_i_fixed_chi_fixed_phi"]
        # Update the reference-constraint target via the public API:
        # rebuild the constraint set with the requested alpha_i value
        # (the YAML defaults to 0.0; we want non-zero targets too).
        from ad_hoc_diffractometer import REQUIRED
        from ad_hoc_diffractometer import ConstraintSet
        from ad_hoc_diffractometer import ReferenceConstraint
        from ad_hoc_diffractometer import SampleConstraint

        g.modes["fixed_alpha_i_fixed_chi_fixed_phi"] = ConstraintSet(
            [
                SampleConstraint("chi", 0.0),
                SampleConstraint("phi", 0.0),
                ReferenceConstraint("alpha_i", alpha_target),
            ],
            computed=cs.computed,
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        )
        g.mode_name = "fixed_alpha_i_fixed_chi_fixed_phi"

        sols = g.forward(h, k, l)
        assert len(sols) > 0

        mode = g.modes["fixed_alpha_i_fixed_chi_fixed_phi"]
        assert isinstance(mode.extras["alpha_i"], list)
        assert isinstance(mode.extras["beta_out"], list)
        assert len(mode.extras["alpha_i"]) == len(sols)
        assert len(mode.extras["beta_out"]) == len(sols)

        # Each populated alpha_i value matches the per-solution recomputation
        # and (per the mode's constraint) equals the target.
        for ai_stored, bo_stored, sol in zip(
            mode.extras["alpha_i"], mode.extras["beta_out"], sols, strict=True
        ):
            assert ai_stored == pytest.approx(incidence_angle(g, angles=sol), abs=1e-8)
            assert bo_stored == pytest.approx(exit_angle(g, angles=sol), abs=1e-8)
            assert ai_stored == pytest.approx(alpha_target, abs=1e-3)


# ---------------------------------------------------------------------------
# omega slot: psic / fixed_omega_vertical and fixed_omega_horizontal
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, h, k, l, context",
    [
        pytest.param(
            "fixed_omega_vertical", 0, 1, 1, does_not_raise(), id="omega-vert-011"
        ),
        pytest.param(
            "fixed_omega_vertical", 1, 1, 1, does_not_raise(), id="omega-vert-111"
        ),
        pytest.param(
            "fixed_omega_horizontal",
            0,
            1,
            1,
            does_not_raise(),
            id="omega-horiz-011",
        ),
    ],
)
def test_psic_fixed_omega_populates_omega_extra(
    mode_name,
    h,
    k,
    l,  # noqa: E741
    context,
):
    """``fixed_omega_*`` populates ``extras['omega']`` per solution."""
    with context:
        g = _setup_cubic("psic")
        g.mode_name = mode_name
        sols = g.forward(h, k, l)
        assert len(sols) > 0

        mode = g.modes[mode_name]
        assert isinstance(mode.extras["omega"], list)
        assert len(mode.extras["omega"]) == len(sols)
        for stored, sol in zip(mode.extras["omega"], sols, strict=True):
            assert stored == pytest.approx(omega_pseudo(g, angles=sol), abs=1e-8)


# ---------------------------------------------------------------------------
# psi slot: fourcv / fixed_psi (the always-implemented psi validation mode).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "h, k, l, ref, context",
    [
        pytest.param(0, 1, 0, (0, 0, 1), does_not_raise(), id="fourcv-010-ref001"),
        pytest.param(1, 1, 1, (0, 0, 1), does_not_raise(), id="fourcv-111-ref001"),
    ],
)
def test_fourcv_fixed_psi_populates_psi_extra(
    h,
    k,
    l,  # noqa: E741
    ref,
    context,
):
    """``fixed_psi`` populates ``extras['psi']`` per solution."""
    with context:
        from ad_hoc_diffractometer import REQUIRED
        from ad_hoc_diffractometer import ConstraintSet
        from ad_hoc_diffractometer import ReferenceConstraint
        from ad_hoc_diffractometer.forward import _compute_natural_psi

        g = _setup_cubic("fourcv")
        g.azimuthal_reference = ref
        Q_phi = g.sample.UB @ np.array([h, k, l], dtype=float)
        natural = _compute_natural_psi(g, Q_phi)
        assert natural is not None

        g.modes["fixed_psi"] = ConstraintSet(
            [ReferenceConstraint("psi", natural)],
            computed=g.modes["fixed_psi"].computed,
            extras={"n_hat": REQUIRED, "psi": None},
        )
        g.mode_name = "fixed_psi"
        sols = g.forward(h, k, l)
        assert len(sols) > 0

        mode = g.modes["fixed_psi"]
        assert isinstance(mode.extras["psi"], list)
        assert len(mode.extras["psi"]) == len(sols)
        for stored, sol in zip(mode.extras["psi"], sols, strict=True):
            assert stored == pytest.approx(psi_angle(g, angles=sol), abs=1e-6)
            # And by the validation-filter property, every solution must
            # have the natural psi (modulo 360).
            assert (
                stored == pytest.approx(natural, abs=1e-3)
                or stored == pytest.approx(natural - 360.0, abs=1e-3)
                or stored == pytest.approx(natural + 360.0, abs=1e-3)
            )


# ---------------------------------------------------------------------------
# Empty-solutions case: slot is reset to [].
# ---------------------------------------------------------------------------


def test_extras_reset_when_forward_returns_no_solutions():
    """When ``forward()`` returns ``[]`` the output slot is reset to ``[]``."""
    from ad_hoc_diffractometer import REQUIRED
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import ReferenceConstraint
    from ad_hoc_diffractometer.forward import _compute_natural_psi

    g = _setup_cubic("fourcv")
    g.azimuthal_reference = (0, 0, 1)

    # Choose a psi target that differs from the natural psi by 45° so the
    # validation filter rejects every solution and returns [].
    Q_phi = g.sample.UB @ np.array([0.0, 1.0, 0.0])
    natural = _compute_natural_psi(g, Q_phi)
    assert natural is not None
    wrong_target = natural + 45.0

    g.modes["fixed_psi"] = ConstraintSet(
        [ReferenceConstraint("psi", wrong_target)],
        computed=g.modes["fixed_psi"].computed,
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = "fixed_psi"
    with pytest.warns(UserWarning):
        sols = g.forward(0, 1, 0)
    assert sols == []
    assert g.modes["fixed_psi"].extras["psi"] == []


# ---------------------------------------------------------------------------
# Overwrite: a second forward() call replaces the first's contents.
# ---------------------------------------------------------------------------


def test_extras_overwritten_on_subsequent_forward_calls():
    """Each successful ``forward()`` replaces the prior slot contents."""
    g = _setup_cubic("psic")
    g.mode_name = "fixed_omega_vertical"

    sols_1 = g.forward(0, 1, 1)
    assert len(sols_1) > 0
    stored_1 = list(g.modes["fixed_omega_vertical"].extras["omega"])

    sols_2 = g.forward(1, 1, 1)
    assert len(sols_2) > 0
    stored_2 = list(g.modes["fixed_omega_vertical"].extras["omega"])

    # Distinct hkls produce distinct OMEGA pseudoangle sets, so the second
    # call's stored values are not a prefix of the first.
    assert stored_2 != stored_1


# ---------------------------------------------------------------------------
# Modes without any output slot are untouched.
# ---------------------------------------------------------------------------


def test_modes_without_output_slots_are_untouched():
    """``bisecting_vertical`` declares no output-slot extras: nothing is added."""
    g = _setup_cubic("psic")
    g.mode_name = "bisecting_vertical"
    sols = g.forward(0, 1, 1)
    assert len(sols) > 0
    extras = g.modes["bisecting_vertical"].extras
    for key in ("alpha_i", "beta_out", "psi", "omega"):
        assert key not in extras, (
            f"populate hook must not add {key!r} to a mode whose YAML did not "
            f"declare it; bisecting_vertical extras: {sorted(extras)}"
        )


# ---------------------------------------------------------------------------
# Soft-fail: if the underlying reference helper raises (e.g. surface_normal
# unset at the last instant, or psi-undefined geometry), the slot is left as
# ``None`` and a debug-level message is logged.  Triggered here by monkey-
# patching the reference helper so the rest of the solver stays valid.
# ---------------------------------------------------------------------------


def test_populate_extras_soft_fails_when_reference_helper_raises(monkeypatch, caplog):
    """A raising reference helper leaves the slot ``None`` and logs at DEBUG."""
    import logging

    from ad_hoc_diffractometer import forward as forward_mod

    g = _setup_cubic("psic")
    g.mode_name = "fixed_omega_vertical"

    def _explode(*_args, **_kwargs):
        raise ValueError("synthetic failure for issue #292 coverage")

    monkeypatch.setattr(
        forward_mod, "_populate_output_extras", forward_mod._populate_output_extras
    )
    monkeypatch.setattr(
        "ad_hoc_diffractometer.reference.omega_pseudo",
        _explode,
    )

    with caplog.at_level(logging.DEBUG, logger="ad_hoc_diffractometer.forward"):
        sols = g.forward(0, 1, 1)

    assert len(sols) > 0, "solver should still return solutions"
    assert g.modes["fixed_omega_vertical"].extras["omega"] is None
    assert any(
        "leaving mode.extras['omega'] as None" in rec.getMessage()
        for rec in caplog.records
    ), "expected a debug log record naming the soft-failed slot"
