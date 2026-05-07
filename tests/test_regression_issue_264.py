# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
Cross-module regression tests for issue #264.

Issue #264 introduced two distinct sets of psic-mode changes:

1. The SPEC ``OMEGA`` pseudo-angle (``def OMEGA 'Q[6]'`` —
   "the angle between Q and the plane of the chi circle") was added as
   a new :class:`~ad_hoc_diffractometer.mode.ReferenceConstraint` name,
   together with two new psic modes
   (``fixed_omega_vertical`` / ``fixed_omega_horizontal``) and the
   solver dispatch that handles them.

2. Three additional psic modes were introduced or revised per the
   @jwkim-anl review: ``fixed_alpha_i_fixed_chi_fixed_phi`` (B3),
   ``lifting_detector_eta`` (B4), plus the
   ``lifting_detector_phi`` / ``lifting_detector_mu`` revisions
   (C3/C4) that drop the ``qaz = 90`` constraint and fix every sample
   stage except the named one.

The per-module unit tests in ``tests/test_reference.py``,
``tests/test_forward.py`` and ``tests/test_regression_issue_267.py``
already cover the individual pieces; this file collects the
cross-module *invariants* that must hold across the whole #264 patch:

- OMEGA = 0 ⇔ bisecting (the central physical equivalence claimed by
  @jwkim-anl in the issue thread).
- Every mode named by issue #264 is present in the registry, has the
  expected ``is_implemented`` status, and produces solutions that
  satisfy the Bragg condition end-to-end.
- The dispatcher routes each new/revised mode to its intended solver
  branch (``_solve_omega_mode`` for ``fixed_omega_*``,
  ``_solve_free_detectors`` for ``fixed_alpha_i_fixed_chi_fixed_phi``,
  ``lifting_detector_eta``, and the revised
  ``lifting_detector_phi`` / ``lifting_detector_mu``).

C1 (``fixed_psi_vertical``) and C2 (``fixed_psi_horizontal``) drop
the previous BisectConstraint and pin a detector stage instead per
the @jwkim-anl review on the issue thread:

- ``fixed_psi_vertical``   = ``nu = 0`` (vertical) + ``mu`` fixed +
  ``psi`` target.
- ``fixed_psi_horizontal`` = ``delta = 0`` (horizontal) + ``eta``
  fixed + ``psi`` target.

Both modes route through ``_solve_psi_mode`` for the ψ-validation
filter, then to ``_solve_fixed_sample`` to solve the remaining free
sample stages plus the active detector at ``2θ``.
"""

from __future__ import annotations

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from helpers import psic

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import REQUIRED
from ad_hoc_diffractometer import ConstraintSet
from ad_hoc_diffractometer import DetectorConstraint
from ad_hoc_diffractometer import ReferenceConstraint
from ad_hoc_diffractometer import SampleConstraint
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer.forward import _is_free_detectors_mode
from ad_hoc_diffractometer.forward import _is_omega_mode
from ad_hoc_diffractometer.reference import incidence_angle
from ad_hoc_diffractometer.reference import omega_pseudo

WAVELENGTH = 1.5406  # Cu Kα

# All modes added or revised by #264.
_ISSUE_264_NEW_MODES = {
    "fixed_omega_vertical",
    "fixed_omega_horizontal",
    "fixed_alpha_i_fixed_chi_fixed_phi",
    "lifting_detector_eta",
}
_ISSUE_264_REVISED_MODES = {
    "lifting_detector_phi",
    "lifting_detector_mu",
    "fixed_psi_vertical",
    "fixed_psi_horizontal",
}


def _setup_psic_cubic(a: float = 4.0):
    """Return a psic geometry with UB=B for a cubic lattice of side ``a``."""
    g = psic()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=a)
    ub_identity(g.sample)
    return g


# ---------------------------------------------------------------------------
# Mode registry: every #264 mode is present
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name",
    sorted(_ISSUE_264_NEW_MODES | _ISSUE_264_REVISED_MODES),
)
def test_issue_264_mode_present(mode_name):
    """Every #264 mode is registered in the psic mode dict."""
    g = psic()
    assert mode_name in g.modes


# ---------------------------------------------------------------------------
# OMEGA = 0 ⇔ bisecting (the @jwkim-anl equivalence)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "h, k, l, context",
    [
        pytest.param(1, 0, 0, does_not_raise(), id="100"),
        pytest.param(0, 1, 1, does_not_raise(), id="011"),
        pytest.param(1, 1, 1, does_not_raise(), id="111"),
    ],
)
def test_omega_zero_equals_bisecting_vertical(h, k, l, context):  # noqa: E741
    """OMEGA = 0 ⇒ bisecting_vertical (vertical scattering plane).

    @jwkim-anl wrote on issue #264:
        "Yes. This is including bisecting mode. If omega is fixed at 0,
        it is bisecting."

    Verifies the equivalence numerically: the solution sets returned
    by ``bisecting_vertical`` and ``fixed_omega_vertical`` (target 0)
    must agree motor-for-motor for every reachable reflection.
    """
    with context:
        g = _setup_psic_cubic()

        g.mode_name = "bisecting_vertical"
        bisect_sols = g.forward(h, k, l)
        g.mode_name = "fixed_omega_vertical"
        omega_sols = g.forward(h, k, l)

        assert len(bisect_sols) == len(omega_sols), (
            f"({h},{k},{l}): bisecting returned {len(bisect_sols)} sols, "
            f"omega=0 returned {len(omega_sols)}"
        )

        bisect_sorted = sorted(bisect_sols, key=lambda s: (s["eta"], s["chi"]))
        omega_sorted = sorted(omega_sols, key=lambda s: (s["eta"], s["chi"]))
        for b, o in zip(bisect_sorted, omega_sorted, strict=False):
            for stage in ("mu", "eta", "chi", "phi", "nu", "delta"):
                assert b[stage] == pytest.approx(o[stage], abs=1e-6), (
                    f"({h},{k},{l}) stage {stage}: "
                    f"bisecting={b[stage]}, omega=0={o[stage]}"
                )
            # Independent confirmation: omega_pseudo evaluates to 0 in
            # both solution sets.
            assert omega_pseudo(g, angles=b) == pytest.approx(0.0, abs=1e-7)
            assert omega_pseudo(g, angles=o) == pytest.approx(0.0, abs=1e-7)


@pytest.mark.parametrize(
    "h, k, l, context",
    [
        pytest.param(0, 0, 1, does_not_raise(), id="001"),
        pytest.param(1, 0, 1, does_not_raise(), id="101"),
    ],
)
def test_omega_zero_equals_bisecting_horizontal(h, k, l, context):  # noqa: E741
    """OMEGA = 0 ⇒ bisecting_horizontal (horizontal scattering plane)."""
    with context:
        g = _setup_psic_cubic()

        g.mode_name = "bisecting_horizontal"
        bisect_sols = g.forward(h, k, l)
        g.mode_name = "fixed_omega_horizontal"
        omega_sols = g.forward(h, k, l)

        assert len(bisect_sols) == len(omega_sols)
        bisect_sorted = sorted(bisect_sols, key=lambda s: s["mu"])
        omega_sorted = sorted(omega_sols, key=lambda s: s["mu"])
        for b, o in zip(bisect_sorted, omega_sorted, strict=False):
            for stage in ("mu", "eta", "chi", "phi", "nu", "delta"):
                assert b[stage] == pytest.approx(o[stage], abs=1e-6)
            assert omega_pseudo(g, angles=b) == pytest.approx(0.0, abs=1e-7)
            assert omega_pseudo(g, angles=o) == pytest.approx(0.0, abs=1e-7)


# ---------------------------------------------------------------------------
# OMEGA pseudo-angle is the SPEC Q[6] definition
# ---------------------------------------------------------------------------


def test_omega_is_angle_between_q_and_chi_circle_plane():
    """OMEGA matches the SPEC definition: angle between Q and the
    plane of the chi circle.

    The chi circle plane is the plane perpendicular to the chi
    rotation axis after the outer (mu, eta) sample stages have been
    applied.  This test computes both the
    :func:`~ad_hoc_diffractometer.reference.omega_pseudo` value and
    the geometric angle from first principles, then asserts they
    agree.
    """
    g = _setup_psic_cubic()

    angles = {
        "mu": 7.0,
        "eta": 13.0,
        "chi": 25.0,
        "phi": 41.0,
        "nu": 0.0,
        "delta": 30.0,
    }

    om = omega_pseudo(g, angles=angles)

    # First-principles computation of the chi-circle plane normal in
    # the lab frame:  apply mu and eta rotations to the chi axis vector.
    from ad_hoc_diffractometer.rotation import _rotation_matrix_normalized

    chi_stage = g.stage("chi")
    mu_stage = g.stage("mu")
    eta_stage = g.stage("eta")
    R_mu = _rotation_matrix_normalized(mu_stage._axis_hat, angles["mu"])  # noqa: SLF001
    R_eta = _rotation_matrix_normalized(eta_stage._axis_hat, angles["eta"])  # noqa: SLF001
    chi_axis_lab = R_eta @ R_mu @ chi_stage._axis_hat  # noqa: SLF001
    chi_axis_lab /= np.linalg.norm(chi_axis_lab)

    # Lab-frame Q vector.
    nu_stage = g.stage("nu")
    delta_stage = g.stage("delta")
    R_nu = _rotation_matrix_normalized(nu_stage._axis_hat, angles["nu"])  # noqa: SLF001
    R_delta = _rotation_matrix_normalized(delta_stage._axis_hat, angles["delta"])  # noqa: SLF001
    D = R_delta @ R_nu
    y_hat = np.asarray(g.basis["longitudinal"], dtype=float)
    y_hat /= np.linalg.norm(y_hat)
    Q_lab = D @ y_hat - y_hat
    Q_hat = Q_lab / np.linalg.norm(Q_lab)

    # OMEGA = arcsin( Q_hat · n_chi ).  Sign is positive when Q has
    # positive projection on the chi axis.
    sin_om = float(np.dot(Q_hat, chi_axis_lab))
    expected_om = np.degrees(np.arcsin(np.clip(sin_om, -1.0, 1.0)))
    assert om == pytest.approx(expected_om, abs=1e-9), (
        f"omega_pseudo() = {om}, expected {expected_om}"
    )


# ---------------------------------------------------------------------------
# Dispatch routing: every #264 mode lands on the intended solver
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, predicate",
    [
        pytest.param(
            "fixed_omega_vertical",
            _is_omega_mode,
            id="fixed_omega_vertical-omega_solver",
        ),
        pytest.param(
            "fixed_omega_horizontal",
            _is_omega_mode,
            id="fixed_omega_horizontal-omega_solver",
        ),
        pytest.param(
            "lifting_detector_eta",
            _is_free_detectors_mode,
            id="lifting_detector_eta-free_detectors",
        ),
        pytest.param(
            "lifting_detector_phi",
            _is_free_detectors_mode,
            id="lifting_detector_phi-free_detectors",
        ),
        pytest.param(
            "lifting_detector_mu",
            _is_free_detectors_mode,
            id="lifting_detector_mu-free_detectors",
        ),
    ],
)
def test_issue_264_mode_routes_to_intended_solver(mode_name, predicate):
    """Each #264 mode is detected by its intended dispatcher predicate.

    Catches dispatch regressions that would silently route a mode
    through a generic solver and break the constraint semantics.
    """
    g = _setup_psic_cubic()
    cs = g.modes[mode_name]
    assert predicate(g, cs) is True


def test_fixed_alpha_i_fixed_chi_fixed_phi_routes_to_free_detectors():
    """B3 routes to ``_solve_free_detectors`` only after ``surface_normal``
    is set (otherwise the mode is a stub)."""
    g = _setup_psic_cubic()
    cs = g.modes["fixed_alpha_i_fixed_chi_fixed_phi"]
    # Without surface_normal, the alpha_i ReferenceConstraint reports
    # is_implemented=False, so the mode is a stub regardless of dispatch.
    assert cs.is_implemented(g) is False
    g.surface_normal = (0, 0, 1)
    assert cs.is_implemented(g) is True
    assert _is_free_detectors_mode(g, cs) is True


# ---------------------------------------------------------------------------
# End-to-end sanity: every implemented #264 mode round-trips for at least
# one canonical reflection.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, h, k, l, needs_surface_normal",
    [
        pytest.param("fixed_omega_vertical", 1, 0, 0, False, id="fov-100"),
        pytest.param("fixed_omega_horizontal", 0, 0, 1, False, id="foh-001"),
        pytest.param("fixed_alpha_i_fixed_chi_fixed_phi", 0, 1, 1, True, id="b3-011"),
        pytest.param("lifting_detector_eta", 1, 1, 0, False, id="le-110"),
        pytest.param("lifting_detector_phi", 1, 0, 0, False, id="lp-100"),
        pytest.param("lifting_detector_mu", 0, 1, 0, False, id="lm-010"),
    ],
)
def test_issue_264_mode_round_trip(
    mode_name,
    h,
    k,
    l,  # noqa: E741
    needs_surface_normal,
):
    """Every #264 implemented mode produces solutions that round-trip
    back to the requested (h, k, l) under ``inverse()``."""
    g = _setup_psic_cubic()
    if needs_surface_normal:
        g.surface_normal = (0, 0, 1)
    g.mode_name = mode_name
    sols = g.forward(h, k, l)
    assert len(sols) > 0, f"{mode_name} ({h},{k},{l}): no solutions"
    for sol in sols:
        hkl_back = g.inverse(sol)
        assert np.allclose(hkl_back, [h, k, l], atol=1e-5), (
            f"{mode_name} ({h},{k},{l}): inverse mismatch {hkl_back}"
        )


# ---------------------------------------------------------------------------
# Issue #264 design invariant: dropping qaz from C3/C4 freed both detectors.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name",
    ["lifting_detector_phi", "lifting_detector_mu"],
)
def test_revised_lifting_detector_has_no_qaz_constraint(mode_name):
    """The C3/C4 revision dropped the qaz=90 ``DetectorConstraint``."""
    g = psic()
    cs = g.modes[mode_name]
    # No DetectorConstraint at all in the revised modes.
    detector_constraints = [
        c for c in cs.constraints if isinstance(c, DetectorConstraint)
    ]
    assert detector_constraints == [], (
        f"{mode_name}: expected no DetectorConstraint after C3/C4 revision; "
        f"found {detector_constraints!r}"
    )
    # Three SampleConstraints — every sample stage except the free one.
    sample_constraints = [c for c in cs.constraints if isinstance(c, SampleConstraint)]
    assert len(sample_constraints) == 3, (
        f"{mode_name}: expected 3 SampleConstraints; got {len(sample_constraints)}"
    )


# ---------------------------------------------------------------------------
# Issue #264 design invariant: revised fixed_psi_* modes drop the bisect
# and pin a detector stage instead (C1/C2).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, expected_detector_stage, expected_sample_stage",
    [
        pytest.param("fixed_psi_vertical", "nu", "mu", id="fixed_psi_vertical-nu-mu"),
        pytest.param(
            "fixed_psi_horizontal",
            "delta",
            "eta",
            id="fixed_psi_horizontal-delta-eta",
        ),
    ],
)
def test_revised_fixed_psi_constraint_shape(
    mode_name, expected_detector_stage, expected_sample_stage
):
    """The C1/C2 revision swapped the BisectConstraint for a
    DetectorConstraint pinning the scattering plane.
    """
    from ad_hoc_diffractometer.mode import BisectConstraint

    g = psic()
    cs = g.modes[mode_name]
    # No BisectConstraint after the revision.
    bisects = [c for c in cs.constraints if isinstance(c, BisectConstraint)]
    assert bisects == [], (
        f"{mode_name}: expected no BisectConstraint after C1/C2 revision; "
        f"found {bisects!r}"
    )
    # Exactly one DetectorConstraint pinning the scattering plane.
    det = cs.detector_constraint
    assert det is not None
    assert det.name == expected_detector_stage
    assert det.value == 0.0
    # Exactly one SampleConstraint at the named stage (default value 0).
    samples = [c for c in cs.constraints if isinstance(c, SampleConstraint)]
    assert len(samples) == 1
    assert samples[0].name == expected_sample_stage
    # The psi reference is still present.
    ref = cs.reference_constraint
    assert ref is not None
    assert ref.name == "psi"


def _natural_psi(g, h, k, l):  # noqa: E741
    """Return the natural psi for (h, k, l) on geometry ``g``."""
    from ad_hoc_diffractometer.forward import _compute_natural_psi

    Q_phi = g.sample.UB @ np.array([h, k, l], float)
    return _compute_natural_psi(g, Q_phi)


@pytest.mark.parametrize(
    "mode_name, h, k, l",
    [
        pytest.param("fixed_psi_vertical", 1, 0, 0, id="fpv-100"),
        pytest.param("fixed_psi_vertical", 1, 1, 0, id="fpv-110"),
        pytest.param("fixed_psi_vertical", 1, 1, 1, id="fpv-111"),
        pytest.param("fixed_psi_horizontal", 1, 0, 1, id="fph-101"),
        pytest.param("fixed_psi_horizontal", 1, 0, 0, id="fph-100"),
    ],
)
def test_revised_fixed_psi_round_trip(
    mode_name,
    h,
    k,
    l,  # noqa: E741
):
    """C1/C2 revised fixed_psi modes round-trip and satisfy psi target."""
    from ad_hoc_diffractometer.reference import psi_angle

    g = _setup_psic_cubic()
    g.azimuthal_reference = (0, 0, 1)
    natural = _natural_psi(g, h, k, l)
    assert natural is not None, (
        f"{mode_name} ({h},{k},{l}): natural psi is undefined for this hkl"
    )

    # Patch the mode to use the natural psi target so the validator passes.
    old_mode = g.modes[mode_name]
    new_constraints = []
    for c in old_mode.constraints:
        if isinstance(c, ReferenceConstraint) and c.name == "psi":
            new_constraints.append(ReferenceConstraint("psi", natural))
        else:
            new_constraints.append(c)
    g.modes[mode_name] = ConstraintSet(
        new_constraints,
        computed=old_mode.computed,
        extras=dict(old_mode.extras),
    )
    g.mode_name = mode_name

    sols = g.forward(h, k, l)
    assert len(sols) > 0, f"{mode_name} ({h},{k},{l}): no solutions"
    for sol in sols:
        # psi target satisfied
        psi = psi_angle(g, angles=sol)
        assert psi == pytest.approx(natural, abs=1e-3), (
            f"{mode_name} ({h},{k},{l}): expected psi={natural}, got {psi}"
        )
        # Bragg round-trip
        hkl_back = g.inverse(sol)
        assert np.allclose(hkl_back, [h, k, l], atol=1e-5)


def test_revised_fixed_psi_wrong_target_returns_empty():
    """C1/C2: psi-validation rejects requests whose natural psi differs
    from the stored target."""
    g = _setup_psic_cubic()
    g.azimuthal_reference = (0, 0, 1)
    natural = _natural_psi(g, 1, 0, 0)
    assert natural is not None

    g.modes["fixed_psi_vertical"] = ConstraintSet(
        [
            DetectorConstraint("nu", 0.0),
            SampleConstraint("mu", 0.0),
            ReferenceConstraint("psi", natural + 45.0),
        ],
        computed=["eta", "chi", "phi", "delta"],
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = "fixed_psi_vertical"
    assert g.forward(1, 0, 0) == []


# ---------------------------------------------------------------------------
# Issue #264 design invariant: OMEGA pseudo-angle constraint accepts only
# the new "omega" reference name and is gated on the chi sample stage.
# ---------------------------------------------------------------------------


def test_omega_reference_constraint_validation():
    """``ReferenceConstraint('omega', value)`` is now valid; nonsense
    names still raise."""
    rc = ReferenceConstraint("omega", 5.0)
    assert rc.name == "omega"
    assert rc.value == 5.0
    with pytest.raises(ValueError, match=re.escape("ReferenceConstraint name")):
        ReferenceConstraint("not_a_pseudo_angle", 0.0)


def test_omega_constraint_unimplemented_without_chi_stage():
    """``ReferenceConstraint('omega', ...)`` is implemented only on
    geometries with a ``chi`` sample stage."""
    # Build a synthetic geometry without a chi stage (zaxis has none).
    g = ahd.make_geometry("zaxis")
    g.wavelength = WAVELENGTH
    rc = ReferenceConstraint("omega", 0.0)
    assert rc.has_reference_vector(g) is True  # no reference needed
    assert rc.is_implemented(g) is False  # but the geometry can't run it


# ---------------------------------------------------------------------------
# Issue #264 design invariant: B3 produces solutions whose alpha_i matches
# the requested target.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "alpha_target",
    [0.0, 3.0, 7.5],
)
def test_b3_alpha_i_target_satisfied(alpha_target):
    """B3 ``fixed_alpha_i_fixed_chi_fixed_phi`` solutions satisfy
    alpha_i == target within tolerance."""
    g = _setup_psic_cubic()
    g.surface_normal = (0, 0, 1)
    cs = g.modes["fixed_alpha_i_fixed_chi_fixed_phi"]
    # Override alpha_i target with the parametrized value.
    g.modes["__b3_test"] = ConstraintSet(
        [
            SampleConstraint("chi", 0.0),
            SampleConstraint("phi", 0.0),
            ReferenceConstraint("alpha_i", alpha_target),
        ],
        computed=cs.computed,
        extras=dict(cs.extras),
    )
    g.mode_name = "__b3_test"
    sols = g.forward(0, 1, 1)
    assert len(sols) > 0
    for sol in sols:
        ai = incidence_angle(g, angles=sol)
        assert ai == pytest.approx(alpha_target, abs=1e-3), (
            f"B3 alpha_i target {alpha_target}: got {ai}"
        )
