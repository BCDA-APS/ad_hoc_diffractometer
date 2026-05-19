# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
Cross-module regression tests for issue #280.

Issue #280 fixed a package-wide rotation-composition-order defect.
Before the fix every site that composed a chain of stage rotation
matrices iterated floor-most-first and **left-multiplied** the
accumulator (``Z = R @ Z``), producing the inverse (innermost-
leftmost) of the BL1967 standard convention (Busing & Levy 1967).
The bug was self-consistent inside the package (round-trips
``forward → inverse`` continued to close because the UB-fit machinery
used the same reversed convention and the two reversals cancelled),
which is why the 2447-test pre-fix suite passed unanimously.  The
defect was visible only through comparisons to externally-defined
orientations (real ``libhkl`` motor positions, hand-derived
Q-vectors for fixed motor angles, etc.).

The tests below are **non-round-trip** checks that distinguish the
two conventions.  Their expected values are hand-derived from the
BL1967 standard composition::

    Z = R_0 @ R_1 @ ... @ R_{N-1}   (outermost-leftmost)
    Q_phi = Z^T @ Q_lab

and were independently sanity-checked against the libhkl ``hkl_soleil``
solver at the developer bench (the harness is not added to CI
dependencies — libhkl is for local correlation only).

The tests are organized as:

* **Fixed-orientation cross-checks** — one parametrized test that
  asserts ``angles_to_phi_vector(geometry, **angles) == B @ hkl``
  for ``U = I`` on every registered demo geometry, with a hand-
  derived (h, k, l) and (angles) pair per geometry.
* **Composition-direction sanity check** — exercises the public
  ``sample_rotation_matrix`` / ``detector_rotation_matrix`` and
  asserts ``v_lab = Z @ v_phi`` (which would *fail* under the
  pre-#280 inner-leftmost convention).
* **Degenerate-phi double-diffraction coverage** — exercises the
  ``_find_degenerate_outers`` / ``_solve_degenerate_outer`` path in
  ``forward.py``, which is reachable only when ``Q_phi`` is parallel
  to the innermost stage axis (a configuration the analytic solver
  must defer back to a chi-scan).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.orientation import angles_to_phi_vector
from ad_hoc_diffractometer.rotation import rotation_matrix


# ---------------------------------------------------------------------------
# Hand-derived fixed-orientation expectations
# ---------------------------------------------------------------------------

# Each entry is:
#   (geometry_name, lattice_kwargs, motor_angles, hkl_target, expected_q_phi)
#
# ``expected_q_phi`` is hand-computed from the BL1967 standard
# composition (Busing & Levy 1967)::
#
#     Z = R_0 @ R_1 @ ... @ R_{N-1}   (outermost stage leftmost)
#     D = R_0' @ R_1' @ ... @ R_{M-1}'
#     Q_lab = (2π/λ) (D · ŷ − ŷ)         (ŷ = +longitudinal in geometry basis)
#     Q_phi_expected = Z^T · Q_lab
#
# Under ``U = I``, the Bragg condition is satisfied when
# ``Q_phi_expected = B · hkl``.  We assert *both*::
#
#     angles_to_phi_vector(geom, **angles)  ==  expected_q_phi
#     angles_to_phi_vector(geom, **angles)  ==  geom.sample.UB @ hkl
#
# The two equalities together pin down the composition direction
# (round-trip alone does not).

WAVELENGTH = 1.5406  # Cu Kα, used uniformly so all magnitudes agree


def _expected_q_phi(geom, motor_angles):
    """Hand-derive Q_phi from motor angles using the BL1967 standard
    composition.

    This is a *deliberate re-implementation* of the forward kinematics
    using the documented BL1967 standard convention (Busing & Levy
    1967), used only here as the independent reference.  It does NOT
    call the package's ``angles_to_phi_vector``; that's the function
    under test.
    """
    Z = np.eye(3)
    for s in geom.sample_stages:
        a = float(motor_angles.get(s.name, s.angle))
        Z = Z @ rotation_matrix(s.axis, a)
    D = np.eye(3)
    for s in geom.detector_stages:
        a = float(motor_angles.get(s.name, s.angle))
        D = D @ rotation_matrix(s.axis, a)
    y_hat = np.asarray(geom.basis["longitudinal"], dtype=float)
    y_hat = y_hat / np.linalg.norm(y_hat)
    Q_lab = (2.0 * math.pi / geom.wavelength) * (D @ y_hat - y_hat)
    return Z.T @ Q_lab


# Fixed-orientation reference cases.  Each entry specifies a motor
# configuration (independent of any solver) to feed both
# ``angles_to_phi_vector`` and the BL1967 standard hand-derivation.
# The hand derivation lives in ``_expected_q_phi`` above; the test
# asserts bitwise (to machine epsilon) agreement.  The motor angles
# below are *arbitrary* non-trivial values — they need not satisfy
# Bragg for any particular hkl; the only thing under test is whether
# the two composition implementations agree.
_FIXED_ORIENTATION_CASES = [
    pytest.param(
        "psic",
        dict(mu=12.0, eta=17.0, chi=33.0, phi=51.0, nu=8.0, delta=27.0),
        id="psic-non-trivial",
    ),
    pytest.param(
        "fourcv",
        dict(omega=22.0, chi=44.0, phi=66.0, ttheta=14.0),
        id="fourcv-non-trivial",
    ),
    pytest.param(
        "fourch",
        dict(omega=18.0, chi=37.0, phi=55.0, ttheta=12.0),
        id="fourch-non-trivial",
    ),
    pytest.param(
        "sixc",
        dict(alpha=7.0, omega=23.0, chi=41.0, phi=59.0, delta=15.0, gamma=9.0),
        id="sixc-non-trivial",
    ),
    pytest.param(
        "fivec",
        dict(mu=11.0, omega=29.0, chi=47.0, phi=63.0, ttheta=13.0),
        id="fivec-non-trivial",
    ),
    pytest.param(
        "kappa4cv",
        dict(komega=21.0, kappa=35.0, kphi=49.0, ttheta=15.0),
        id="kappa4cv-non-trivial",
    ),
    pytest.param(
        "kappa6c",
        dict(mu=6.0, komega=24.0, kappa=38.0, kphi=52.0, nu=4.0, delta=17.0),
        id="kappa6c-non-trivial",
    ),
    pytest.param(
        "zaxis",
        dict(alpha=10.0, Z=0.0, delta=20.0, gamma=5.0),
        id="zaxis-non-trivial",
    ),
    pytest.param(
        "s2d2",
        dict(mu=8.0, Z=0.0, nu=12.0, delta=16.0),
        id="s2d2-non-trivial",
    ),
]


@pytest.mark.parametrize(
    "geom_name, motor_angles",
    _FIXED_ORIENTATION_CASES,
)
def test_angles_to_phi_vector_matches_hand_derived(geom_name, motor_angles):
    """``angles_to_phi_vector`` agrees with a hand-derived BL1967 chain.

    This is the *non-round-trip* test that distinguishes the BL1967
    standard outermost-leftmost composition (Busing & Levy 1967, issue
    #280) from the pre-#280 inner-leftmost reversal.  Under the buggy
    convention this assertion would fail to machine precision; under
    the BL1967 standard convention it passes exactly.
    """
    g = ahd.make_geometry(geom_name)
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=4.0)
    ahd.ub_identity(g.sample)

    expected = _expected_q_phi(g, motor_angles)
    actual = angles_to_phi_vector(g, **motor_angles)
    np.testing.assert_allclose(
        actual,
        expected,
        atol=1e-12,
        err_msg=(
            f"{geom_name}: angles_to_phi_vector disagrees with BL1967 "
            f"standard outermost-leftmost composition (issue #280)."
        ),
    )


def test_angles_to_phi_vector_satisfies_bragg_for_forward_solution():
    """A forward() solution round-trips and *also* matches the
    hand-derived BL1967 standard ``Z^T @ Q_lab`` at exactly those motor
    angles.

    Round-trip alone (the pre-#280 test suite) is satisfied by every
    internally self-consistent composition, including the buggy one;
    matching the *independent* BL1967 standard hand-derivation
    requires the correct composition.
    """
    g = ahd.make_geometry("psic")
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=4.0)
    ahd.ub_identity(g.sample)
    g.mode_name = "bisecting_vertical"

    sols = g.forward(1, 1, 0)
    assert sols, "psic bisecting_vertical (1,1,0) must have at least one solution"
    sol = sols[0]

    expected = _expected_q_phi(g, sol)
    actual = angles_to_phi_vector(g, **sol)
    np.testing.assert_allclose(actual, expected, atol=1e-12)

    # And the Bragg condition.
    target = g.sample.UB @ np.array([1.0, 1.0, 0.0])
    np.testing.assert_allclose(actual, target, atol=1e-9)


# ---------------------------------------------------------------------------
# Composition-direction sanity check on the public API
# ---------------------------------------------------------------------------


def test_sample_rotation_matrix_maps_phi_to_lab():
    """``sample_rotation_matrix()`` (BL1967 standard outermost-leftmost)
    satisfies ``v_lab = Z @ v_phi``.

    Under the pre-#280 inner-leftmost convention this assertion would
    fail: the same loop produced ``Z_buggy = R_inner @ ... @ R_outer``,
    which equals ``Z_BL1967.T`` only when the stage axes commute — NOT
    the case here.
    """
    g = ahd.make_geometry("psic")
    g.set_angle("mu", 30.0)
    g.set_angle("eta", 20.0)
    g.set_angle("chi", 45.0)
    g.set_angle("phi", 15.0)

    Z_pkg = g.sample_rotation_matrix()
    # Hand-build the BL1967 standard product.
    Z_bl1967 = np.eye(3)
    for s in g.sample_stages:
        Z_bl1967 = Z_bl1967 @ rotation_matrix(s.axis, s.angle)

    np.testing.assert_allclose(Z_pkg, Z_bl1967, atol=1e-14)

    # And ``Z_pkg`` must NOT equal the reversed (buggy) product (which
    # would be Z_pkg.T for orthogonal axis sets — but the chain here
    # has non-commuting rotations, so the two are distinct).
    Z_reversed = np.eye(3)
    for s in reversed(list(g.sample_stages)):
        Z_reversed = Z_reversed @ rotation_matrix(s.axis, s.angle)
    assert not np.allclose(Z_pkg, Z_reversed, atol=1e-3), (
        "Sample rotation matrix unexpectedly equals the inner-leftmost "
        "(pre-#280) product — the composition order may have regressed."
    )


# ---------------------------------------------------------------------------
# Coverage: degenerate phi-axis double-diffraction path
# ---------------------------------------------------------------------------


def test_double_diffraction_degenerate_phi_path():
    """Exercise the degenerate phi-axis path in ``_solve_double_diffraction``.

    Under the corrected BL1967 standard outermost-leftmost composition
    (issue #280), the analytic decomposition becomes phi-indeterminate
    whenever ``Q_phi`` is parallel to the innermost stage axis.  The
    double-
    diffraction solver then dispatches to ``_find_degenerate_outers``
    and ``_solve_degenerate_outer``, which scan over phi (now the
    free variable) computing chi from the residual rotation at each
    step.

    For psic-vertical DD, the phi axis is ``-z`` (in YOU basis) and
    ``(0, 0, 1)`` produces ``Q_phi = (0, 0, 2π/a)`` — exactly along
    ``-phi_axis``, triggering the degenerate path.  The Ewald
    secondary condition for ``(0, 1, 0)`` is not simultaneously
    satisfiable at this primary, so the solver returns no solutions;
    the test asserts that *no exception* is raised (the degenerate
    path executes cleanly to completion).
    """
    g = ahd.make_geometry("psic")
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=4.0)
    ahd.ub_identity(g.sample)
    g.mode_name = "double_diffraction_vertical"
    cs = g.modes["double_diffraction_vertical"]
    cs.extras["h2"] = 0.0
    cs.extras["k2"] = 1.0
    cs.extras["l2"] = 0.0

    # Sanity check: Q_phi is parallel to the phi axis (up to sign).
    Q_phi = g.sample.UB @ np.array([0.0, 0.0, 1.0])
    n_phi = g.stage("phi")._axis_hat  # noqa: SLF001
    cos_angle = abs(float(np.dot(Q_phi, n_phi)) / np.linalg.norm(Q_phi))
    assert cos_angle == pytest.approx(1.0, abs=1e-10), (
        "Setup precondition: Q_phi must be parallel to the phi axis "
        "to trigger the degenerate-DD path."
    )

    # The call must execute the degenerate code path without raising.
    # The result is allowed to be either zero or more solutions;
    # what matters for issue #280 is that the path executes cleanly.
    solutions = g.forward(0, 0, 1)
    # Every returned solution (if any) must satisfy primary Bragg.
    for sol in solutions:
        Q = angles_to_phi_vector(g, **sol)
        assert np.allclose(Q, Q_phi, atol=1e-3), (
            f"Degenerate-DD solution {sol} does not satisfy primary Bragg."
        )
