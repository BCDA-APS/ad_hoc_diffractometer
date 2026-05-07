# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
Cross-module regression tests for issue #237.

The B matrix returned by :func:`ad_hoc_diffractometer.lattice.b_matrix`
was built with a stray ``.T`` that placed the reciprocal lattice
vectors as the **rows** of B instead of the **columns** required by
the Busing & Levy (1967) eq. 3 convention.  The bug was silent for
cubic, tetragonal, and orthorhombic cells (whose reciprocal vectors
are mutually orthogonal Cartesian axes, so ``M`` and ``M.T`` are
equal) and silently wrong for hexagonal, trigonal, monoclinic, and
triclinic cells whenever a reflection mixed two or more reciprocal
vectors along the same Cartesian axis (e.g. ``(h, k, 0)``-type
reflections in a hexagonal cell).

The canonical reproducer in the issue is hexagonal sapphire
(α-Al₂O₃, a = 4.785 Å, c = 12.991 Å, γ = 120°) with reflection
(1, 0, 0).  The pre-fix B gave ``|B @ (1, 0, 0)| = 1.3131 Å⁻¹``
(the projection of b1 onto x̂), whereas the true magnitude is
``|b1| = 2π / (a · sin γ) ≈ 1.5162 Å⁻¹``.

Cross-module: the bug originated in ``lattice.b_matrix`` but the
visible symptom appeared downstream in ``orientation`` (UB matrix),
``forward`` (motor angles for a given hkl), and any independent
cross-validation against external solvers (hkl_soleil/hklpy2).  This
test file lives at the cross-module level rather than in
``tests/test_lattice.py`` so the regression case is preserved against
any future refactor that moves the responsibility for B-matrix
construction.
"""

from __future__ import annotations

import math

import numpy as np
from helpers import fourcv

import ad_hoc_diffractometer as ahd

# ---------------------------------------------------------------------------
# Sapphire crystal — canonical reproducer from the issue
# ---------------------------------------------------------------------------

# α-Al₂O₃ (hexagonal): a = 4.785 Å, c = 12.991 Å, γ = 120°.
# Reciprocal lattice constants computed from first principles:
#   |b1| = |b2| = 2π / (a · sin 120°) = 1.5162383 Å⁻¹
#   |b3|       = 2π /  c              = 0.4836568 Å⁻¹
SAPPHIRE_KWARGS = dict(a=4.785, c=12.991, gamma=120.0)
SAPPHIRE_B1_MAG_EXPECTED = 2.0 * np.pi / (4.785 * math.sin(math.radians(120.0)))
SAPPHIRE_B3_MAG_EXPECTED = 2.0 * np.pi / 12.991


def test_sapphire_b1_magnitude():
    """``|b1| == 2π / (a · sin γ)`` for hexagonal sapphire."""
    lat = ahd.Lattice(**SAPPHIRE_KWARGS)
    b1, b2, b3 = lat.reciprocal_lattice_vectors
    np.testing.assert_allclose(np.linalg.norm(b1), SAPPHIRE_B1_MAG_EXPECTED, rtol=1e-10)
    np.testing.assert_allclose(np.linalg.norm(b2), SAPPHIRE_B1_MAG_EXPECTED, rtol=1e-10)
    np.testing.assert_allclose(np.linalg.norm(b3), SAPPHIRE_B3_MAG_EXPECTED, rtol=1e-10)


def test_sapphire_B_at_100_equals_b1():
    """``B @ (1, 0, 0)`` returns ``b1`` (not ``b1[0] · x̂``).

    This is the direct reproducer of the bug.  Pre-fix, the stray ``.T`` in
    :func:`ad_hoc_diffractometer.lattice.b_matrix` made ``B @ (1, 0, 0)``
    return ``(b1·x̂, b2·x̂, b3·x̂) = (1.31310, 0, 0)``, with magnitude
    1.31310 Å⁻¹.  The correct value is ``b1 = (1.31310, 0.75812, 0)`` with
    magnitude 1.51624 Å⁻¹.
    """
    lat = ahd.Lattice(**SAPPHIRE_KWARGS)
    b1, _, _ = lat.reciprocal_lattice_vectors
    Q = lat.B @ np.array([1.0, 0.0, 0.0])

    np.testing.assert_allclose(Q, b1, atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(Q), SAPPHIRE_B1_MAG_EXPECTED, rtol=1e-10)


def test_sapphire_B_at_006_unchanged_by_fix():
    """``B @ (0, 0, 6)`` was correct *by accident* both pre- and post-fix.

    For hexagonal sapphire, ``b3 = (0, 0, |b3|)`` is the only reciprocal
    vector with a nonzero z-component, so both layouts of B map
    ``(0, 0, 6)`` to ``(0, 0, 6 · |b3|)``.  This is why the issue's
    cross-check against hkl_soleil agreed for (0, 0, 6) but disagreed
    for (1, 0, 0).
    """
    lat = ahd.Lattice(**SAPPHIRE_KWARGS)
    Q = lat.B @ np.array([0.0, 0.0, 6.0])
    np.testing.assert_allclose(
        np.linalg.norm(Q), 6.0 * SAPPHIRE_B3_MAG_EXPECTED, rtol=1e-10
    )
    # Direction is along z (within the BL1967 crystal Cartesian frame)
    np.testing.assert_allclose(Q[0], 0.0, atol=1e-12)
    np.testing.assert_allclose(Q[1], 0.0, atol=1e-12)


def test_sapphire_bragg_2theta_at_100():
    """The Bragg 2θ for sapphire (1, 0, 0) at λ = 1.5498 Å is ~21.555°.

    From first principles: d_{100} = 1 / |b1|·2π = a · sin γ ≈ 4.144 Å,
    so 2θ = 2 · arcsin(λ / (2 d)) ≈ 21.555°.

    Pre-fix the package returned 2θ ≈ 18.640° (computed from the wrong
    |Q| = 1.3131 Å⁻¹).  Post-fix it returns the true value, agreeing
    with the independent ``hkl_soleil`` reference quoted in the issue.
    """
    lat = ahd.Lattice(**SAPPHIRE_KWARGS)
    b1, _, _ = lat.reciprocal_lattice_vectors
    Q_mag = float(np.linalg.norm(b1))

    wavelength = 1.5498  # Å (the value used in the issue's hkl_soleil cross-check)
    sin_theta = wavelength * Q_mag / (4.0 * math.pi)
    two_theta = 2.0 * math.degrees(math.asin(sin_theta))

    np.testing.assert_allclose(two_theta, 21.5551, atol=1e-3)


def test_sapphire_forward_2theta_matches_hkl_soleil():
    """Full forward() for sapphire (1, 0, 0) on fourcv matches the hkl_soleil reference.

    The issue tabulates the hkl_soleil cross-check for sapphire on E4CV /
    bissector mode at λ = 1.5498 Å, with the diffractometer oriented from
    two observed reflections r1 = (0, 0, 6) and r2 = (1, 0, 0).  Pre-fix,
    ad_hoc_diffractometer reported 2θ = 18.6395° for (1, 0, 0); hkl_soleil
    reported 21.5551°.  Post-fix, ad_hoc_diffractometer agrees with
    hkl_soleil to within 1e-3°.

    This test does not orient against r1/r2 (which would require numerical
    solvers); it uses ``ub_identity`` and verifies that the resulting 2θ
    for (1, 0, 0) matches the analytical Bragg value.
    """
    g = fourcv()
    g.wavelength = 1.5498
    g.sample.lattice = ahd.Lattice(**SAPPHIRE_KWARGS)
    ahd.ub_identity(g.sample)
    g.mode_name = "bisecting"

    solutions = g.forward(1, 0, 0)
    assert len(solutions) > 0, "forward(1, 0, 0) must return at least one solution"

    # All solutions must report the analytical Bragg 2θ.
    for sol in solutions:
        np.testing.assert_allclose(sol["ttheta"], 21.5551, atol=1e-3)
