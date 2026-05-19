# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
Cross-module regression tests for issue #252.

The kappa axis vector for ``kappa4cv``, ``kappa4ch``, and ``kappa6c``
was incorrect under the v0.9.1 convention introduced in #241/#247.
The corrected convention follows the published literature:

* ``kappa4cv``: kappa axis in the transverse-vertical plane, between
  +T and +V (Walko 2016 Fig. 3; Thorkildsen et al. 2006 Table 1).
* ``kappa4ch``: kappa axis in the vertical-longitudinal plane,
  between +V and +L (Wyckoff 1985 Fig. 2(b) on p. 334).
* ``kappa6c``:  kappa axis as for kappa4cv (Sønsteby et al. 2013;
  Thorkildsen 2006 §3).

The fix lives entirely in ``presets.py``; the closed-form solver in
``kappa.py`` already worked from the four signed stage axes and
required no math change.

This file lives at the cross-module level because the bug spans
``presets.py`` (kappa-axis construction), ``kappa.py`` (pseudoangle
conversions consumed by the solver), ``mode.py``
(``VirtualBisectConstraint.evaluate``), ``forward.py`` (the kappa
virtual-mode dispatcher), and ``scan.py`` (``_kappa_from_Z``).
"""

from __future__ import annotations

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from helpers import kappa4ch
from helpers import kappa4cv
from helpers import kappa6c

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.kappa import eulerian_to_kappa_axes
from ad_hoc_diffractometer.kappa import kappa_to_eulerian_axes

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _setup_cubic(factory, a=4.0, wavelength=1.5406):
    """Cubic crystal geometry with UB = identity."""
    g = factory()
    g.wavelength = wavelength
    g.sample.lattice = ahd.Lattice(a=a)
    ahd.ub_identity(g.sample)
    return g


# ---------------------------------------------------------------------------
# Kappa axis vector matches published literature
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, plane_directions, alpha_deg, context",
    [
        pytest.param(
            kappa4cv,
            ("transverse", "vertical"),
            50.0,
            does_not_raise(),
            id="kappa4cv-walko2016-fig3",
        ),
        pytest.param(
            kappa4ch,
            ("vertical", "longitudinal"),
            50.0,
            does_not_raise(),
            id="kappa4ch-wyckoff1985-fig2b",
        ),
        pytest.param(
            kappa6c,
            ("transverse", "vertical"),
            50.0,
            does_not_raise(),
            id="kappa6c-walko2016-fig3",
        ),
        # Custom alpha values — the kappa axis must still lie in the
        # documented basis-direction plane for any α in (0, 90°).
        pytest.param(
            kappa4cv,
            ("transverse", "vertical"),
            45.0,
            does_not_raise(),
            id="kappa4cv-alpha45",
        ),
        pytest.param(
            kappa4ch,
            ("vertical", "longitudinal"),
            35.0,
            does_not_raise(),
            id="kappa4ch-alpha35",
        ),
        pytest.param(
            kappa6c,
            ("transverse", "vertical"),
            55.0,
            does_not_raise(),
            id="kappa6c-alpha55",
        ),
    ],
)
def test_kappa_axis_matches_published_literature(
    factory, plane_directions, alpha_deg, context
):
    """The kappa stage axis matches ``cos(α)·v_from + sin(α)·v_to``,
    where ``v_from`` and ``v_to`` are the two physical-direction basis
    vectors that span the documented kappa-arm plane.

    This is the single, most direct regression target for issue #252:
    if any future refactor breaks the published kappa-axis convention
    we documented here, this test fires immediately.
    """
    with context:
        g = factory(alpha_deg=alpha_deg)
        from_name, to_name = plane_directions
        from_vec = np.asarray(g.basis[from_name], dtype=float)
        to_vec = np.asarray(g.basis[to_name], dtype=float)
        cos_a = np.cos(np.deg2rad(alpha_deg))
        sin_a = np.sin(np.deg2rad(alpha_deg))
        expected = cos_a * from_vec + sin_a * to_vec
        np.testing.assert_allclose(
            g.stage("kappa").axis,
            expected,
            atol=1e-12,
            err_msg=(
                f"{factory.__name__}: kappa axis must lie in the "
                f"({from_name}, {to_name}) plane, between +{from_name} "
                f"and +{to_name}, tilted α={alpha_deg}° from +{from_name}."
            ),
        )


@pytest.mark.parametrize(
    "factory, plane_directions, normal_direction, context",
    [
        pytest.param(
            kappa4cv,
            ("transverse", "vertical"),
            "longitudinal",
            does_not_raise(),
            id="kappa4cv-normal-to-L",
        ),
        pytest.param(
            kappa4ch,
            ("vertical", "longitudinal"),
            "transverse",
            does_not_raise(),
            id="kappa4ch-normal-to-T",
        ),
        pytest.param(
            kappa6c,
            ("transverse", "vertical"),
            "longitudinal",
            does_not_raise(),
            id="kappa6c-normal-to-L",
        ),
    ],
)
def test_kappa_axis_orthogonal_to_third_basis_direction(
    factory, plane_directions, normal_direction, context
):
    """Confirms the kappa axis has *exactly zero* component along the
    basis direction perpendicular to its containing plane.  This
    catches the v0.9.1 (#247) regression where the kappa axis had a
    non-zero longitudinal component for kappa4cv/kappa6c.
    """
    with context:
        g = factory()
        kappa = np.asarray(g.stage("kappa").axis, dtype=float)
        normal_vec = np.asarray(g.basis[normal_direction], dtype=float)
        # The kappa axis must be perpendicular to the normal direction
        # within machine precision.
        dot = float(np.dot(kappa, normal_vec))
        assert abs(dot) < 1e-12, (
            f"{factory.__name__}: kappa axis has a component "
            f"({dot:.3e}) along {normal_direction}, but should lie "
            f"in the ({plane_directions[0]}, {plane_directions[1]}) plane."
        )


# ---------------------------------------------------------------------------
# Round-trip: eulerian_to_kappa_axes ∘ kappa_to_eulerian_axes = identity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(kappa4cv, id="kappa4cv"),
        pytest.param(kappa4ch, id="kappa4ch"),
        pytest.param(kappa6c, id="kappa6c"),
    ],
)
@pytest.mark.parametrize(
    "omega, chi, phi, branch",
    [
        pytest.param(10, 30, 45, +1, id="general-pos"),
        pytest.param(20, 60, 70, +1, id="general-2-pos"),
        pytest.param(0, 89, 0, +1, id="near-limit-pos"),
        pytest.param(15, -45, 25, +1, id="negative-chi-pos-branch"),
        pytest.param(0, 0, 0, +1, id="origin-pos"),
        pytest.param(11.10, 0, 90, +1, id="chi0-degenerate-pos"),
    ],
)
def test_eulerian_kappa_round_trip(factory, omega, chi, phi, branch):
    """Round-trip ``kappa_to_eulerian_axes ∘ eulerian_to_kappa_axes``
    is the identity on the +1 branch for every kappa preset.

    The closed-form solver in ``kappa.py`` works directly from the
    four signed stage axes stored on the geometry; the round-trip
    invariant must hold regardless of the kappa-axis vector
    convention.  This is the strongest internal-consistency guard
    for the kappa pseudoangle layer.

    The −1 branch round-trips onto the chi-mirrored point and is
    therefore not tested as an identity here (the geometry-aware
    decomposition encodes the chi sign in the kappa branch
    parameter).
    """
    g = factory()
    convention = g.kappa_pseudo_angle_convention
    ko, k, kp = eulerian_to_kappa_axes(omega, chi, phi, convention, branch=branch)
    om, ch, ph = kappa_to_eulerian_axes(ko, k, kp, convention)
    assert om == pytest.approx(omega, abs=1e-10)
    if abs(chi) < 1e-10:
        assert ch == pytest.approx(0.0, abs=1e-10)
        assert ph == pytest.approx(phi, abs=1e-10)
    else:
        assert ch == pytest.approx(chi, abs=1e-10)
        assert ph == pytest.approx(phi, abs=1e-10)


# ---------------------------------------------------------------------------
# Reachability: every reflection reachable on the bisecting kappa
# preset round-trips back to the input (h, k, l).
#
# This replaces the pre-fix "kappa solver returned no solutions for
# physically reachable reflections" diagnostic from issue #241.  The
# specific reflections that v0.9.1 mishandled are not the same set
# the corrected convention mishandles (the "reachable" set itself
# shifts when the kappa-axis convention changes), so the test
# exercises a representative sample of cubic reflections rather than
# the original sapphire reproducer.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, mode_name",
    [
        pytest.param(kappa4cv, "bisecting", id="kappa4cv-bisecting"),
        pytest.param(kappa4ch, "bisecting", id="kappa4ch-bisecting"),
        pytest.param(
            kappa6c,
            "bisecting_vertical",
            id="kappa6c-bisecting_vertical",
        ),
    ],
)
@pytest.mark.parametrize(
    "hkl",
    [
        pytest.param((1, 0, 0), id="100"),
        pytest.param((0, 1, 0), id="010"),
        pytest.param((0, 0, 1), id="001"),
        pytest.param((1, 1, 0), id="110"),
        pytest.param((1, 0, 1), id="101"),
        pytest.param((0, 1, 1), id="011"),
        pytest.param((1, 1, 1), id="111"),
        pytest.param((2, 0, 0), id="200"),
        pytest.param((2, 1, 0), id="210"),
        pytest.param((2, 1, 1), id="211"),
    ],
)
def test_bisecting_reachable_reflections_round_trip(factory, mode_name, hkl):
    """Every reflection reachable in the bisecting mode round-trips
    through ``forward`` → ``inverse``.

    A reflection that the bisecting solver cannot reach (because of
    kappa-arm geometry, motor limits, or the Ewald sphere) is skipped
    rather than failing — those are properties of the geometry, not
    bugs.
    """
    g = _setup_cubic(factory)
    g.mode_name = mode_name
    sols = g.forward(*hkl)
    if not sols:
        pytest.skip(
            f"{factory.__name__}/{mode_name} cannot reach {hkl} "
            f"(kappa-arm or motor-limit geometry)."
        )
    for sol in sols:
        hkl_back = g.inverse(sol)
        np.testing.assert_allclose(hkl_back, hkl, atol=1e-8)


# ---------------------------------------------------------------------------
# Sapphire reachability — historical issue-#241 reproducer.
#
# The original diagnostic in #241 was that several sapphire
# reflections returned "No solutions" on kappa4cv/bisecting although
# fourcv/bisecting accepted them.  Under the corrected convention
# (#252) the kappa↔fourcv equivalence is no longer 1:1 — the kappa
# arm physically extends in a different plane, so "reachable on
# fourcv" and "reachable on kappa4cv" describe different subsets of
# reciprocal space.  This test therefore checks only:
#   (i)  the sapphire (0,0,2) reflection that originally triggered
#        #241 still solves on kappa4cv, and
#   (ii) every solution round-trips.
# ---------------------------------------------------------------------------


def test_sapphire_100_kappa4cv_bisecting_solves():
    """A sapphire reflection in the issue-#241 family still solves on
    kappa4cv/bisecting after the issue-#252 axis correction and the
    issue-#280 ub_identity update.

    Sapphire: a=4.785 Å, c=12.991 Å, γ=120°; λ=1.5498 Å; UB=identity.

    The original #241 reproducer was sapphire ``(0, 0, 2)``.  Under
    issue #280 ub_identity (crystal a* axis physically along the beam,
    c* axis physically along transverse), ``(0, 0, 2)`` puts Q_phi
    along the kappa4cv kphi-axis direction and is physically
    unreachable in bisecting on this geometry.  ``(1, 0, 0)`` is the
    equivalent reachable sapphire reflection under the new
    ub_identity and exercises the same kappa decomposition code path.
    """
    g = kappa4cv()
    g.wavelength = 1.5498
    g.sample.lattice = ahd.Lattice(
        a=4.785, b=4.785, c=12.991, alpha=90.0, beta=90.0, gamma=120.0
    )
    ahd.ub_identity(g.sample)
    g.mode_name = "bisecting"
    sols = g.forward(1, 0, 0)
    assert sols, (
        "sapphire (1,0,0) must solve on kappa4cv/bisecting — exercises "
        "the same code path as the original issue-#241 reproducer."
    )
    for sol in sols:
        np.testing.assert_allclose(g.inverse(sol), [1, 0, 0], atol=1e-8)
