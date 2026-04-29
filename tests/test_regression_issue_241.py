# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Cross-module regression tests for issue #241.

The kappa virtual-bisecting solver returned ``"No solutions"`` for
several physically valid reflections that the equivalent Eulerian
``fourcv/bisecting`` solver handled cleanly.  The diagnostic that
uncovered the bug demonstrated two independent failures:

1. ``eulerian_to_kappa(omega, chi, phi)`` (the textbook Walko (2016)
   eq. [16] formula) does **not** preserve the scattering vector
   when the kappa preset's outer-stage axis is not aligned with the
   textbook ``+vertical`` convention.  ``kappa4cv`` (BL,
   ``-TRANSVERSE``), ``kappa4ch`` (BL, ``-VERTICAL``), and
   ``kappa6c`` (You, ``-TRANSVERSE``) all encode mixed-handedness
   stacks that the textbook formula cannot describe.

2. The kappa-axis vector was historically computed as
   ``vertical·cos(α) + transverse·sin(α)`` regardless of which
   physical axis the preset's outer komega rotates about.  The
   correct definition is ``n_kappa = cos(α)·n_komega +
   sin(α)·n_chi_eq`` — that is, tilted from the *actual* komega axis
   toward the *actual* equivalent-Eulerian chi axis.

The fix introduces a per-geometry
:class:`~ad_hoc_diffractometer.kappa.KappaPseudoAngleConvention` and
two geometry-aware functions
:func:`~ad_hoc_diffractometer.kappa.eulerian_to_kappa_axes` and
:func:`~ad_hoc_diffractometer.kappa.kappa_to_eulerian_axes` that
derive the pseudoangle relations directly from the preset's signed
stage axes.  Each kappa preset declares its convention.

This test file lives at the cross-module level because the bug spans
``kappa.py`` (pseudoangle conversions), ``presets.py`` (kappa-axis
construction), ``mode.py`` (``VirtualBisectConstraint.evaluate``),
``forward.py`` (the kappa virtual-mode dispatcher), and ``scan.py``
(``_kappa_from_Z``).
"""

from __future__ import annotations

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.kappa import eulerian_to_kappa_axes
from ad_hoc_diffractometer.kappa import kappa_to_eulerian_axes
from ad_hoc_diffractometer.orientation import angles_to_phi_vector

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _setup_cubic_eulerian(factory, a=4.0, wavelength=1.5406):
    """Cubic crystal geometry for direct comparison with kappa equivalent."""
    g = factory()
    g.wavelength = wavelength
    g.sample.lattice = ahd.Lattice(a=a)
    ahd.ub_identity(g.sample)
    return g


# ---------------------------------------------------------------------------
# Q-equivalence: kappa-stack matches Eulerian-stack across the χ range
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory_pair",
    [
        pytest.param((ahd.presets.fourcv, ahd.presets.kappa4cv), id="fourcv-kappa4cv"),
        pytest.param((ahd.presets.fourch, ahd.presets.kappa4ch), id="fourch-kappa4ch"),
        pytest.param((ahd.presets.psic, ahd.presets.kappa6c), id="psic-kappa6c"),
    ],
)
@pytest.mark.parametrize(
    "omega, chi, phi, ttheta, context",
    [
        pytest.param(0, 0, 0, 30, does_not_raise(), id="origin"),
        pytest.param(10, 0, 0, 30, does_not_raise(), id="omega-only"),
        pytest.param(0, 0, 10, 30, does_not_raise(), id="phi-only"),
        pytest.param(0, 5, 0, 30, does_not_raise(), id="small-chi"),
        pytest.param(0, 10, 0, 30, does_not_raise(), id="moderate-chi"),
        pytest.param(0, 30, 0, 30, does_not_raise(), id="chi-30"),
        pytest.param(0, 60, 0, 30, does_not_raise(), id="chi-60"),
        pytest.param(0, 80, 0, 30, does_not_raise(), id="chi-80"),
        pytest.param(0, 89, 0, 30, does_not_raise(), id="chi-89"),
        pytest.param(0, 99, 0, 30, does_not_raise(), id="chi-99"),
        pytest.param(0, -30, 0, 30, does_not_raise(), id="chi-negative-30"),
        pytest.param(0, -89, 0, 30, does_not_raise(), id="chi-negative-89"),
        pytest.param(10, 20, 30, 30, does_not_raise(), id="general-1"),
        pytest.param(5, 45, 15, 30, does_not_raise(), id="general-2"),
        pytest.param(45, 30, 60, 30, does_not_raise(), id="general-3"),
        pytest.param(6.85, 84.23, -13.35, 13.7, does_not_raise(), id="sapphire-002"),
    ],
)
@pytest.mark.parametrize("branch", [+1, -1])
def test_eulerian_to_kappa_axes_preserves_q(
    factory_pair,
    omega,
    chi,
    phi,
    ttheta,
    branch,
    context,
):
    """The geometry-aware decomposition preserves the scattering vector.

    For every kappa preset and a wide range of virtual Eulerian
    pseudoangles, ``eulerian_to_kappa_axes`` must return a kappa
    motor triple whose physical Q vector matches the Q vector of the
    sister Eulerian preset at the same pseudoangles.

    Pre-fix (issue #241) this property held only at ``chi = 0``;
    every non-trivial chi produced a Q diff of order unity.
    """
    eulerian_factory, kappa_factory = factory_pair
    g_eul = eulerian_factory()
    g_kap = kappa_factory()
    g_eul.wavelength = 1.5
    g_kap.wavelength = 1.5

    convention = g_kap.kappa_pseudo_angle_convention
    sample_names_eul = [s.name for s in g_eul.sample_stages]
    detector_name = g_eul.detector_stages[-1].name

    # Build the Eulerian motor dict.  The sister Eulerian geometry
    # may have outer/inner stages beyond the (omega, chi, phi) triple
    # (e.g. mu on psic) — set those to zero so the comparison is
    # apples-to-apples.
    eul_angles = {n: 0.0 for n in sample_names_eul}
    eul_angles[detector_name] = ttheta
    # The bottom three sample stages on the Eulerian preset map onto
    # (omega, chi, phi) by name in fourcv/fourch and onto (eta, chi,
    # phi) on psic — locate them positionally.
    eul_angles[sample_names_eul[-3]] = omega
    eul_angles[sample_names_eul[-2]] = chi
    eul_angles[sample_names_eul[-1]] = phi
    # Outer detector stages (e.g. nu on psic) stay at zero.
    for s in g_eul.detector_stages[:-1]:
        eul_angles[s.name] = 0.0

    Q_eul = angles_to_phi_vector(g_eul, **eul_angles)

    with context:
        ko, k, kp = eulerian_to_kappa_axes(omega, chi, phi, convention, branch=branch)

    # If the requested orientation is outside the kappa arm's
    # reachable range the test point is skipped — that is a feature
    # of the kappa geometry, not a bug.
    sample_names_kap = [s.name for s in g_kap.sample_stages]
    detector_name_kap = g_kap.detector_stages[-1].name
    kap_angles = {n: 0.0 for n in sample_names_kap}
    kap_angles[detector_name_kap] = ttheta
    kap_angles[sample_names_kap[-3]] = ko
    kap_angles[sample_names_kap[-2]] = k
    kap_angles[sample_names_kap[-1]] = kp
    for s in g_kap.detector_stages[:-1]:
        kap_angles[s.name] = 0.0

    Q_kap = angles_to_phi_vector(g_kap, **kap_angles)

    np.testing.assert_allclose(Q_kap, Q_eul, atol=1e-10)


# ---------------------------------------------------------------------------
# Round-trip: eulerian_to_kappa_axes ∘ kappa_to_eulerian_axes = identity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(ahd.presets.kappa4cv, id="kappa4cv"),
        pytest.param(ahd.presets.kappa4ch, id="kappa4ch"),
        pytest.param(ahd.presets.kappa6c, id="kappa6c"),
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
    """``eulerian → kappa → eulerian`` is the identity for every preset.

    Verifies the per-geometry decomposition is invertible at the +1
    branch (the natural identity branch).  The −1 branch round-trips
    onto the chi-mirrored point because the geometry-aware
    decomposition encodes the chi sign in the kappa branch parameter.
    """
    g = factory()
    convention = g.kappa_pseudo_angle_convention
    ko, k, kp = eulerian_to_kappa_axes(omega, chi, phi, convention, branch=branch)
    om, ch, ph = kappa_to_eulerian_axes(ko, k, kp, convention)
    # The recovered pseudoangles match to machine precision.
    assert om == pytest.approx(omega, abs=1e-10)
    # When chi=0 the (omega, phi) split is degenerate; the identity
    # still holds via the explicit ``parallel_axes`` short-circuit.
    if abs(chi) < 1e-10:
        assert ch == pytest.approx(0.0, abs=1e-10)
        assert ph == pytest.approx(phi, abs=1e-10)
    else:
        assert ch == pytest.approx(chi, abs=1e-10)
        assert ph == pytest.approx(phi, abs=1e-10)


# ---------------------------------------------------------------------------
# Sapphire reflections from the original issue
# ---------------------------------------------------------------------------

# The reproducer in #241 used hexagonal sapphire (a=4.785, c=12.991,
# γ=120°) at λ=1.5498 Å oriented from r1=(0,0,6) and r2=(1,0,0).
# Pre-fix, several reflections returned ``"No solutions"`` on
# ``kappa4cv/bisecting`` and ``kappa6c/bisecting_vertical``; the
# canonical failures were (0,0,2) and (0,0,4) on both kappa
# geometries (reachable in fourcv) and (1,1,0) and (1,1,3) on
# kappa6c only.

LAMBDA_SAPPHIRE = 1.5498
SAPPHIRE_LATTICE = dict(a=4.785, b=4.785, c=12.991, alpha=90.0, beta=90.0, gamma=120.0)
EULER_R1 = (20.9709, 89.5521, -0.0047, 41.9418)  # (omega, chi, phi, ttheta) for (0,0,6)
EULER_R2 = (9.3197, -7.3547, 2.6644, 21.5551)  # (omega, chi, phi, ttheta) for (1,0,0)


def _setup_sapphire_eulerian(factory, mode_name="bisecting"):
    """Sapphire setup on a 4-circle Eulerian preset.

    Uses the literal ``EULER_R1``/``EULER_R2`` motor tuples
    ``(omega, chi, phi, ttheta)`` from the issue reproducer.  Only
    valid for ``fourcv`` and ``fourch`` (4-circle).  For ``psic`` and
    other extended-stack Eulerian sister geometries, the seed
    reflections must be supplied by name with the outer/inner stages
    set to zero — see :func:`_setup_sapphire_psic` below.
    """
    g = factory()
    g.wavelength = LAMBDA_SAPPHIRE
    g.add_sample("sapphire", lattice=ahd.Lattice(**SAPPHIRE_LATTICE))
    g.sample = "sapphire"
    sample_names = [s.name for s in g.sample_stages]
    det = g.detector_stages[-1].name
    keys = [*sample_names, det]
    g.add_reflection(
        "r1",
        (0, 0, 6),
        dict(zip(keys, EULER_R1, strict=False)),
        wavelength=LAMBDA_SAPPHIRE,
    )
    g.add_reflection(
        "r2",
        (1, 0, 0),
        dict(zip(keys, EULER_R2, strict=False)),
        wavelength=LAMBDA_SAPPHIRE,
    )
    ahd.ub_from_two_reflections_bl1967(g.sample, "r1", "r2")
    g.mode_name = mode_name
    return g


def _setup_sapphire_psic(mode_name="bisecting_vertical"):
    """Sapphire setup on the psic 6-circle preset.

    Maps the 4-circle EULER_R*[0..3] = ``(omega, chi, phi, ttheta)``
    tuples onto psic's ``(mu=0, eta=omega, chi, phi, nu=0,
    delta=ttheta)`` motor dict.  In ``bisecting_vertical`` with
    ``mu=nu=0`` the psic stack is mathematically equivalent to
    fourcv.
    """
    g = ahd.presets.psic()
    g.wavelength = LAMBDA_SAPPHIRE
    g.add_sample("sapphire", lattice=ahd.Lattice(**SAPPHIRE_LATTICE))
    g.sample = "sapphire"
    g.add_reflection(
        "r1",
        (0, 0, 6),
        dict(
            mu=0.0,
            eta=EULER_R1[0],
            chi=EULER_R1[1],
            phi=EULER_R1[2],
            nu=0.0,
            delta=EULER_R1[3],
        ),
        wavelength=LAMBDA_SAPPHIRE,
    )
    g.add_reflection(
        "r2",
        (1, 0, 0),
        dict(
            mu=0.0,
            eta=EULER_R2[0],
            chi=EULER_R2[1],
            phi=EULER_R2[2],
            nu=0.0,
            delta=EULER_R2[3],
        ),
        wavelength=LAMBDA_SAPPHIRE,
    )
    ahd.ub_from_two_reflections_bl1967(g.sample, "r1", "r2")
    g.mode_name = mode_name
    return g


def _setup_sapphire_kappa(factory, mode_name):
    """Sapphire setup on a kappa preset — uses the geometry-aware
    pseudoangle conversion to translate the seed Eulerian motor angles
    into kappa motor angles for the reflection orientations.
    """
    g = factory()
    g.wavelength = LAMBDA_SAPPHIRE
    g.add_sample("sapphire", lattice=ahd.Lattice(**SAPPHIRE_LATTICE))
    g.sample = "sapphire"
    convention = g.kappa_pseudo_angle_convention
    kr1 = eulerian_to_kappa_axes(*EULER_R1[:3], convention, branch=+1)
    kr2 = eulerian_to_kappa_axes(*EULER_R2[:3], convention, branch=+1)
    sample_names = [s.name for s in g.sample_stages]
    det = g.detector_stages[-1].name
    if len(sample_names) == 3:  # kappa4cv / kappa4ch
        r1_angles = dict(zip(sample_names + [det], [*kr1, EULER_R1[3]], strict=False))
        r2_angles = dict(zip(sample_names + [det], [*kr2, EULER_R2[3]], strict=False))
    else:  # kappa6c (mu, komega, kappa, kphi); detector (nu, delta)
        r1_angles = dict(
            mu=0.0, komega=kr1[0], kappa=kr1[1], kphi=kr1[2], nu=0.0, delta=EULER_R1[3]
        )
        r2_angles = dict(
            mu=0.0, komega=kr2[0], kappa=kr2[1], kphi=kr2[2], nu=0.0, delta=EULER_R2[3]
        )
    g.add_reflection("r1", (0, 0, 6), r1_angles, wavelength=LAMBDA_SAPPHIRE)
    g.add_reflection("r2", (1, 0, 0), r2_angles, wavelength=LAMBDA_SAPPHIRE)
    ahd.ub_from_two_reflections_bl1967(g.sample, "r1", "r2")
    g.mode_name = mode_name
    return g


@pytest.mark.parametrize(
    "hkl",
    [
        pytest.param((0, 0, 2), id="002"),
        pytest.param((0, 0, 4), id="004"),
        pytest.param((0, 0, 6), id="006"),
        pytest.param((1, 0, 0), id="100"),
        pytest.param((1, 0, 2), id="102"),
        pytest.param((1, 1, 0), id="110"),
        pytest.param((2, 0, 0), id="200"),
        pytest.param((2, 1, 0), id="210"),
        pytest.param((0, 1, 2), id="012"),
    ],
)
def test_sapphire_kappa4cv_bisecting_solves(hkl):
    """Every reflection that fourcv/bisecting solves on sapphire is
    also solved by kappa4cv/bisecting and produces the same |2θ|.

    Pre-fix (issue #241), (0,0,2) and (0,0,4) returned "No solutions"
    on kappa4cv/bisecting although fourcv/bisecting accepted them.
    """
    fc = _setup_sapphire_eulerian(ahd.presets.fourcv, "bisecting")
    fc_sols = fc.forward(*hkl)
    if not fc_sols:
        pytest.skip(f"{hkl} not reachable on fourcv/bisecting; comparison vacuous")

    k4 = _setup_sapphire_kappa(ahd.presets.kappa4cv, "bisecting")
    k4_sols = k4.forward(*hkl)
    assert k4_sols, f"kappa4cv/bisecting returned no solutions for {hkl} (issue #241)"
    # Compare |2θ| magnitudes; both libraries-independent choices.
    fc_tt = abs(fc_sols[0]["ttheta"])
    k4_tt = abs(k4_sols[0]["ttheta"])
    assert k4_tt == pytest.approx(fc_tt, abs=1e-9)


@pytest.mark.parametrize(
    "hkl",
    [
        pytest.param((0, 0, 2), id="002"),
        pytest.param((0, 0, 4), id="004"),
        pytest.param((0, 0, 6), id="006"),
        pytest.param((1, 0, 0), id="100"),
        pytest.param((2, 0, 0), id="200"),
        pytest.param((2, 1, 0), id="210"),
        pytest.param((0, 1, 2), id="012"),
    ],
)
def test_sapphire_kappa6c_bisecting_vertical_solves(hkl):
    """Every reflection that psic/bisecting_vertical solves on sapphire
    is also solved by kappa6c/bisecting_vertical and produces the
    same |delta|.

    Pre-fix (issue #241), (0,0,2), (0,0,4) and (1,1,0), (1,1,3)
    returned "No solutions" on kappa6c/bisecting_vertical although
    the equivalent psic bisecting solver accepted them.
    """
    psic = _setup_sapphire_psic("bisecting_vertical")
    psic_sols = psic.forward(*hkl)
    if not psic_sols:
        pytest.skip(
            f"{hkl} not reachable on psic/bisecting_vertical; comparison vacuous"
        )

    k6 = _setup_sapphire_kappa(ahd.presets.kappa6c, "bisecting_vertical")
    k6_sols = k6.forward(*hkl)
    assert k6_sols, (
        f"kappa6c/bisecting_vertical returned no solutions for {hkl} (issue #241)"
    )
    psic_dd = abs(psic_sols[0]["delta"])
    k6_dd = abs(k6_sols[0]["delta"])
    assert k6_dd == pytest.approx(psic_dd, abs=1e-9)


# ---------------------------------------------------------------------------
# Fixed-virtual-angle modes — exercised by the same Eulerian-equivalent
# infrastructure as the bisecting mode (issue #241).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(ahd.presets.kappa4cv, id="kappa4cv"),
        pytest.param(ahd.presets.kappa4ch, id="kappa4ch"),
    ],
)
@pytest.mark.parametrize(
    "mode_name, virtual_angle, target_value, hkl",
    [
        pytest.param("fixed_omega", "omega", 0.0, (0, 0, 1), id="fixed_omega"),
        pytest.param("fixed_chi", "chi", 90.0, (0, 0, 1), id="fixed_chi"),
        pytest.param("fixed_phi", "phi", 0.0, (0, 1, 0), id="fixed_phi"),
    ],
)
def test_fixed_virtual_angle_modes(
    factory, mode_name, virtual_angle, target_value, hkl
):
    """Each fixed virtual-angle mode returns at least one solution
    whose virtual angle (recovered via the geometry-aware
    decomposition) matches the constraint value to 1e-4°.

    Issue #241: pre-fix, the textbook ``kappa_to_eulerian`` validation
    in the dispatcher silently rejected solutions whose Walko-frame
    virtual angle disagreed with the geometry-aware angle.
    """
    g = _setup_cubic_eulerian(factory, a=4.0)
    g.mode_name = mode_name
    sols = g.forward(*hkl)
    assert sols, f"{factory.__name__}/{mode_name} returned no solutions for {hkl}"
    convention = g.kappa_pseudo_angle_convention
    for sol in sols:
        om, ch, ph = kappa_to_eulerian_axes(
            sol["komega"], sol["kappa"], sol["kphi"], convention
        )
        recovered = {"omega": om, "chi": ch, "phi": ph}[virtual_angle]
        assert recovered == pytest.approx(target_value, abs=1e-4)
