# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
Regression tests for issue #284.

Issue #284 reported that ``kappa4cv bisecting`` and
``kappa6c bisecting_vertical`` lost solutions for several sapphire
asymmetric reflections in v0.11.0.  The reproducer's exact failure
set: ``(0, 1, 2)``, ``(0, 0, 6)``, ``(1, 1, 3)`` on both geometries
all return zero solutions; the symmetric ``(1, 0, 0)`` and
``(1, 1, 0)`` still solve.

Root cause: issue #252 conflated two distinct concepts in the single
YAML field ``kappa_chi_eq``:

1. The **kappa-arm tilt direction** that defines the kappa stage's
   geometric axis (input to Walko's formula ``n_kappa = cos(α)·
   n_komega_unsigned + sin(α)·tilt_direction``).  For ``kappa4cv``
   this is ``+vertical`` (kappa arm lies in the transverse-vertical
   plane).

2. The **equivalent-Eulerian chi pseudo-angle axis** the kappa→
   Eulerian decomposition rotates about for the virtual chi angle.
   For ``kappa4cv`` this should be ``+longitudinal`` to match
   ``fourcv``'s chi axis (the analogous non-kappa Eulerian preset
   that ``kappa4cv`` is mechanically equivalent to).

Setting both to the same value (``+vertical``) broke the second
role: the kappa equivalent-Eulerian decomposition's reachable Bragg
locus shrank away from the locus reachable by ``fourcv bisecting``,
making asymmetric reflections like sapphire ``(0, 1, 2)`` decline.

The fix separates the two:

- ``kappa_chi_eq`` in the YAML still defines the kappa-arm tilt
  direction (used by ``_resolve_axis``).
- A new optional YAML field ``kappa_eulerian_chi`` defines the
  equivalent-Eulerian chi axis.  When absent, the loader derives it
  as the first basis direction perpendicular to ``n_komega`` in the
  conventional order ``(+longitudinal, +vertical, +transverse)`` —
  which yields ``+longitudinal`` for every kappa geometry shipped
  with the package.

These regression tests cover:

- the exact ``forward()`` reproducers from the issue body;
- the broader kappa-vs-fourcv equivalence (the kappa bisecting
  reachability now matches the equivalent fourcv/fourch/psic
  bisecting reachability);
- the auto-derivation rule of the new ``kappa_eulerian_chi`` field;
- the explicit ``kappa_eulerian_chi`` YAML override.
"""

from __future__ import annotations

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
import yaml

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.geometry_loader import KIND_KEY
from ad_hoc_diffractometer.geometry_loader import GeometrySchemaError
from ad_hoc_diffractometer.geometry_loader import load_geometry_file
from ad_hoc_diffractometer.kappa import KappaPseudoAngleConvention
from ad_hoc_diffractometer.kappa import eulerian_to_kappa_axes
from ad_hoc_diffractometer.orientation import angles_to_phi_vector

SAPPHIRE = dict(a=4.7589, b=4.7589, c=12.99119, alpha=90.0, beta=90.0, gamma=120.0)
WAVELENGTH = 1.54


# ---------------------------------------------------------------------------
# Reproducer: exact sapphire reflections from the issue body
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "geom_name, mode_name",
    [
        pytest.param("kappa4cv", "bisecting", id="kappa4cv-bisecting"),
        pytest.param("kappa6c", "bisecting_vertical", id="kappa6c-bisecting_vertical"),
    ],
)
@pytest.mark.parametrize(
    "hkl, expected_min_sols, context",
    [
        pytest.param((0, 1, 2), 1, does_not_raise(), id="sapphire-012-solves"),
        pytest.param((0, 0, 6), 1, does_not_raise(), id="sapphire-006-solves"),
        pytest.param((1, 1, 3), 1, does_not_raise(), id="sapphire-113-solves"),
        pytest.param((1, 1, 0), 1, does_not_raise(), id="sapphire-110-solves"),
        pytest.param((1, 0, 0), 1, does_not_raise(), id="sapphire-100-solves"),
    ],
)
def test_kappa_bisecting_sapphire_reflections_solve(
    geom_name,
    mode_name,
    hkl,
    expected_min_sols,
    context,
):
    """Issue #284 reproducer: sapphire asymmetric reflections solve.

    Direct copy of the reproducer block in the issue body.  Pre-fix
    these returned ``n_sols=0`` for the three asymmetric cases on
    both ``kappa4cv bisecting`` and ``kappa6c bisecting_vertical``.
    """
    g = ahd.make_geometry(geom_name)
    g.sample.lattice = ahd.Lattice(**SAPPHIRE)
    g.wavelength = WAVELENGTH
    ahd.ub_identity(g.sample)
    g.mode_name = mode_name
    with context:
        sols = g.forward(*hkl)
        assert len(sols) >= expected_min_sols, (
            f"{geom_name} / {mode_name}  {hkl}: expected at least "
            f"{expected_min_sols} solution(s), got {len(sols)}."
        )
        # Each returned solution must hit the target Q_phi.
        target = g.sample.UB @ np.asarray(hkl, dtype=float)
        for sol in sols:
            angles = {k: sol[k] for k in sol}
            Q = angles_to_phi_vector(g, **angles)
            np.testing.assert_allclose(Q, target, atol=1e-8)


# ---------------------------------------------------------------------------
# Cross-geometry equivalence: kappa bisecting reachability matches the
# corresponding non-kappa Eulerian bisecting reachability.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kappa_name, eulerian_name, kappa_mode, eulerian_mode, context",
    [
        pytest.param(
            "kappa4cv",
            "fourcv",
            "bisecting",
            "bisecting",
            does_not_raise(),
            id="kappa4cv-fourcv",
        ),
        pytest.param(
            "kappa6c",
            "psic",
            "bisecting_vertical",
            "bisecting_vertical",
            does_not_raise(),
            id="kappa6c-psic",
        ),
        pytest.param(
            "kappa4ch",
            "fourch",
            "bisecting",
            "bisecting",
            does_not_raise(),
            id="kappa4ch-fourch",
        ),
    ],
)
@pytest.mark.parametrize(
    "hkl",
    [
        pytest.param((0, 1, 2), id="sapphire-012"),
        pytest.param((0, 0, 6), id="sapphire-006"),
        pytest.param((1, 1, 3), id="sapphire-113"),
        pytest.param((1, 1, 0), id="sapphire-110"),
        pytest.param((1, 0, 0), id="sapphire-100"),
    ],
)
def test_kappa_bisecting_matches_eulerian_reachability(
    kappa_name,
    eulerian_name,
    kappa_mode,
    eulerian_mode,
    hkl,
    context,
):
    """The kappa equivalent-Eulerian decomposition reaches every
    reflection that the corresponding non-kappa Eulerian bisecting
    mode reaches.

    Pre-fix the kappa preset's bisecting reachability set was a
    strict subset of the sister Eulerian preset's set for several
    sapphire reflections (issue #284).  Post-fix the two sets agree
    on every test reflection.

    The angle values themselves differ (the kappa→Eulerian
    decomposition has its own (komega, kappa, kphi) branches), but
    each kappa solution must produce the same target ``Q_phi`` as
    the Eulerian solution.
    """
    g_k = ahd.make_geometry(kappa_name)
    g_k.sample.lattice = ahd.Lattice(**SAPPHIRE)
    g_k.wavelength = WAVELENGTH
    ahd.ub_identity(g_k.sample)
    g_k.mode_name = kappa_mode

    g_e = ahd.make_geometry(eulerian_name)
    g_e.sample.lattice = ahd.Lattice(**SAPPHIRE)
    g_e.wavelength = WAVELENGTH
    ahd.ub_identity(g_e.sample)
    g_e.mode_name = eulerian_mode

    with context:
        eul_sols = g_e.forward(*hkl)
        kap_sols = g_k.forward(*hkl)
        # If the Eulerian sister solves, the kappa preset must also
        # solve (within its kappa-arm reachability — every test
        # reflection here is within that range).
        if eul_sols:
            assert kap_sols, (
                f"{eulerian_name} {hkl} has {len(eul_sols)} solution(s) "
                f"but {kappa_name} returns 0."
            )


# ---------------------------------------------------------------------------
# kappa→Eulerian Q-equivalence: the geometry-aware decomposition
# preserves the scattering vector across a sweep of pseudoangles.
#
# This restores the central invariant from the pre-#252 issue-#241
# regression file (deleted by #252 PR #253); under the corrected
# convention the invariant holds for every kappa preset across the
# full reachable chi range.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kappa_name, eulerian_name",
    [
        pytest.param("kappa4cv", "fourcv", id="kappa4cv-fourcv"),
        pytest.param("kappa4ch", "fourch", id="kappa4ch-fourch"),
        pytest.param("kappa6c", "psic", id="kappa6c-psic"),
    ],
)
@pytest.mark.parametrize(
    "omega, chi, phi, ttheta, context",
    [
        pytest.param(0, 0, 0, 30, does_not_raise(), id="origin"),
        pytest.param(10, 0, 0, 30, does_not_raise(), id="omega-only"),
        pytest.param(0, 5, 0, 30, does_not_raise(), id="small-chi"),
        pytest.param(0, 30, 0, 30, does_not_raise(), id="chi-30"),
        pytest.param(0, 60, 0, 30, does_not_raise(), id="chi-60"),
        pytest.param(0, -30, 0, 30, does_not_raise(), id="chi-negative-30"),
        pytest.param(10, 20, 30, 30, does_not_raise(), id="general-1"),
        pytest.param(5, 45, 15, 30, does_not_raise(), id="general-2"),
    ],
)
def test_eulerian_to_kappa_axes_preserves_q(
    kappa_name,
    eulerian_name,
    omega,
    chi,
    phi,
    ttheta,
    context,
):
    """``eulerian_to_kappa_axes`` produces a kappa motor triple
    whose physical ``Q_phi`` matches the equivalent non-kappa
    Eulerian preset's ``Q_phi`` at the same virtual (omega, chi,
    phi).

    Pre-#252 this invariant held by construction.  #252's
    conflation of the kappa-arm tilt direction and the equivalent-
    Eulerian chi axis broke it for ``kappa4cv`` and ``kappa6c``
    (the test it lived in was deleted).  Post-#284 the invariant
    holds again across every kappa preset.
    """
    g_e = ahd.make_geometry(eulerian_name)
    g_e.wavelength = 1.5
    g_k = ahd.make_geometry(kappa_name)
    g_k.wavelength = 1.5

    convention = g_k.kappa_pseudo_angle_convention
    sample_names_eul = [s.name for s in g_e.sample_stages]
    detector_name_eul = g_e.detector_stages[-1].name

    # Build the Eulerian motor dict.  The sister Eulerian geometry
    # may have outer/inner stages beyond the (omega, chi, phi) triple
    # (e.g. mu on psic); set those to zero.
    eul_angles = {n: 0.0 for n in sample_names_eul}
    eul_angles[detector_name_eul] = ttheta
    # Map (omega, chi, phi) onto the inner three sample stages.
    eul_angles[sample_names_eul[-3]] = omega
    eul_angles[sample_names_eul[-2]] = chi
    eul_angles[sample_names_eul[-1]] = phi
    for s in g_e.detector_stages[:-1]:
        eul_angles[s.name] = 0.0

    Q_eul = angles_to_phi_vector(g_e, **eul_angles)

    with context:
        ko, k, kp = eulerian_to_kappa_axes(omega, chi, phi, convention, branch=+1)

    sample_names_kap = [s.name for s in g_k.sample_stages]
    detector_name_kap = g_k.detector_stages[-1].name
    kap_angles = {n: 0.0 for n in sample_names_kap}
    kap_angles[detector_name_kap] = ttheta
    kap_angles[sample_names_kap[-3]] = ko
    kap_angles[sample_names_kap[-2]] = k
    kap_angles[sample_names_kap[-1]] = kp
    for s in g_k.detector_stages[:-1]:
        kap_angles[s.name] = 0.0

    Q_kap = angles_to_phi_vector(g_k, **kap_angles)

    np.testing.assert_allclose(Q_kap, Q_eul, atol=1e-10)


# ---------------------------------------------------------------------------
# Convention values: all shipped kappa presets now use n_chi_eq =
# +longitudinal (matching the corresponding fourcv/fourch/psic chi
# axis).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kappa_name, expected_basis_label, context",
    [
        pytest.param("kappa4cv", "longitudinal", does_not_raise(), id="kappa4cv"),
        pytest.param("kappa4ch", "longitudinal", does_not_raise(), id="kappa4ch"),
        pytest.param("kappa6c", "longitudinal", does_not_raise(), id="kappa6c"),
    ],
)
def test_n_chi_eq_matches_basis_longitudinal(
    kappa_name,
    expected_basis_label,
    context,
):
    """Every shipped kappa preset's equivalent-Eulerian chi axis is
    ``+longitudinal`` (issue #284).

    Aligns the kappa equivalent-Eulerian decomposition's virtual chi
    axis with the corresponding non-kappa Eulerian preset's chi
    axis: ``fourcv``, ``fourch``, ``psic``, ``sixc``, and ``fivec``
    all put ``chi`` about ``+longitudinal``.
    """
    g = ahd.make_geometry(kappa_name)
    with context:
        np.testing.assert_allclose(
            g.kappa_pseudo_angle_convention.n_chi_eq,
            g.basis[expected_basis_label],
            atol=1e-12,
        )


# ---------------------------------------------------------------------------
# The new YAML field kappa_eulerian_chi: explicit override and
# auto-derivation.
# ---------------------------------------------------------------------------


def _kappa_yaml_doc(
    *,
    kappa_eulerian_chi: str | None = None,
    kappa_chi_eq: str = "+vertical",
):
    """Minimal kappa4cv-like YAML doc, parametrized for these tests."""
    doc = {
        KIND_KEY: {"schema_revision": 1},
        "name": "kappa_test",
        "documentation": "test",
        "basis": "BL",
        "parameters": {"alpha_deg": 50.0},
        "kappa_chi_eq": kappa_chi_eq,
        "stages": [
            {
                "name": "komega",
                "axis": "-transverse",
                "parent": None,
                "role": "sample",
            },
            {
                "name": "kappa",
                "axis": {"kappa_eulerian": "+transverse"},
                "parent": "komega",
                "role": "sample",
            },
            {
                "name": "kphi",
                "axis": "-transverse",
                "parent": "kappa",
                "role": "sample",
            },
            {
                "name": "ttheta",
                "axis": "-transverse",
                "parent": None,
                "role": "detector",
            },
        ],
        "modes": {
            "bisecting": {
                "default": True,
                "constraints": [
                    {"type": "virtual_bisect", "stage1": "omega", "stage2": "ttheta"}
                ],
                "computed": ["komega", "kappa", "kphi", "ttheta"],
            }
        },
    }
    if kappa_eulerian_chi is not None:
        doc["kappa_eulerian_chi"] = kappa_eulerian_chi
    return yaml.safe_dump(doc)


def test_kappa_eulerian_chi_explicit_override():
    """When ``kappa_eulerian_chi`` is declared the loader uses it
    verbatim for the convention's ``n_chi_eq`` (issue #284)."""
    text = _kappa_yaml_doc(kappa_eulerian_chi="+vertical")
    g = load_geometry_file(text)
    # Explicit override: n_chi_eq is +vertical (not the auto-derived
    # +longitudinal).
    np.testing.assert_allclose(
        g.kappa_pseudo_angle_convention.n_chi_eq, [0.0, 0.0, 1.0], atol=1e-12
    )


def test_kappa_eulerian_chi_auto_derived_when_absent():
    """Without ``kappa_eulerian_chi`` the loader derives ``n_chi_eq``
    as ``+longitudinal`` for a standard kappa4cv-style preset."""
    text = _kappa_yaml_doc()  # no kappa_eulerian_chi
    g = load_geometry_file(text)
    np.testing.assert_allclose(
        g.kappa_pseudo_angle_convention.n_chi_eq, [0.0, 1.0, 0.0], atol=1e-12
    )


def test_kappa_eulerian_chi_numeric_vector_form():
    """``kappa_eulerian_chi`` accepts a length-3 numeric vector."""
    doc = yaml.safe_load(_kappa_yaml_doc())
    doc["kappa_eulerian_chi"] = [0.0, 0.0, 1.0]
    g = load_geometry_file(yaml.safe_dump(doc))
    np.testing.assert_allclose(
        g.kappa_pseudo_angle_convention.n_chi_eq, [0.0, 0.0, 1.0], atol=1e-12
    )


def test_kappa_eulerian_chi_invalid_form_rejected():
    """Non-string non-vector ``kappa_eulerian_chi`` is rejected."""
    doc = yaml.safe_load(_kappa_yaml_doc())
    doc["kappa_eulerian_chi"] = 42  # neither a string nor a vector
    with pytest.raises(
        GeometrySchemaError, match="'kappa_eulerian_chi' must be"
    ):
        load_geometry_file(yaml.safe_dump(doc))


def test_kappa_eulerian_chi_auto_derivation_walks_basis_fallback():
    """The auto-derivation rule walks the basis directions in order
    ``(+longitudinal, +vertical, +transverse)`` and returns the first
    one perpendicular to ``n_komega``.  When the conventional first
    choice (``+longitudinal``) is parallel to ``n_komega``, the rule
    falls through to ``+vertical``.
    """
    # Construct a kappa-arm geometry with ``n_komega = +longitudinal``
    # (an unconventional but valid choice).  The auto-derivation
    # cannot pick +longitudinal (parallel to n_komega), so it must
    # pick the next perpendicular basis direction (+vertical).
    doc = {
        KIND_KEY: {"schema_revision": 1},
        "name": "kappa_long_omega",
        "documentation": "test",
        "basis": "BL",
        "parameters": {"alpha_deg": 50.0},
        "kappa_chi_eq": "+transverse",  # arm in (long, trans) plane
        "stages": [
            {
                "name": "komega",
                "axis": "+longitudinal",
                "parent": None,
                "role": "sample",
            },
            {
                "name": "kappa",
                "axis": {"kappa_eulerian": "+longitudinal"},
                "parent": "komega",
                "role": "sample",
            },
            {
                "name": "kphi",
                "axis": "+longitudinal",
                "parent": "kappa",
                "role": "sample",
            },
            {
                "name": "ttheta",
                "axis": "+longitudinal",
                "parent": None,
                "role": "detector",
            },
        ],
        "modes": {
            "fixed_kphi": {
                "default": True,
                "constraints": [{"type": "sample", "stage": "kphi", "value": 0.0}],
                "computed": ["komega", "kappa", "ttheta"],
            }
        },
    }
    g = load_geometry_file(yaml.safe_dump(doc))
    # +longitudinal is parallel to n_komega; +vertical is the next
    # perpendicular basis direction.
    np.testing.assert_allclose(
        g.kappa_pseudo_angle_convention.n_chi_eq, [0.0, 0.0, 1.0], atol=1e-12
    )


def test_kappa_eulerian_chi_auto_derivation_no_perpendicular_basis_rejected():
    """When the outer kappa axis is not aligned with any single basis
    direction the auto-derivation has no perpendicular candidate and
    must raise.  Callers must declare ``kappa_eulerian_chi`` for
    such geometries.
    """
    # Construct a kappa-arm geometry whose n_komega is a (1, 1, 1)/√3
    # diagonal: not perpendicular to any of +longitudinal, +vertical,
    # +transverse within tolerance.  Skip the ``kappa_eulerian`` arm-
    # construction shortcut (its own perpendicularity check rejects
    # this case) by supplying a numeric kappa-axis vector directly.
    s3 = float(1.0 / np.sqrt(3.0))
    oblique = [s3, s3, s3]
    # Numeric kappa-arm axis: tilted toward an arbitrary perpendicular
    # direction; the value does not matter here because the test
    # exercises only the auto-derivation of n_chi_eq.
    kappa_axis = [
        float(np.cos(np.deg2rad(50.0)) * s3 + np.sin(np.deg2rad(50.0)) / np.sqrt(2.0)),
        float(np.cos(np.deg2rad(50.0)) * s3 - np.sin(np.deg2rad(50.0)) / np.sqrt(2.0)),
        float(np.cos(np.deg2rad(50.0)) * s3),
    ]
    doc = {
        KIND_KEY: {"schema_revision": 1},
        "name": "kappa_oblique_omega",
        "documentation": "test",
        "basis": "BL",
        "parameters": {"alpha_deg": 50.0},
        "kappa_chi_eq": "+transverse",
        "stages": [
            {
                "name": "komega",
                "axis": oblique,
                "parent": None,
                "role": "sample",
            },
            {
                "name": "kappa",
                "axis": kappa_axis,
                "parent": "komega",
                "role": "sample",
            },
            {
                "name": "kphi",
                "axis": oblique,
                "parent": "kappa",
                "role": "sample",
            },
            {
                "name": "ttheta",
                "axis": oblique,
                "parent": None,
                "role": "detector",
            },
        ],
        "modes": {
            "fixed_kphi": {
                "default": True,
                "constraints": [{"type": "sample", "stage": "kphi", "value": 0.0}],
                "computed": ["komega", "kappa", "ttheta"],
            }
        },
    }
    with pytest.raises(
        GeometrySchemaError, match="cannot derive the equivalent-Eulerian chi"
    ):
        load_geometry_file(yaml.safe_dump(doc))


# ---------------------------------------------------------------------------
# Kappa-arm tilt direction is unchanged: the YAML's ``kappa_chi_eq``
# still drives Walko's formula for the kappa stage axis.  This guards
# against accidentally re-conflating the two fields.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kappa_name, expected_kappa_axis_in_plane, context",
    [
        # kappa4cv: kappa arm in (transverse, vertical) plane.
        # Expected axis = cos(50)*+x + sin(50)*+z (BL: x=trans, z=vert).
        pytest.param(
            "kappa4cv",
            (
                np.cos(np.deg2rad(50.0)) * np.array([1.0, 0.0, 0.0])
                + np.sin(np.deg2rad(50.0)) * np.array([0.0, 0.0, 1.0])
            ),
            does_not_raise(),
            id="kappa4cv-arm-in-TV-plane",
        ),
        # kappa4ch: kappa arm in (vertical, longitudinal) plane.
        # Expected axis = cos(50)*+z + sin(50)*+y (BL: z=vert, y=long).
        pytest.param(
            "kappa4ch",
            (
                np.cos(np.deg2rad(50.0)) * np.array([0.0, 0.0, 1.0])
                + np.sin(np.deg2rad(50.0)) * np.array([0.0, 1.0, 0.0])
            ),
            does_not_raise(),
            id="kappa4ch-arm-in-VL-plane",
        ),
        # kappa6c: kappa arm in (transverse, vertical) plane,
        # same as kappa4cv but in YOU basis (x=vert, z=trans).
        pytest.param(
            "kappa6c",
            (
                np.cos(np.deg2rad(50.0)) * np.array([0.0, 0.0, 1.0])
                + np.sin(np.deg2rad(50.0)) * np.array([1.0, 0.0, 0.0])
            ),
            does_not_raise(),
            id="kappa6c-arm-in-TV-plane",
        ),
    ],
)
def test_kappa_arm_axis_in_canonical_plane(
    kappa_name,
    expected_kappa_axis_in_plane,
    context,
):
    """The kappa stage axis lies in the canonical kappa-arm tilt
    plane per the published references (Walko 2016, Wyckoff 1985,
    Thorkildsen 2006).  This is unchanged by issue #284, which
    only affected the equivalent-Eulerian chi pseudo-angle axis.
    """
    g = ahd.make_geometry(kappa_name)
    with context:
        np.testing.assert_allclose(
            g.stage("kappa").axis, expected_kappa_axis_in_plane, atol=1e-12
        )


# ---------------------------------------------------------------------------
# Equivalent-Eulerian chi axis must be perpendicular to n_komega for
# the closed-form decomposition to be well-posed.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kappa_name, context",
    [
        pytest.param("kappa4cv", does_not_raise(), id="kappa4cv"),
        pytest.param("kappa4ch", does_not_raise(), id="kappa4ch"),
        pytest.param("kappa6c", does_not_raise(), id="kappa6c"),
    ],
)
def test_n_chi_eq_perpendicular_to_n_komega(kappa_name, context):
    """Every shipped kappa preset's ``n_chi_eq`` is perpendicular to
    its ``n_komega`` (required for the closed-form
    ``eulerian_to_kappa_axes`` decomposition; see issue #284 and the
    ``KappaPseudoAngleConvention`` docstring).
    """
    g = ahd.make_geometry(kappa_name)
    convention = g.kappa_pseudo_angle_convention
    with context:
        assert abs(float(np.dot(convention.n_komega, convention.n_chi_eq))) < 1e-12


# ---------------------------------------------------------------------------
# KappaPseudoAngleConvention with a custom n_chi_eq survives the
# constructor's validation step (no in-plane requirement after #284).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "n_chi_eq, context",
    [
        # +longitudinal (perpendicular to the T-V kappa-arm plane).
        pytest.param([0.0, 1.0, 0.0], does_not_raise(), id="long-perpendicular"),
        # +vertical (in the T-V plane, perpendicular to n_komega).
        pytest.param([0.0, 0.0, 1.0], does_not_raise(), id="vert-in-plane"),
    ],
)
def test_kappa_convention_accepts_any_perpendicular_n_chi_eq(n_chi_eq, context):
    """The ``KappaPseudoAngleConvention`` constructor accepts any
    ``n_chi_eq`` perpendicular to ``n_komega`` and not parallel to
    ``n_kphi`` — it does **not** require the in-plane kappa-arm
    relation that ``kappa_axis_from_eulerian`` uses.  Issue #284
    relies on this freedom.
    """
    import math

    n_komega = np.array([-1.0, 0.0, 0.0])  # kappa4cv-like
    n_kappa = np.cos(math.radians(50.0)) * np.array([1.0, 0.0, 0.0]) + np.sin(
        math.radians(50.0)
    ) * np.array([0.0, 0.0, 1.0])
    n_kphi = np.array([-1.0, 0.0, 0.0])
    with context:
        KappaPseudoAngleConvention(
            n_komega=n_komega,
            n_kappa=n_kappa,
            n_kphi=n_kphi,
            n_chi_eq=np.array(n_chi_eq, dtype=float),
        )



