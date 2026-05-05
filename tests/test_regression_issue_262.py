# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Cross-module regression tests for issue #262 — psic and kappa6c zone modes.

Zone mode (You 1999 §6, SPEC ``setmode 5``) confines the scattering vector
Q to the plane spanned by two reciprocal-lattice vectors ``z0`` and
``z1``.  The implementation is "extras-driven", parallel to the existing
``double_diffraction`` modes:

- The mode definition lives in :mod:`presets` (zone_vertical and
  zone_horizontal on both psic and kappa6c).
- The dispatcher :func:`forward._is_zone_mode` recognises modes whose
  ``extras`` carry both ``z0`` and ``z1`` keys.
- :func:`forward._solve_zone` validates inputs, computes the zone-plane
  normal in the φ frame, applies the in-plane prefilter, and delegates
  to the existing bisecting solver (which itself dispatches to the
  kappa-virtual variant for kappa6c).

This file lives at the cross-module level because the change spans
``presets.py`` (mode registration), ``forward.py`` (dispatcher and
solver), ``mode.py`` (extras-driven dispatch documentation), and
``benchmark.py`` (default ``z0``/``z1`` for the slow-benchmark sweep).
"""

from __future__ import annotations

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.forward import _is_zone_mode

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _setup_cubic(factory, a=4.0, wavelength=1.5406):
    """Return a fresh diffractometer with a cubic lattice and UB = identity."""
    g = factory()
    g.wavelength = wavelength
    g.sample.lattice = ahd.Lattice(a=a)
    ahd.ub_identity(g.sample)
    return g


def _activate_zone(geom, mode_name, z0, z1):
    """Activate ``mode_name`` and set the zone-plane vectors."""
    geom.mode_name = mode_name
    cs = geom.modes[mode_name]
    cs.extras["z0"] = z0
    cs.extras["z1"] = z1
    return cs


# Zone-plane combinations that work for each (geometry, mode_name) pair.
# psic horizontal cannot reach the (h, k, 0) plane in the You basis
# (vertical = +x), so its tests use the (h, 0, l) plane instead.
_GEOM_MODES = [
    pytest.param(
        ahd.presets.psic,
        "zone_vertical",
        (1, 0, 0),
        (0, 1, 0),
        id="psic-vertical-hk0",
    ),
    pytest.param(
        ahd.presets.psic,
        "zone_horizontal",
        (1, 0, 0),
        (0, 0, 1),
        id="psic-horizontal-h0l",
    ),
    pytest.param(
        ahd.presets.kappa6c,
        "zone_vertical",
        (1, 0, 0),
        (0, 1, 0),
        id="kappa6c-vertical-hk0",
    ),
    pytest.param(
        ahd.presets.kappa6c,
        "zone_horizontal",
        (1, 0, 0),
        (0, 1, 0),
        id="kappa6c-horizontal-hk0",
    ),
]


# ---------------------------------------------------------------------------
# Mode registration
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, mode_name, context",
    [
        pytest.param(
            ahd.presets.psic,
            "zone_vertical",
            does_not_raise(),
            id="psic-zone_vertical-registered",
        ),
        pytest.param(
            ahd.presets.psic,
            "zone_horizontal",
            does_not_raise(),
            id="psic-zone_horizontal-registered",
        ),
        pytest.param(
            ahd.presets.kappa6c,
            "zone_vertical",
            does_not_raise(),
            id="kappa6c-zone_vertical-registered",
        ),
        pytest.param(
            ahd.presets.kappa6c,
            "zone_horizontal",
            does_not_raise(),
            id="kappa6c-zone_horizontal-registered",
        ),
    ],
)
def test_zone_mode_registered(factory, mode_name, context):
    """Both zone modes are registered on psic and kappa6c, with the
    expected extras schema (``z0`` and ``z1`` REQUIRED, ``in_plane_residual``
    None) and dispatcher recognition.
    """
    with context:
        from ad_hoc_diffractometer.mode import REQUIRED

        g = factory()
        cs = g.modes[mode_name]
        assert cs.extras["z0"] is REQUIRED
        assert cs.extras["z1"] is REQUIRED
        assert cs.extras["in_plane_residual"] is None
        assert _is_zone_mode(g, cs)
        assert cs.is_implemented(g)


# ---------------------------------------------------------------------------
# Extras validation: missing, malformed, zero, parallel
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "z0, z1, context",
    [
        # z0/z1 both still REQUIRED → ValueError mentioning the extras keys
        pytest.param(
            None,
            None,
            pytest.raises(
                ValueError,
                match=re.escape("zone mode requires z0 and z1"),
            ),
            id="both-required",
        ),
        # z0 zero vector
        pytest.param(
            (0, 0, 0),
            (0, 1, 0),
            pytest.raises(
                ValueError,
                match=re.escape("non-zero reciprocal-lattice vectors"),
            ),
            id="z0-zero",
        ),
        # z1 zero vector
        pytest.param(
            (1, 0, 0),
            (0, 0, 0),
            pytest.raises(
                ValueError,
                match=re.escape("non-zero reciprocal-lattice vectors"),
            ),
            id="z1-zero",
        ),
        # parallel z0 and z1 (collinear) → undefined plane normal
        pytest.param(
            (1, 0, 0),
            (2, 0, 0),
            pytest.raises(
                ValueError,
                match=re.escape("parallel"),
            ),
            id="parallel",
        ),
        # malformed: not a 3-vector
        pytest.param(
            (1, 0),
            (0, 1, 0),
            pytest.raises(
                ValueError,
                match=re.escape("3-element sequence"),
            ),
            id="z0-too-short",
        ),
    ],
)
def test_zone_extras_validation(z0, z1, context):
    """The zone solver rejects malformed, zero, or parallel z0/z1
    inputs with a descriptive ``ValueError``.
    """
    g = _setup_cubic(ahd.presets.psic)
    g.mode_name = "zone_vertical"
    if z0 is not None:
        g.modes["zone_vertical"].extras["z0"] = z0
    if z1 is not None:
        g.modes["zone_vertical"].extras["z1"] = z1
    with context:
        g.forward(1, 0, 0)


# ---------------------------------------------------------------------------
# In-plane reflections: forward returns ≥ 1 solution; round-trip closes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("factory, mode_name, z0, z1", _GEOM_MODES)
@pytest.mark.parametrize(
    "hkl",
    [
        pytest.param((1, 0, 0), id="100"),
        pytest.param((0, 1, 0), id="010"),
        pytest.param((1, 1, 0), id="110"),
        pytest.param((2, 0, 0), id="200"),
    ],
)
def test_zone_in_plane_round_trip(factory, mode_name, z0, z1, hkl):
    """Reflections that lie in the zone plane:

    1. Return at least one solution (when geometrically reachable).
    2. Round-trip ``forward`` → ``inverse`` to within machine precision.
    3. Have ``in_plane_residual`` recorded near zero in mode.extras.

    A reflection that is in-plane mathematically but unreachable in the
    chosen scattering geometry (e.g. ``(0,1,0)`` in psic's horizontal
    family using the ``(h,0,l)`` zone plane) is skipped, not failed —
    that is a property of the geometry, not the zone solver.
    """
    g = _setup_cubic(factory)

    # Skip hkl values not in the active zone plane (some test cases
    # cross-combine).  Compute Q-component along z-normal in the
    # cubic-orthonormal reciprocal frame.
    z0_arr = np.asarray(z0, dtype=float)
    z1_arr = np.asarray(z1, dtype=float)
    n_zone = np.cross(z0_arr, z1_arr)
    if abs(float(np.dot(np.asarray(hkl, dtype=float), n_zone))) > 1e-9:
        pytest.skip(f"{hkl} is not in the zone plane spanned by {z0} and {z1}")

    cs = _activate_zone(g, mode_name, z0, z1)
    sols = g.forward(*hkl)
    if not sols:
        pytest.skip(
            f"{factory.__name__}/{mode_name}: {hkl} is in the zone plane "
            f"but unreachable in this scattering geometry."
        )
    assert cs.extras["in_plane_residual"] == pytest.approx(0.0, abs=1e-10)
    for sol in sols:
        np.testing.assert_allclose(g.inverse(sol), hkl, atol=1e-8)


# ---------------------------------------------------------------------------
# Off-plane reflections: empty list + warning
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, mode_name, z0, z1, off_plane_hkl",
    [
        pytest.param(
            ahd.presets.psic,
            "zone_vertical",
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            id="psic-vertical-off-001",
        ),
        pytest.param(
            ahd.presets.psic,
            "zone_horizontal",
            (1, 0, 0),
            (0, 0, 1),
            (0, 1, 0),
            id="psic-horizontal-off-010",
        ),
        pytest.param(
            ahd.presets.kappa6c,
            "zone_vertical",
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            id="kappa6c-vertical-off-001",
        ),
        pytest.param(
            ahd.presets.kappa6c,
            "zone_horizontal",
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            id="kappa6c-horizontal-off-001",
        ),
    ],
)
def test_zone_off_plane_returns_empty_with_warning(
    factory, mode_name, z0, z1, off_plane_hkl, caplog
):
    """An (h, k, l) that does not lie in the zone plane returns an empty
    solution list and logs a warning identifying the in-plane residual
    (matching the soft-failure pattern used for ``fixed_psi`` in #176).
    """
    import logging

    g = _setup_cubic(factory)
    cs = _activate_zone(g, mode_name, z0, z1)

    with caplog.at_level(logging.WARNING, logger="ad_hoc_diffractometer.forward"):
        sols = g.forward(*off_plane_hkl)

    assert sols == []
    assert cs.extras["in_plane_residual"] > 1e-6
    assert any(
        "zone mode" in record.message and "not in the zone plane" in record.message
        for record in caplog.records
    ), f"expected zone-mode warning; got: {[r.message for r in caplog.records]}"


# ---------------------------------------------------------------------------
# Plane-lock invariants: the constant_stages of each mode are honored
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, mode_name, z0, z1, in_plane_hkl, constant_stage_values",
    [
        pytest.param(
            ahd.presets.psic,
            "zone_vertical",
            (1, 0, 0),
            (0, 1, 0),
            (1, 1, 0),
            {"mu": 0.0, "nu": 0.0},
            id="psic-vertical-mu0-nu0",
        ),
        pytest.param(
            ahd.presets.psic,
            "zone_horizontal",
            (1, 0, 0),
            (0, 0, 1),
            (1, 0, 1),
            {"eta": 0.0, "delta": 0.0},
            id="psic-horizontal-eta0-delta0",
        ),
        pytest.param(
            ahd.presets.kappa6c,
            "zone_vertical",
            (1, 0, 0),
            (0, 1, 0),
            (1, 1, 0),
            {"mu": 0.0, "nu": 0.0},
            id="kappa6c-vertical-mu0-nu0",
        ),
        pytest.param(
            ahd.presets.kappa6c,
            "zone_horizontal",
            (1, 0, 0),
            (0, 1, 0),
            (1, 1, 0),
            {"komega": 0.0, "delta": 0.0},
            id="kappa6c-horizontal-komega0-delta0",
        ),
    ],
)
def test_zone_constant_stages(
    factory, mode_name, z0, z1, in_plane_hkl, constant_stage_values
):
    """Every solution has the constant stages frozen at their declared
    values (mu/nu for vertical, eta/delta for psic horizontal,
    komega/delta for kappa6c horizontal).  This is the structural
    counterpart of the in-plane prefilter — the bisecting solver may
    not violate the SampleConstraint/DetectorConstraint values that
    accompany the zone constraint.
    """
    g = _setup_cubic(factory)
    _activate_zone(g, mode_name, z0, z1)
    sols = g.forward(*in_plane_hkl)
    if not sols:  # pragma: no cover — covered by round-trip test above
        pytest.skip(f"{factory.__name__}/{mode_name}: {in_plane_hkl} unreachable.")
    for sol in sols:
        for stage, expected in constant_stage_values.items():
            assert sol[stage] == pytest.approx(expected, abs=1e-10), (
                f"{factory.__name__}/{mode_name}: stage {stage!r} should "
                f"be {expected}, got {sol[stage]} for hkl={in_plane_hkl}."
            )
