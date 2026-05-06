# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Basis-invariance tests for all preset geometries.

The choice of internal Cartesian basis (BL vs You) is an
implementation detail and must not be visible in any user-facing
physical quantity.  Switching a preset's ``basis`` argument from
BL to You (or vice versa) produces an internal representation
where each Cartesian axis carries a different physical-direction
label, but every basis-invariant quantity (``|Q|``, ``|2θ|``,
``d``-spacing, …) must come out the same.

These tests guard against any future leak of internal-basis
choice into user-visible behavior, including from work in issue
#252 that re-derived the kappa-axis convention.
"""

from __future__ import annotations

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.factories import BASIS_BL
from ad_hoc_diffractometer.factories import BASIS_YOU

# Preset registry — all 10 geometries.  The detector-stage name (last
# detector stage) and the default mode used for the |2θ| comparison
# are listed alongside each preset.  Presets whose default mode is
# not yet implemented for ``forward()`` are flagged with
# ``forward_implemented=False`` and exercise only the |Q| invariance.

_PRESETS = [
    ("fourcv", "ttheta", None, True),
    ("fourch", "ttheta", None, True),
    ("psic", "delta", None, True),
    ("sixc", "delta", None, True),
    ("kappa4cv", "ttheta", None, True),
    ("kappa4ch", "ttheta", None, True),
    ("kappa6c", "delta", None, True),
    ("zaxis", "delta", "zaxis", False),
    ("s2d2", "delta", "fixed_mu", False),
    ("fivec", "ttheta", None, True),
]


def _setup(preset_name, basis, mode_name=None, a=4.0, wavelength=1.0):
    """Build a preset under the given basis with a cubic UB=identity."""
    g = ahd.make_geometry(preset_name, basis=basis)
    g.wavelength = wavelength
    g.sample.lattice = ahd.Lattice(a=a)
    ahd.ub_identity(g.sample)
    if mode_name is not None:
        g.mode_name = mode_name
    return g


# ---------------------------------------------------------------------------
# |Q| invariance — must hold for every preset and every (h, k, l)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "preset_name",
    [pytest.param(p[0], id=p[0]) for p in _PRESETS],
)
@pytest.mark.parametrize(
    "h, k, l, context",
    [
        pytest.param(1, 0, 0, does_not_raise(), id="100"),
        pytest.param(0, 1, 0, does_not_raise(), id="010"),
        pytest.param(0, 0, 1, does_not_raise(), id="001"),
        pytest.param(1, 1, 0, does_not_raise(), id="110"),
        pytest.param(1, 1, 1, does_not_raise(), id="111"),
        pytest.param(2, 1, 0, does_not_raise(), id="210"),
    ],
)
def test_q_magnitude_is_basis_invariant(preset_name, h, k, l, context):  # noqa: E741
    """``|UB · hkl|`` is identical when the preset is constructed with
    BL vs You basis (same lattice, same ``ub_identity``)."""
    with context:
        g_bl = _setup(preset_name, BASIS_BL)
        g_yo = _setup(preset_name, BASIS_YOU)
        hkl = np.array([h, k, l], dtype=float)
        q_bl = float(np.linalg.norm(g_bl.sample.UB @ hkl))
        q_yo = float(np.linalg.norm(g_yo.sample.UB @ hkl))
        assert q_bl == pytest.approx(q_yo, abs=1e-12), (
            f"{preset_name}: |Q| differs between BL ({q_bl:.6e}) and "
            f"You ({q_yo:.6e}) for ({h},{k},{l})"
        )


# ---------------------------------------------------------------------------
# |2θ| invariance — for presets whose default mode supports forward()
# ---------------------------------------------------------------------------


_FWD_PRESETS = [(p[0], p[1], p[2]) for p in _PRESETS if p[3]]


@pytest.mark.parametrize(
    "preset_name, det_name, mode_name",
    [pytest.param(*p, id=p[0]) for p in _FWD_PRESETS],
)
@pytest.mark.parametrize(
    "h, k, l, context",
    [
        pytest.param(1, 0, 0, does_not_raise(), id="100"),
        pytest.param(0, 1, 0, does_not_raise(), id="010"),
        pytest.param(0, 0, 1, does_not_raise(), id="001"),
    ],
)
def test_two_theta_is_basis_invariant(
    preset_name,
    det_name,
    mode_name,
    h,
    k,
    l,  # noqa: E741
    context,
):
    """The detector angle (``|2θ|`` or ``|delta|``) returned by
    ``forward(h, k, l)`` is identical when the preset is constructed
    with BL vs You basis.

    Reflections that are unreachable in either basis are skipped
    (a reflection that fails to solve in *both* bases is not a
    basis-invariance failure; a reflection that solves in one but
    not the other *is* a failure and is asserted explicitly).
    """
    with context:
        g_bl = _setup(preset_name, BASIS_BL, mode_name=mode_name)
        g_yo = _setup(preset_name, BASIS_YOU, mode_name=mode_name)
        sols_bl = g_bl.forward(h, k, l)
        sols_yo = g_yo.forward(h, k, l)
        if not sols_bl and not sols_yo:
            pytest.skip(f"{preset_name}: ({h},{k},{l}) unreachable in both bases")
        assert bool(sols_bl) == bool(sols_yo), (
            f"{preset_name}: ({h},{k},{l}) reachable in one basis but "
            f"not the other — basis invariance violated"
        )
        tt_bl = abs(sols_bl[0][det_name])
        tt_yo = abs(sols_yo[0][det_name])
        assert tt_bl == pytest.approx(tt_yo, abs=1e-9), (
            f"{preset_name}: |{det_name}| differs between BL ({tt_bl:.6e}) "
            f"and You ({tt_yo:.6e}) for ({h},{k},{l})"
        )
