# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.drawing.

Covers:
  - _physical_label() helper function
  - arc direction handedness for all preset geometry stages
"""

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.drawing import _physical_label
from ad_hoc_diffractometer.factories import BASIS_BL
from ad_hoc_diffractometer.factories import BASIS_YOU

# ---------------------------------------------------------------------------
# _physical_label
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "axis_vec, basis, expected, context",
    [
        pytest.param(
            np.array([1.0, 0.0, 0.0]),
            BASIS_BL,
            "+transverse",
            does_not_raise(),
            id="positive-transverse-BL",
        ),
        pytest.param(
            np.array([-1.0, 0.0, 0.0]),
            BASIS_BL,
            "-transverse",
            does_not_raise(),
            id="negative-transverse-BL",
        ),
        pytest.param(
            np.array([0.0, 0.0, 1.0]),
            BASIS_BL,
            "+vertical",
            does_not_raise(),
            id="positive-vertical-BL",
        ),
        pytest.param(
            np.array([1.0, 0.0, 0.0]),
            BASIS_YOU,
            "+vertical",
            does_not_raise(),
            id="positive-vertical-YOU",
        ),
        pytest.param(
            np.array([0.5, 0.0, 0.5]),
            BASIS_BL,
            "[0.5, 0, 0.5]",  # falls back to axis_label for non-standard vectors
            does_not_raise(),
            id="non-standard-axis",
        ),
        pytest.param(
            np.array([1.0, 0.0, 0.0]),
            None,
            "+x",
            does_not_raise(),
            id="no-basis",
        ),
    ],
)
def test_physical_label(axis_vec, basis, expected, context):
    """_physical_label returns expected direction labels."""
    with context:
        result = _physical_label(axis_vec, basis)
        assert result == expected


# ---------------------------------------------------------------------------
# Arc direction handedness — reproduces the math from _draw_axis and
# _draw_axis_arc to verify the perpendicular frame has the correct
# handedness for every stage of every preset geometry.
#
# The invariant:
#   direction=+1 (right-handed stage) => dot(cross(perp1, perp2), axis_d) > 0
#   direction=-1 (left-handed stage)  => dot(cross(perp1, perp2), axis_d) < 0
#
# This ensures a sweep from -π/3 to +π/3 draws the arc in the
# physically correct rotation sense.  See issue #207.
# ---------------------------------------------------------------------------


def _build_arc_params():
    """Generate (geometry_name, stage_name) pairs for all presets."""
    params = []
    for geom_name in ahd.list_geometries():
        geometry = ahd.make_geometry(geom_name)
        for stage_name in geometry._stages:
            params.append(
                pytest.param(
                    geom_name,
                    stage_name,
                    does_not_raise(),
                    id=f"{geom_name}-{stage_name}",
                )
            )
    return params


@pytest.mark.parametrize("geometry_name, stage_name, context", _build_arc_params())
def test_arc_direction_handedness(geometry_name, stage_name, context):
    """Perpendicular frame handedness matches stage rotation direction.

    Reproduces the exact math from _draw_axis / _draw_axis_arc (without
    Plotly) and checks the handedness invariant for every preset stage.
    """
    with context:
        geometry = ahd.make_geometry(geometry_name)
        stage = geometry._stages[stage_name]
        basis = geometry.basis

        # -- _draw_basis_vectors: compute display rotation R --
        ver = np.asarray(basis["vertical"], dtype=float)
        lon = np.asarray(basis["longitudinal"], dtype=float)
        lat = np.asarray(basis["transverse"], dtype=float)
        R = np.array([lat, -lon, ver])

        # -- _draw_axis_vector: normalize and flip to positive --
        stage_axis = np.asarray(stage.axis, dtype=float)
        axis_norm = stage_axis / np.linalg.norm(stage_axis)
        for component in axis_norm:
            if abs(component) > 1e-9:
                if component < 0:
                    axis_norm = -axis_norm
                break

        # -- _draw_axis: compute direction --
        stage_norm = stage_axis / np.linalg.norm(stage_axis)
        direction = +1 if np.dot(stage_norm, axis_norm) > 0 else -1

        # -- _draw_axis_arc: build perpendicular frame --
        axis_d = R @ axis_norm

        eye = np.array([0.92, 1.18, -0.10])
        view_dir = -eye / np.linalg.norm(eye)
        cam_up = np.array([0.0, 0.0, 1.0])
        screen_right = np.cross(cam_up, view_dir)
        screen_right /= np.linalg.norm(screen_right)
        screen_up = np.cross(view_dir, screen_right)

        perp1_d = screen_right - np.dot(screen_right, axis_d) * axis_d
        if np.linalg.norm(perp1_d) < 1e-6:
            perp1_d = screen_up - np.dot(screen_up, axis_d) * axis_d
        perp1_d /= np.linalg.norm(perp1_d)

        perp2_d = np.cross(axis_d, perp1_d)
        if direction < 0:
            perp2_d = -perp2_d

        # -- invariant: frame handedness must match direction --
        handedness = np.dot(np.cross(perp1_d, perp2_d), axis_d)
        assert handedness == pytest.approx(float(direction), abs=1e-10), (
            f"{geometry_name}/{stage_name}: expected {direction}, got {handedness}"
        )
