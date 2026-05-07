# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
Unit tests for ad_hoc_diffractometer.rotation.

Covers:
  - rotation_matrix()
"""

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from helpers import Rx
from helpers import Ry
from helpers import Rz

from ad_hoc_diffractometer.constants import XHAT
from ad_hoc_diffractometer.constants import YHAT
from ad_hoc_diffractometer.constants import ZHAT
from ad_hoc_diffractometer.rotation import _rotation_matrix_and_derivative_normalized
from ad_hoc_diffractometer.rotation import rotation_matrix

# ---------------------------------------------------------------------------
# rotation_matrix() — correctness
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "axis, angle_deg, expected, context",
    [
        pytest.param(
            XHAT, 0.0, np.eye(3), does_not_raise(), id="zero-angle-x-is-identity"
        ),
        pytest.param(
            YHAT, 0.0, np.eye(3), does_not_raise(), id="zero-angle-y-is-identity"
        ),
        pytest.param(
            ZHAT, 0.0, np.eye(3), does_not_raise(), id="zero-angle-z-is-identity"
        ),
        pytest.param(XHAT, 90.0, Rx(90), does_not_raise(), id="90deg-about-x"),
        pytest.param(YHAT, 90.0, Ry(90), does_not_raise(), id="90deg-about-y"),
        pytest.param(ZHAT, 90.0, Rz(90), does_not_raise(), id="90deg-about-z"),
        pytest.param(XHAT, 180.0, Rx(180), does_not_raise(), id="180deg-about-x"),
        pytest.param(YHAT, 180.0, Ry(180), does_not_raise(), id="180deg-about-y"),
        pytest.param(ZHAT, 180.0, Rz(180), does_not_raise(), id="180deg-about-z"),
        pytest.param(
            XHAT, -45.0, Rx(-45), does_not_raise(), id="negative-angle-about-x"
        ),
        pytest.param(
            ZHAT, -90.0, Rz(-90), does_not_raise(), id="negative-angle-about-z"
        ),
        pytest.param(
            -ZHAT,
            30.0,
            Rz(-30),
            does_not_raise(),
            id="negated-axis-equals-negated-angle",
        ),
        pytest.param(
            -XHAT,
            45.0,
            Rx(-45),
            does_not_raise(),
            id="negated-x-axis-equals-negated-angle",
        ),
        pytest.param(
            2 * XHAT, 90.0, Rx(90), does_not_raise(), id="non-unit-axis-is-normalized"
        ),
        pytest.param(YHAT, 360.0, np.eye(3), does_not_raise(), id="360deg-is-identity"),
    ],
)
def test_rotation_matrix(axis, angle_deg, expected, context):
    with context:
        R = rotation_matrix(axis, angle_deg)
        assert R.shape == (3, 3)
        np.testing.assert_allclose(R, expected, atol=1e-10)


# ---------------------------------------------------------------------------
# rotation_matrix() — orthogonality and determinant
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "axis, angle_deg, context",
    [
        pytest.param(XHAT, 37.0, does_not_raise(), id="x-axis-37deg"),
        pytest.param(YHAT, -73.0, does_not_raise(), id="y-axis-neg73deg"),
        pytest.param(ZHAT, 120.0, does_not_raise(), id="z-axis-120deg"),
        pytest.param(-ZHAT, 55.0, does_not_raise(), id="neg-z-axis-55deg"),
        pytest.param(
            np.array([1, 1, 1]), 60.0, does_not_raise(), id="diagonal-axis-60deg"
        ),
    ],
)
def test_rotation_matrix_orthogonal(axis, angle_deg, context):
    with context:
        R = rotation_matrix(axis, angle_deg)
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-10)
        assert abs(np.linalg.det(R) - 1.0) < 1e-10


# ---------------------------------------------------------------------------
# _rotation_matrix_and_derivative_normalized() — R matches, dR is correct
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "axis, angle_deg, context",
    [
        pytest.param(XHAT, 0.0, does_not_raise(), id="x-0deg"),
        pytest.param(XHAT, 90.0, does_not_raise(), id="x-90deg"),
        pytest.param(YHAT, 45.0, does_not_raise(), id="y-45deg"),
        pytest.param(ZHAT, 180.0, does_not_raise(), id="z-180deg"),
        pytest.param(-ZHAT, 30.0, does_not_raise(), id="neg-z-30deg"),
        pytest.param(XHAT, -60.0, does_not_raise(), id="x-neg60deg"),
        pytest.param(
            np.array([1, 1, 1]) / np.sqrt(3),
            72.0,
            does_not_raise(),
            id="diagonal-72deg",
        ),
        pytest.param(YHAT, 360.0, does_not_raise(), id="y-360deg"),
    ],
)
def test_rotation_matrix_and_derivative_R_matches(axis, angle_deg, context):
    """The R returned by the combined function matches rotation_matrix()."""
    with context:
        n = np.asarray(axis, dtype=float)
        n = n / np.linalg.norm(n)
        R, _dR = _rotation_matrix_and_derivative_normalized(n, angle_deg)
        R_expected = rotation_matrix(axis, angle_deg)
        np.testing.assert_allclose(R, R_expected, atol=1e-14)


@pytest.mark.parametrize(
    "axis, angle_deg, context",
    [
        pytest.param(XHAT, 0.0, does_not_raise(), id="x-0deg"),
        pytest.param(XHAT, 90.0, does_not_raise(), id="x-90deg"),
        pytest.param(YHAT, 45.0, does_not_raise(), id="y-45deg"),
        pytest.param(ZHAT, 180.0, does_not_raise(), id="z-180deg"),
        pytest.param(-ZHAT, 30.0, does_not_raise(), id="neg-z-30deg"),
        pytest.param(XHAT, -60.0, does_not_raise(), id="x-neg60deg"),
        pytest.param(
            np.array([1, 1, 1]) / np.sqrt(3),
            72.0,
            does_not_raise(),
            id="diagonal-72deg",
        ),
        pytest.param(YHAT, 360.0, does_not_raise(), id="y-360deg"),
    ],
)
def test_rotation_matrix_and_derivative_dR_vs_finite_difference(
    axis, angle_deg, context
):
    """dR/dtheta agrees with a central finite-difference estimate."""
    with context:
        n = np.asarray(axis, dtype=float)
        n = n / np.linalg.norm(n)
        _R, dR = _rotation_matrix_and_derivative_normalized(n, angle_deg)

        # Central finite-difference in radians
        h_deg = 1e-6
        R_plus = rotation_matrix(n, angle_deg + h_deg)
        R_minus = rotation_matrix(n, angle_deg - h_deg)
        h_rad = np.deg2rad(h_deg)
        dR_fd = (R_plus - R_minus) / (2 * h_rad)

        np.testing.assert_allclose(dR, dR_fd, atol=1e-6)


@pytest.mark.parametrize(
    "axis, angle_deg, context",
    [
        pytest.param(XHAT, 37.0, does_not_raise(), id="x-37deg"),
        pytest.param(YHAT, -73.0, does_not_raise(), id="y-neg73deg"),
        pytest.param(ZHAT, 120.0, does_not_raise(), id="z-120deg"),
        pytest.param(-ZHAT, 55.0, does_not_raise(), id="neg-z-55deg"),
        pytest.param(
            np.array([1, 1, 1]) / np.sqrt(3),
            60.0,
            does_not_raise(),
            id="diagonal-60deg",
        ),
    ],
)
def test_rotation_matrix_and_derivative_dR_antisymmetric_product(
    axis, angle_deg, context
):
    """R^T @ dR/dtheta is antisymmetric (a rotation derivative property)."""
    with context:
        n = np.asarray(axis, dtype=float)
        n = n / np.linalg.norm(n)
        R, dR = _rotation_matrix_and_derivative_normalized(n, angle_deg)
        product = R.T @ dR
        # R^T dR should be antisymmetric: product + product.T == 0
        np.testing.assert_allclose(product + product.T, np.zeros((3, 3)), atol=1e-12)
