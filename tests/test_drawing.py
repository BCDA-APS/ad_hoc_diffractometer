# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.drawing.

Covers:
  - _physical_label() helper function
"""

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

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
