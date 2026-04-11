# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
constants.py — shared Cartesian basis vectors.

These are the internal numpy representations used throughout the package.
The caller-facing notation (+x, -z, etc.) is handled in axes.py.

Convention follows You (1999):
    XHAT = vertical (out of floor)
    YHAT = longitudinal (along beam / toward equipment)
    ZHAT = lateral (to our left)
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)

logger = logging.getLogger(__name__)

XHAT = np.array([1.0, 0.0, 0.0])
YHAT = np.array([0.0, 1.0, 0.0])
ZHAT = np.array([0.0, 0.0, 1.0])
