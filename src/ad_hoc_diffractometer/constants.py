"""
constants.py — shared Cartesian basis vectors.

These are the internal numpy representations used throughout the package.
The caller-facing notation (+x, -z, etc.) is handled in axes.py.

Convention follows You (1999):
    XHAT = vertical (out of floor)
    YHAT = longitudinal (along beam / toward equipment)
    ZHAT = lateral (to our left)
"""

import numpy as np

XHAT = np.array([1.0, 0.0, 0.0])
YHAT = np.array([0.0, 1.0, 0.0])
ZHAT = np.array([0.0, 0.0, 1.0])
