"""
helpers.py — shared test utilities for the test suite.

Pure functions and constants used in parametrize() decorators, where
pytest fixtures cannot be used (collection-time evaluation).

Import explicitly in any test file that needs them:
    from helpers import Rx, Ry, Rz, STANDARD_BASIS
"""

import numpy as np

from ad_hoc_diffractometer.constants import XHAT
from ad_hoc_diffractometer.constants import YHAT
from ad_hoc_diffractometer.constants import ZHAT

# ---------------------------------------------------------------------------
# Reference rotation matrices (right-handed, about the standard axes)
# ---------------------------------------------------------------------------


def Rx(deg: float) -> np.ndarray:
    """Right-handed rotation matrix about +x by deg degrees."""
    t = np.deg2rad(deg)
    c, s = np.cos(t), np.sin(t)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def Ry(deg: float) -> np.ndarray:
    """Right-handed rotation matrix about +y by deg degrees."""
    t = np.deg2rad(deg)
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def Rz(deg: float) -> np.ndarray:
    """Right-handed rotation matrix about +z by deg degrees."""
    t = np.deg2rad(deg)
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


# ---------------------------------------------------------------------------
# Standard basis (You 1999 convention)
# ---------------------------------------------------------------------------

STANDARD_BASIS = {
    "vertical": XHAT,
    "longitudinal": YHAT,
    "transverse": ZHAT,
}
