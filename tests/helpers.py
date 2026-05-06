"""
helpers.py — shared test utilities for the test suite.

Pure functions and constants used in parametrize() decorators, where
pytest fixtures cannot be used (collection-time evaluation).

Import explicitly in any test file that needs them:
    from helpers import Rx, Ry, Rz, STANDARD_BASIS
    from helpers import fourcv, fourch, psic, sixc      # demo geometries
    from helpers import kappa4cv, kappa4ch, kappa6c     # kappa demos
    from helpers import zaxis, s2d2, fivec              # special demos
"""

import numpy as np

from ad_hoc_diffractometer import make_geometry as _make_geometry
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


# ---------------------------------------------------------------------------
# Demo-geometry helper wrappers.
#
# Issue #267 deleted ``ad_hoc_diffractometer.presets`` after migrating
# every demo geometry to a declarative YAML file.  These wrappers
# provide a one-line import path that the test suite can use in place
# of the old ``from helpers import <name>``
# statement; each is a thin call to ``ahd.make_geometry``.
# ---------------------------------------------------------------------------


def fourcv(**kwargs):
    """Return the declarative ``fourcv`` demo geometry."""
    return _make_geometry("fourcv", **kwargs)


def fourch(**kwargs):
    """Return the declarative ``fourch`` demo geometry."""
    return _make_geometry("fourch", **kwargs)


def fivec(**kwargs):
    """Return the declarative ``fivec`` demo geometry."""
    return _make_geometry("fivec", **kwargs)


def psic(**kwargs):
    """Return the declarative ``psic`` demo geometry."""
    return _make_geometry("psic", **kwargs)


def sixc(**kwargs):
    """Return the declarative ``sixc`` demo geometry."""
    return _make_geometry("sixc", **kwargs)


def kappa4cv(**kwargs):
    """Return the declarative ``kappa4cv`` demo geometry."""
    return _make_geometry("kappa4cv", **kwargs)


def kappa4ch(**kwargs):
    """Return the declarative ``kappa4ch`` demo geometry."""
    return _make_geometry("kappa4ch", **kwargs)


def kappa6c(**kwargs):
    """Return the declarative ``kappa6c`` demo geometry."""
    return _make_geometry("kappa6c", **kwargs)


def zaxis(**kwargs):
    """Return the declarative ``zaxis`` demo geometry."""
    return _make_geometry("zaxis", **kwargs)


def s2d2(**kwargs):
    """Return the declarative ``s2d2`` demo geometry."""
    return _make_geometry("s2d2", **kwargs)
