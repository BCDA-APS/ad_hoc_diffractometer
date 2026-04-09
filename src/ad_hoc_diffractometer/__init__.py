"""
ad_hoc_diffractometer — multi-circle diffractometer geometry calculations.

Suggested import alias:  import ad_hoc_diffractometer as ahd

Based on (chronological):
  Busing & Levy, Acta Cryst. 22, 457-464 (1967)
  J.M. Bloch, J. Appl. Cryst. 18, 33-36 (1985).
  E. Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987).
  M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
  K.W. Evans-Lutterodt & M.-T. Tang, J. Appl. Cryst. 28, 318-326 (1995).
  H. You, J. Appl. Cryst. 32, 614-623 (1999). DOI: 10.1107/S0021889899001223
  International Tables for Crystallography, Vol. C, Section 2.2.6 (2006).
    DOI: 10.1107/97809553602060000577
  D.A. Walko, Reference Module in Materials Science and Materials Engineering (2016).
"""

from .axes import axis_from_physical
from .axes import axis_label
from .axes import kappa_axis
from .axes import parse_axis
from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT
from .display import fmt
from .display import get_precision
from .display import set_precision
from .factories import KAPPA_ALPHA_DEFAULT
from .factories import fivec
from .factories import fourc_h
from .factories import fourc_v
from .factories import get_geometry
from .factories import kappa4c
from .factories import kappa4c_h
from .factories import kappa6c
from .factories import list_geometries
from .factories import make_geometry
from .factories import psic
from .factories import s2d2
from .factories import sixc
from .factories import zaxis
from .geometry import AdHocDiffractometer
from .lattice import CRYSTAL_SYSTEMS
from .lattice import Lattice
from .lattice import b_matrix
from .lattice import lattice_vectors
from .lattice import reciprocal_vectors
from .rotation import rotation_matrix
from .stage import Stage

__all__ = [
    # constants
    "XHAT",
    "YHAT",
    "ZHAT",
    # display precision
    "get_precision",
    "set_precision",
    "fmt",
    # axes
    "parse_axis",
    "axis_label",
    "axis_from_physical",
    "kappa_axis",
    # rotation
    "rotation_matrix",
    # stage
    "Stage",
    # geometry
    "AdHocDiffractometer",
    # factories
    "list_geometries",
    "get_geometry",
    "make_geometry",
    "KAPPA_ALPHA_DEFAULT",
    "psic",
    "fourc_v",
    "fourc_h",
    "sixc",
    "kappa4c",
    "kappa4c_h",
    "kappa6c",
    "zaxis",
    "s2d2",
    "fivec",
    # lattice
    "CRYSTAL_SYSTEMS",
    "Lattice",
    "lattice_vectors",
    "reciprocal_vectors",
    "b_matrix",
]
