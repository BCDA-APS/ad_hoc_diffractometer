"""
ad_hoc_diffractometer — multi-circle diffractometer geometry calculations.

Suggested import alias:  import ad_hoc_diffractometer as ahd

Based on:
  Busing & Levy, Acta Cryst. 22, 457-464 (1967)
  H. You, J. Appl. Cryst. 32, 614-623 (1999). DOI: 10.1107/S0021889899001223
  M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
  D.A. Walko, Reference Module in Materials Science and Materials Engineering (2016).
"""

from .axes import axis_from_physical
from .axes import axis_label
from .axes import parse_axis
from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT
from .factories import geometry_fourc
from .factories import geometry_psic
from .factories import geometry_sixc
from .geometry import AdHocDiffractometer
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
    # axes
    "parse_axis",
    "axis_label",
    "axis_from_physical",
    # rotation
    "rotation_matrix",
    # stage
    "Stage",
    # geometry
    "AdHocDiffractometer",
    # factories
    "geometry_psic",
    "geometry_fourc",
    "geometry_sixc",
    # lattice
    "lattice_vectors",
    "reciprocal_vectors",
    "b_matrix",
]
