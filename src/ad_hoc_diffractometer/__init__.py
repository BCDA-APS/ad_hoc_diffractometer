"""
ad_hoc_diffractometer — Multi-circle diffractometer geometry and related calculations.
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
from .factories import fourch
from .factories import fourcv
from .factories import get_geometry
from .factories import kappa4ch
from .factories import kappa4cv
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
from .reflection import Reflection
from .reflection import ReflectionList
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
    "Reflection",
    "ReflectionList",
    # factories
    "list_geometries",
    "get_geometry",
    "make_geometry",
    "KAPPA_ALPHA_DEFAULT",
    "psic",
    "fourcv",
    "fourch",
    "sixc",
    "kappa4cv",
    "kappa4ch",
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
