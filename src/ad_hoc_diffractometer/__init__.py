"""
ad_hoc_diffractometer — Multi-circle diffractometer geometry and related calculations.
"""

import logging

from ._version import __version__
from .axes import axis_from_physical
from .axes import axis_label
from .axes import kappa_axis
from .axes import parse_axis
from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT
from .display import allclose
from .display import fmt
from .display import get_precision
from .display import precision_atol
from .display import set_precision
from .factories import BASIS_BL
from .factories import BASIS_YOU
from .factories import GEOMETRY_ENTRY_POINT_GROUP
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
from .mode import BisectingMode
from .mode import DiffractionMode
from .mode import FixedAngleMode
from .mode import ModeDict
from .orientation import angles_to_phi_vector
from .orientation import ub_from_one_reflection
from .orientation import ub_from_three_reflections_bl1967
from .orientation import ub_from_two_reflections_bl1967
from .orientation import ub_identity
from .refinement import refine_lattice_bl1967
from .refinement import refine_lattice_simplex
from .reflection import Reflection
from .reflection import ReflectionList
from .rotation import rotation_matrix
from .sample import Sample
from .sample import SampleDict
from .spec import FourcG1
from .spec import emit_fourc_g1
from .spec import g1_to_sample
from .spec import parse_fourc_g1
from .spec import sample_to_g1
from .stage import Stage

logger = logging.getLogger(__name__)

logger = logging.getLogger(__name__)

__all__ = [
    # version
    "__version__",
    # constants
    "XHAT",
    "YHAT",
    "ZHAT",
    # display precision
    "get_precision",
    "set_precision",
    "fmt",
    "precision_atol",
    "allclose",
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
    "Sample",
    "SampleDict",
    # modes
    "DiffractionMode",
    "FixedAngleMode",
    "BisectingMode",
    "ModeDict",
    # orientation
    "angles_to_phi_vector",
    "ub_identity",
    "ub_from_one_reflection",
    "ub_from_two_reflections_bl1967",
    "ub_from_three_reflections_bl1967",
    # factories
    "BASIS_BL",
    "BASIS_YOU",
    "GEOMETRY_ENTRY_POINT_GROUP",
    "KAPPA_ALPHA_DEFAULT",
    "list_geometries",
    "get_geometry",
    "make_geometry",
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
    # refinement
    "refine_lattice_bl1967",
    "refine_lattice_simplex",
    # spec
    "FourcG1",
    "parse_fourc_g1",
    "emit_fourc_g1",
    "g1_to_sample",
    "sample_to_g1",
]
