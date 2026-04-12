# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
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
from .engines import Q_mag_to_two_theta
from .engines import Q_to_d
from .engines import Q_to_hkl
from .engines import d_to_Q_mag
from .engines import d_to_two_theta
from .engines import hkl_to_d
from .engines import hkl_to_Q
from .engines import hkl_to_two_theta
from .engines import two_theta_to_d
from .engines import two_theta_to_Q_mag
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
from .forward import compute_forward
from .geometry import AdHocDiffractometer
from .geometry import pa
from .geometry import wh
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
from .radiation import HC_KEV_ANGSTROM
from .radiation import HC_KEV_ANGSTROM_UNCERTAINTY
from .radiation import NEUTRON_MEV_ANGSTROM2
from .radiation import NEUTRON_MEV_ANGSTROM2_UNCERTAINTY
from .radiation import SOURCE_TYPES
from .radiation import XRAY_LINES
from .radiation import energy_to_wavelength
from .radiation import neutron_energy_to_wavelength
from .radiation import neutron_wavelength_to_energy
from .radiation import wavelength_to_energy
from .radiation import wavelength_to_wavenumber
from .radiation import wavenumber_to_wavelength
from .refinement import refine_lattice_bl1967
from .refinement import refine_lattice_simplex
from .reflection import Reflection
from .reflection import ReflectionList
from .rotation import rotation_matrix
from .sample import Sample
from .sample import SampleDict
from .scan import NEAREST_ANGLES
from .scan import hkl_trajectory
from .scan import psi_trajectory
from .scan import trajectory_plan
from .stage import Stage
from .surface import alpha_f
from .surface import alpha_i
from .surface import is_evanescent
from .surface import is_specular
from .surface import q_components

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
    "pa",
    "wh",
    "Reflection",
    "ReflectionList",
    "Sample",
    "SampleDict",
    # modes
    "DiffractionMode",
    "FixedAngleMode",
    "BisectingMode",
    "ModeDict",
    # alternative calculation engines
    "hkl_to_Q",
    "Q_to_hkl",
    "Q_to_d",
    "d_to_Q_mag",
    "hkl_to_d",
    "d_to_two_theta",
    "two_theta_to_d",
    "hkl_to_two_theta",
    "two_theta_to_Q_mag",
    "Q_mag_to_two_theta",
    # forward calculation
    "compute_forward",
    # radiation — constants, X-ray and neutron conversions
    "HC_KEV_ANGSTROM",
    "HC_KEV_ANGSTROM_UNCERTAINTY",
    "NEUTRON_MEV_ANGSTROM2",
    "NEUTRON_MEV_ANGSTROM2_UNCERTAINTY",
    "SOURCE_TYPES",
    "XRAY_LINES",
    "wavelength_to_energy",
    "energy_to_wavelength",
    "wavelength_to_wavenumber",
    "wavenumber_to_wavelength",
    "neutron_wavelength_to_energy",
    "neutron_energy_to_wavelength",
    # scan / trajectory
    "NEAREST_ANGLES",
    "hkl_trajectory",
    "psi_trajectory",
    "trajectory_plan",
    # surface geometry
    "alpha_i",
    "alpha_f",
    "q_components",
    "is_specular",
    "is_evanescent",
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
]
