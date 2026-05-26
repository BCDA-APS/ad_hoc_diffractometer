# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
ad_hoc_diffractometer — Multi-circle diffractometer geometry and related calculations.

Tier 1 (routine) names are imported directly into this namespace.
All other names are accessible via their submodules::

    import ad_hoc_diffractometer as ahd

    # Tier 1 — always available as ahd.<name>
    g = ahd.make_geometry("fourcv")
    lat = ahd.Lattice(a=5.43)
    ahd.ub_identity(g.sample)

    # Tier 2/3 — access via submodule
    ahd.display.set_precision(4)
    ahd.radiation.energy_to_wavelength(12.0)
    ahd.conversions.hkl_to_d(lat.B, 1, 1, 1)
"""

import logging

from ._version import __version__
from .diffractometer import AdHocDiffractometer
from .diffractometer import pa
from .diffractometer import wh
from .factories import get_geometry
from .factories import list_geometries
from .factories import make_geometry
from .factories import register_geometry
from .geometry_loader import load_geometry_file
from .geometry_loader import register_geometry_file
from .geometry_loader import register_geometry_yaml
from .lattice import Lattice
from .mode import REQUIRED
from .mode import BisectConstraint
from .mode import ConstraintSet
from .mode import ConstraintViolation
from .mode import DetectorConstraint
from .mode import EwaldSphereViolation
from .mode import ReferenceConstraint
from .mode import SampleConstraint
from .mode import VirtualBisectConstraint
from .orientation import ub_from_one_reflection
from .orientation import ub_from_three_reflections_bl1967
from .orientation import ub_from_two_reflections_bl1967
from .orientation import ub_identity
from .reflection import Reflection
from .sample import Sample

logger = logging.getLogger(__name__)

__all__ = [
    # version
    "__version__",
    # core classes
    "AdHocDiffractometer",
    "Lattice",
    "Sample",
    "Reflection",
    # status commands
    "pa",
    "wh",
    # geometry registry
    "list_geometries",
    "get_geometry",
    "make_geometry",
    "register_geometry",
    "register_geometry_file",
    "register_geometry_yaml",
    "load_geometry_file",
    # orientation
    "ub_identity",
    "ub_from_one_reflection",
    "ub_from_two_reflections_bl1967",
    "ub_from_three_reflections_bl1967",
    # modes and constraints
    "ConstraintSet",
    "SampleConstraint",
    "DetectorConstraint",
    "BisectConstraint",
    "VirtualBisectConstraint",
    "ReferenceConstraint",
    # exceptions
    "EwaldSphereViolation",
    "ConstraintViolation",
    # mode sentinel
    "REQUIRED",
]
