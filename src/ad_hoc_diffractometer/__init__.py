"""
ad_hoc_diffractometer — multi-circle diffractometer geometry calculations.

Suggested import alias:  import ad_hoc_diffractometer as ahd

Based on:
  Busing & Levy, Acta Cryst. 22, 457-464 (1967)
  H. You, J. Appl. Cryst. 32, 614-623 (1999). DOI: 10.1107/S0021889899001223
  M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
  D.A. Walko, Reference Module in Materials Science and Materials Engineering (2016).
"""

from .ad_hoc_diffractometer import (
    XHAT,
    YHAT,
    ZHAT,
    rotation_matrix,
    Stage,
    AdHocDiffractometer,
    geometry_psic,
    geometry_fourc,
    geometry_sixc,
    lattice_vectors,
    reciprocal_vectors,
    b_matrix,
)

__all__ = [
    "XHAT",
    "YHAT",
    "ZHAT",
    "rotation_matrix",
    "Stage",
    "AdHocDiffractometer",
    "geometry_psic",
    "geometry_fourc",
    "geometry_sixc",
    "lattice_vectors",
    "reciprocal_vectors",
    "b_matrix",
]
