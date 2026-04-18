# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
presets.py — pre-built diffractometer geometry functions.

This module provides **pre-built geometries**: factory functions that
construct fully configured :class:`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
instances for the most common multi-circle diffractometer designs used in
synchrotron and laboratory X-ray / neutron crystallography.  They also
serve as worked examples for defining custom geometries.

Pre-built geometries
--------------------

.. list-table::
   :header-rows: 1
   :widths: 15 20 50

   * - Name
     - Type
     - Description
   * - :func:`fourcv`
     - Eulerian 4-circle
     - Vertical scattering plane (synchrotron). See :doc:`/geometries/fourcv`.
   * - :func:`fourch`
     - Eulerian 4-circle
     - Horizontal scattering plane (laboratory). See :doc:`/geometries/fourch`.
   * - :func:`fivec`
     - Eulerian 5-circle
     - Vlieg et al. (1987), fourcv on a mu base. See :doc:`/geometries/fivec`.
   * - :func:`psic`
     - Eulerian 6-circle
     - You (1999) 4S+2D. See :doc:`/geometries/psic`.
   * - :func:`sixc`
     - Eulerian 6-circle
     - Lohmeier & Vlieg (1993) surface geometry. See :doc:`/geometries/sixc`.
   * - :func:`kappa4cv`
     - Kappa 4-circle
     - Vertical scattering plane (synchrotron). See :doc:`/geometries/kappa4cv`.
   * - :func:`kappa4ch`
     - Kappa 4-circle
     - Horizontal scattering plane (laboratory). See :doc:`/geometries/kappa4ch`.
   * - :func:`kappa6c`
     - Kappa 6-circle
     - Psic-style outer axes. See :doc:`/geometries/kappa6c`.
   * - :func:`zaxis`
     - Surface / special
     - Bloch (1985) Z-axis geometry. See :doc:`/geometries/zaxis`.
   * - :func:`s2d2`
     - Surface / special
     - Evans-Lutterodt & Tang (1995). See :doc:`/geometries/s2d2`.

Each factory is decorated with
:func:`~ad_hoc_diffractometer.factories.register_geometry`, which registers
it in the geometry registry.  Use
:func:`~ad_hoc_diffractometer.factories.list_geometries` to retrieve all
registered geometries as a ``{name: callable}`` dict.

Usage
-----
::

    import ad_hoc_diffractometer as ahd

    g = ahd.presets.fourcv()
    g.wavelength = 1.54

Naming convention:

- Eulerian geometries: ``fourcv``, ``fourch``, ``fivec``, ``psic``, ``sixc``
- Kappa geometries: ``kappa4cv``, ``kappa4ch``, ``kappa6c``
- Surface / special geometries: ``zaxis``, ``s2d2``

Scattering-plane suffix convention (``v`` / ``h``):

- ``v`` — vertical scattering plane (synchrotron): ttheta rotates about the lateral axis
- ``h`` — horizontal scattering plane (laboratory): ttheta rotates about the vertical axis

Where no suffix is given, the detector convention is unambiguous (psic, sixc)
or the geometry uses a non-standard detector arrangement (zaxis, s2d2, fivec).

    fourch  (horizontal scattering plane) matches Busing & Levy (1967) Fig. 1b.
    fourcv  (vertical   scattering plane) is the synchrotron convention.
    kappa4ch and kappa4cv follow the same convention as fourch and fourcv.

Kappa angle (alpha) convention (Walko 2016; Enraf-Nonius; ITC Vol. C Sec. 2.2.6):
    The kappa axis lies in the vertical-lateral plane, tilted alpha degrees
    from the vertical axis toward the lateral axis.  Typical value: 50 deg.

Walko (2016) designations:
    S3D1      fourcv, fourch, kappa4cv, kappa4ch
    S4D2      psic, kappa6c
    (S3D2)1   sixc   (sample and detector share base stage)
    (S3D1)1   fivec  (fourcv mounted on vertical base)
    (S1D2)1   zaxis  (sample and detector share alpha base stage)
    S2D2      s2d2   (fully decoupled sample/detector pairs)

References
----------
* W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967) — fourcv / fourch
* J.M. Bloch, J. Appl. Cryst. 18, 33-36 (1985) — zaxis
* E. Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987) — fivec
* M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993) — sixc
* K.W. Evans-Lutterodt & M.-T. Tang, J. Appl. Cryst. 28, 318-326 (1995) — s2d2
* H. You, J. Appl. Cryst. 32, 614-623 (1999)
  DOI:10.1107/S0021889899001223 — psic
* ITC Vol. C, Sec. 2.2.6 (2006)
  DOI:10.1107/97809553602060000577 — kappa
* D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016) — kappa, zaxis, s2d2, fivec
"""

from __future__ import annotations

import inspect

from .axes import kappa_axis
from .factories import BASIS_BL
from .factories import BASIS_YOU
from .factories import KAPPA_ALPHA_DEFAULT
from .factories import register_geometry
from .geometry import AdHocDiffractometer
from .mode import REQUIRED
from .mode import BisectConstraint
from .mode import ConstraintSet
from .mode import DetectorConstraint
from .mode import ReferenceConstraint
from .mode import SampleConstraint
from .stage import Stage

__all__ = [
    "fivec",
    "fourch",
    "fourcv",
    "kappa4ch",
    "kappa4cv",
    "kappa6c",
    "psic",
    "s2d2",
    "sixc",
    "zaxis",
]

# ---------------------------------------------------------------------------
# Eulerian geometries
# ---------------------------------------------------------------------------


@register_geometry
def psic(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """
    You (1999) '4S+2D' six-circle diffractometer (psic geometry).

    Walko (2016) designation: S4D2.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, lateral=+z.
    Detector axis: lateral.  Vertical scattering plane.

    Sample stack (floor first):
        mu  : vertical,     right-handed
        eta : lateral,      left-handed
        chi : longitudinal, right-handed
        phi : lateral,      left-handed

    Detector stack (floor first):
        nu    : vertical, right-handed
        delta : lateral,  left-handed

    mu and nu share the same vertical rotation axis; mechanically independent.

    Reference: H. You, J. Appl. Cryst. 32, 614-623 (1999).
               DOI: 10.1107/S0021889899001223
    """
    VERTICAL = basis["vertical"]
    LATERAL = basis["lateral"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("eta", -LATERAL, parent="mu", role="sample"),
        Stage("chi", +LONGITUDINAL, parent="eta", role="sample"),
        Stage("phi", -LATERAL, parent="chi", role="sample"),
        Stage("nu", +VERTICAL, parent=None, role="detector"),
        Stage("delta", -LATERAL, parent="nu", role="detector"),
    ]
    # psic: 6 DOF, N-3=3 constraints needed per mode (You 1999, S4D2).
    # Vertical bisect pair: eta(lateral) <-> delta(lateral)  => eta = delta/2
    # Horizontal bisect pair: mu(vertical) <-> nu(vertical)  => mu = nu/2
    modes = {
        # ── Implemented analytic modes ──────────────────────────────────────
        "bisecting_vertical": ConstraintSet(
            [
                BisectConstraint("eta", "delta"),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
        ),
        "fixed_chi": ConstraintSet(
            [
                SampleConstraint("chi", 90.0),
                BisectConstraint("eta", "delta"),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "phi", "delta"],
        ),
        "fixed_phi": ConstraintSet(
            [
                SampleConstraint("phi", 0.0),
                BisectConstraint("eta", "delta"),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "delta"],
        ),
        "fixed_mu": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                BisectConstraint("eta", "delta"),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
        ),
        "bisecting_horizontal": ConstraintSet(
            [
                BisectConstraint("mu", "nu"),
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
        ),
        "fixed_nu": ConstraintSet(
            [
                DetectorConstraint("nu", 0.0),
                BisectConstraint("eta", "delta"),
                SampleConstraint("mu", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
        ),
        "double_diffraction_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
        "double_diffraction_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
        "lifting_detector_mu": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("eta", 0.0),
                DetectorConstraint("qaz", 90.0),
            ],
            computed=["mu", "nu", "delta"],
        ),
        "lifting_detector_phi": ConstraintSet(
            [
                SampleConstraint("phi", 0.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("qaz", 90.0),
            ],
            computed=["phi", "nu", "delta"],
        ),
        "psi_constant_vertical": ConstraintSet(
            [
                BisectConstraint("eta", "delta"),
                SampleConstraint("mu", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "psi_constant_horizontal": ConstraintSet(
            [
                BisectConstraint("mu", "nu"),
                SampleConstraint("eta", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=BASIS_YOU,
        description=(
            "You (1999) 4S+2D six-circle diffractometer "
            "(lateral detector, vertical scattering plane, synchrotron)"
        ),
        modes=modes,
        default_mode="bisecting_vertical",
    )


@register_geometry
def fourcv(basis: dict = BASIS_BL) -> AdHocDiffractometer:
    """
    Four-circle Eulerian diffractometer, vertical scattering plane (synchrotron).

    Walko (2016) designation: S3D1.

    Synchrotron configuration: omega and ttheta both rotate about the
    lateral axis, so the scattering plane is vertical.  This exploits the
    s-polarisation and tighter vertical collimation of synchrotron radiation
    (Walko 2016).

    Default basis: Busing & Levy (1967) — lateral=+x, longitudinal=+y, vertical=+z.

    Sample stack (floor first):
        omega  : lateral,      left-handed
        chi    : longitudinal, right-handed
        phi    : lateral,      left-handed

    Detector (floor, mechanically independent of sample stack):
        ttheta : lateral,      left-handed

    omega and ttheta share the same lateral axis; mechanically independent.

    References: W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
                D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
    """
    LATERAL = basis["lateral"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("omega", -LATERAL, parent=None, role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -LATERAL, parent="chi", role="sample"),
        Stage("ttheta", -LATERAL, parent=None, role="detector"),
    ]
    # fourcv: 4 DOF, N-3=1 constraint needed per mode.
    modes = {
        "bisecting": ConstraintSet(
            [BisectConstraint("omega", "ttheta")],
            computed=["omega", "chi", "phi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["omega", "phi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["omega", "chi", "ttheta"],
        ),
        "constant_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["chi", "phi", "ttheta"],
        ),
        "psi_constant": ConstraintSet(
            [ReferenceConstraint("psi", 0.0)],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction": ConstraintSet(
            [],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=BASIS_BL,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(vertical scattering plane, lateral ttheta, synchrotron)"
        ),
        modes=modes,
        default_mode="bisecting",
    )


@register_geometry
def fourch(basis: dict = BASIS_BL) -> AdHocDiffractometer:
    """
    Four-circle Eulerian diffractometer, horizontal scattering plane (laboratory).

    Walko (2016) designation: S3D1.

    Laboratory / default configuration: omega and ttheta both rotate about
    the vertical axis, so the scattering plane is horizontal.  This is the
    geometry described in Busing & Levy (1967), Fig. 1b.

    Default basis: Busing & Levy (1967) — lateral=+x, longitudinal=+y, vertical=+z.

    Sample stack (floor first):
        omega  : vertical, left-handed
        chi    : lateral,  right-handed
        phi    : vertical, left-handed

    Detector (floor, mechanically independent):
        ttheta : vertical, left-handed

    omega and ttheta share the same vertical axis; mechanically independent.

    Reference: W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
               D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
    """
    VERTICAL = basis["vertical"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("omega", -VERTICAL, parent=None, role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -VERTICAL, parent="chi", role="sample"),
        Stage("ttheta", -VERTICAL, parent=None, role="detector"),
    ]
    # fourch: 4 DOF, N-3=1 constraint needed per mode.
    modes = {
        "bisecting": ConstraintSet(
            [BisectConstraint("omega", "ttheta")],
            computed=["omega", "chi", "phi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["omega", "phi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["omega", "chi", "ttheta"],
        ),
        "constant_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["chi", "phi", "ttheta"],
        ),
        "psi_constant": ConstraintSet(
            [ReferenceConstraint("psi", 0.0)],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction": ConstraintSet(
            [],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(horizontal scattering plane, vertical ttheta, laboratory)"
        ),
        modes=modes,
        default_mode="bisecting",
    )


@register_geometry
def sixc(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """
    Lohmeier & Vlieg (1993) six-circle surface diffractometer (sixc geometry).

    Also known as the IUCr six-circle diffractometer.
    Walko (2016) designation: (S3D2)1.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, lateral=+z.

    Sample and detector stacks share the alpha (rotary table) base stage,
    making this a coupled geometry.  Useful for surface diffraction.

    From Fig. 1 and §2.1 of Lohmeier & Vlieg (1993):
    alpha and gamma rotate about the vertical axis (x in LV convention).
    omega, phi, and delta all rotate about the lateral axis (z in LV).
    chi rotates about the longitudinal axis (y in LV).

    Stack (floor first)::

        alpha (shared base): vertical,     right-handed  [rotary table]
          --> omega (sample):  lateral,      left-handed
                --> chi:       longitudinal, right-handed
                      --> phi: lateral,      left-handed
          --> delta (detector): lateral,     left-handed
                --> gamma:      vertical,    right-handed

    Reference: M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
    """
    VERTICAL = basis["vertical"]
    LATERAL = basis["lateral"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("alpha", +VERTICAL, parent=None, role="sample"),
        Stage("omega", -LATERAL, parent="alpha", role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -LATERAL, parent="chi", role="sample"),
        Stage("delta", -LATERAL, parent="alpha", role="detector"),
        Stage("gamma", +VERTICAL, parent="delta", role="detector"),
    ]
    # sixc: 6 DOF, N-3=3 constraints needed per mode.
    # four_circle: alpha=0, gamma=0 frozen; bisect on omega/delta.
    # Reduces to fourcv bisecting — the generic bisecting solver handles it.
    # Surface modes (five_circle_*, zaxis_*) require reference infrastructure (Issue J).
    modes = {
        "bisecting_4c": ConstraintSet(
            [
                SampleConstraint("alpha", 0.0),
                DetectorConstraint("gamma", 0.0),
                BisectConstraint("omega", "delta"),
            ],
            computed=["omega", "chi", "phi", "delta"],
        ),
        "fixed_gamma_5c": ConstraintSet(
            [
                DetectorConstraint("gamma", 0.0),
                SampleConstraint("alpha", 0.0),
                BisectConstraint("omega", "delta"),
            ],
            computed=["omega", "chi", "phi", "delta", "alpha"],
        ),
        "fixed_alpha_5c": ConstraintSet(
            [
                SampleConstraint("alpha", 0.0),
                BisectConstraint("omega", "delta"),
                DetectorConstraint("gamma", 0.0),
            ],
            computed=["omega", "chi", "phi", "delta", "gamma"],
        ),
        "fixed_alpha_zaxis": ConstraintSet(
            [
                SampleConstraint("alpha", 0.0),
                SampleConstraint("chi", 0.0),
                ReferenceConstraint("alpha_i", 0.0),
            ],
            computed=["omega", "delta", "gamma"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_beta_zaxis": ConstraintSet(
            [
                DetectorConstraint("gamma", 0.0),
                SampleConstraint("chi", 0.0),
                ReferenceConstraint("beta_out", 0.0),
            ],
            computed=["omega", "delta", "alpha"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "alpha_eq_beta_zaxis": ConstraintSet(
            [
                SampleConstraint("chi", 0.0),
                SampleConstraint("phi", 0.0),
                ReferenceConstraint("a_eq_b", True),
            ],
            computed=["omega", "delta", "alpha", "gamma"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            "Lohmeier & Vlieg (1993) six-circle surface diffractometer "
            "(IUCr six-circle; Walko S(3D2)1). "
            "Sample and detector share the alpha (rotary table) base stage."
        ),
        modes=modes,
        default_mode="bisecting_4c",
    )


# ---------------------------------------------------------------------------
# Kappa geometries
# ---------------------------------------------------------------------------


@register_geometry
def kappa4cv(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
    basis: dict = BASIS_BL,
) -> AdHocDiffractometer:
    """
    Four-circle kappa diffractometer, vertical scattering plane (synchrotron).

    Walko (2016) designation: S3D1.

    The chi circle of a standard Eulerian fourcv is replaced by a kappa arm.
    The kappa axis lies in the vertical-lateral plane, tilted alpha degrees
    from the vertical toward the lateral axis.

    komega and ttheta both rotate about the lateral axis; the scattering
    plane is vertical (synchrotron convention).

    Default basis: Busing & Levy (1967) — lateral=+x, longitudinal=+y, vertical=+z.

    Sample stack (floor first):
        komega : lateral, left-handed
        kappa  : tilted,  kappa_axis(alpha), right-handed
        kphi   : lateral, left-handed

    Detector (floor, mechanically independent):
        ttheta : lateral, left-handed

    komega and ttheta share the same lateral axis; mechanically independent.

    Parameters
    ----------
    alpha_deg : float
        Kappa tilt angle in degrees (default 50).  Must be in (0, 90).

    References: D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
                W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
                ITC Vol. C, Sec. 2.2.6 (2006).
    """
    LATERAL = basis["lateral"]
    kax = kappa_axis(alpha_deg, basis=basis)
    stages = [
        Stage("komega", -LATERAL, parent=None, role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -LATERAL, parent="kappa", role="sample"),
        Stage("ttheta", -LATERAL, parent=None, role="detector"),
    ]
    # kappa4cv: 4 DOF, N-3=1 constraint needed per mode.
    # Virtual Eulerian angles (omega, chi, phi) are computed from real kappa
    # angles (komega, kappa, kphi) via Walko (2016) eq. [16].  Constraints
    # using virtual angle names are stubs pending Issue I (#153).
    # BisectConstraint('komega','ttheta') approximates bisecting but is
    # physically inaccurate — true bisect requires virtual omega=0.
    modes = {
        "bisecting": ConstraintSet(
            [BisectConstraint("komega", "ttheta")],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_kphi": ConstraintSet(
            [SampleConstraint("kphi", 0.0)],
            computed=["komega", "kappa", "ttheta"],
        ),
        "constant_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "constant_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "constant_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "psi_constant": ConstraintSet(
            [ReferenceConstraint("psi", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction": ConstraintSet(
            [],
            computed=["komega", "kappa", "kphi", "ttheta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            f"Four-circle kappa diffractometer, vertical scattering plane "
            f"(synchrotron). Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
        modes=modes,
        default_mode="bisecting",
    )


@register_geometry
def kappa4ch(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
    basis: dict = BASIS_BL,
) -> AdHocDiffractometer:
    """
    Four-circle kappa diffractometer, horizontal scattering plane (laboratory).

    Walko (2016) designation: S3D1.

    Identical to kappa4cv() but komega and ttheta rotate about the vertical
    axis, giving a horizontal scattering plane (laboratory convention).

    Default basis: Busing & Levy (1967) — lateral=+x, longitudinal=+y, vertical=+z.

    Sample stack (floor first):
        komega : vertical, left-handed
        kappa  : tilted,   kappa_axis(alpha), right-handed
        kphi   : vertical, left-handed

    Detector (floor, mechanically independent):
        ttheta : vertical, left-handed

    komega and ttheta share the same vertical axis; mechanically independent.

    Parameters
    ----------
    alpha_deg : float
        Kappa tilt angle in degrees (default 50).  Must be in (0, 90).

    References: D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
                W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
                ITC Vol. C, Sec. 2.2.6 (2006).
    """
    VERTICAL = basis["vertical"]
    kax = kappa_axis(alpha_deg, basis=basis)
    stages = [
        Stage("komega", -VERTICAL, parent=None, role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -VERTICAL, parent="kappa", role="sample"),
        Stage("ttheta", -VERTICAL, parent=None, role="detector"),
    ]
    # kappa4ch: 4 DOF, N-3=1 constraint needed per mode.
    # Same mode set as kappa4cv; virtual angle stubs pending Issue I (#153).
    modes = {
        "bisecting": ConstraintSet(
            [BisectConstraint("komega", "ttheta")],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_kphi": ConstraintSet(
            [SampleConstraint("kphi", 0.0)],
            computed=["komega", "kappa", "ttheta"],
        ),
        "constant_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "constant_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "constant_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "psi_constant": ConstraintSet(
            [ReferenceConstraint("psi", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            f"Four-circle kappa diffractometer, horizontal scattering plane "
            f"(laboratory). Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
        modes=modes,
        default_mode="bisecting",
    )


@register_geometry
def kappa6c(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
    basis: dict = BASIS_YOU,
) -> AdHocDiffractometer:
    """
    Six-circle kappa diffractometer (psic-style outer axes, kappa inner sample).

    Walko (2016) designation: S4D2.

    Extends the kappa4cv geometry with two additional axes (mu, nu) in
    the style of the psic geometry (You 1999), giving full orientation freedom.
    This is the synchrotron configuration with a lateral detector.

    Default basis: You (1999) — vertical=+x, longitudinal=+y, lateral=+z.

    Sample stack (floor first):
        mu     : vertical,     right-handed   [outermost]
        komega : lateral,      left-handed
        kappa  : tilted,       kappa_axis(alpha), right-handed
        kphi   : lateral,      left-handed

    Detector stack (floor first):
        nu     : vertical,     right-handed
        delta  : lateral,      left-handed

    mu and nu share the same vertical axis; mechanically independent.

    Parameters
    ----------
    alpha_deg : float
        Kappa tilt angle in degrees (default 50).  Must be in (0, 90).

    References: D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
                H. You, J. Appl. Cryst. 32, 614-623 (1999).
                ITC Vol. C, Sec. 2.2.6 (2006).
    """
    VERTICAL = basis["vertical"]
    LATERAL = basis["lateral"]
    kax = kappa_axis(alpha_deg, basis=basis)
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("komega", -LATERAL, parent="mu", role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -LATERAL, parent="kappa", role="sample"),
        Stage("nu", +VERTICAL, parent=None, role="detector"),
        Stage("delta", -LATERAL, parent="nu", role="detector"),
    ]
    # kappa6c: 6 DOF, N-3=3 constraints needed per mode.
    # Vertical bisect pair: komega(lateral) <-> delta(lateral) => komega = delta/2
    #   (approximates virtual omega_euler = delta/2; corrected by Issue I / #153)
    # Horizontal bisect pair: mu(vertical) <-> nu(vertical) => mu = nu/2
    # Virtual Eulerian angles (omega, chi, phi) via Walko (2016) eq. [16].
    modes = {
        # ── Implemented (generic solver) ────────────────────────────────────
        "bisecting_vertical": ConstraintSet(
            [
                BisectConstraint("komega", "delta"),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
        ),
        "bisecting_horizontal": ConstraintSet(
            [
                BisectConstraint("mu", "nu"),
                SampleConstraint("komega", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "kappa", "kphi", "nu"],
        ),
        "fixed_kphi": ConstraintSet(
            [
                SampleConstraint("kphi", 0.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["komega", "kappa", "delta"],
        ),
        "fixed_mu": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                BisectConstraint("komega", "delta"),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
        ),
        "fixed_nu": ConstraintSet(
            [
                DetectorConstraint("nu", 0.0),
                BisectConstraint("komega", "delta"),
                SampleConstraint("mu", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
        ),
        "fixed_delta": ConstraintSet(
            [
                DetectorConstraint("delta", 0.0),
                BisectConstraint("mu", "nu"),
                SampleConstraint("komega", 0.0),
            ],
            computed=["mu", "kappa", "kphi", "nu"],
        ),
        # ── Stubs: solver not yet implemented ───────────────────────────────
        "lifting_detector_mu": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("komega", 0.0),
                DetectorConstraint("qaz", 90.0),
            ],
            computed=["mu", "nu", "delta"],
        ),
        "lifting_detector_kphi": ConstraintSet(
            [
                SampleConstraint("kphi", 0.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("qaz", 90.0),
            ],
            computed=["kphi", "nu", "delta"],
        ),
        "psi_constant_vertical": ConstraintSet(
            [
                BisectConstraint("komega", "delta"),
                SampleConstraint("mu", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "psi_constant_horizontal": ConstraintSet(
            [
                BisectConstraint("mu", "nu"),
                SampleConstraint("komega", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["mu", "kappa", "kphi", "nu"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
        "double_diffraction_horizontal": ConstraintSet(
            [
                SampleConstraint("komega", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "kappa", "kphi", "nu"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            f"Six-circle kappa diffractometer, psic-style outer axes "
            f"(lateral detector, synchrotron). "
            f"Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
        modes=modes,
        default_mode="bisecting_vertical",
    )


# ---------------------------------------------------------------------------
# General-inclination geometries (Walko 2016, Sections 4.2 and 5)
# ---------------------------------------------------------------------------


@register_geometry
def zaxis(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """
    Z-axis four-circle diffractometer (general-inclination geometry).

    Walko (2016) designation: (S1D2)1.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, lateral=+z.

    Designed for surface diffraction.  The sample surface normal is aligned
    parallel to the Z-axis, so the angle of incidence equals the alpha angle.
    The detector and sample both rotate about the shared alpha (base) axis.

    Stack (floor first)::

        alpha (shared base): vertical,     right-handed
          --> Z     (sample)  : longitudinal, right-handed
          --> delta (detector): lateral,      left-handed
                --> gamma :     vertical,     right-handed

    The total scattering angle is a compound of gamma, delta, and alpha
    (Walko 2016, eq. 17)::

        ttheta = arccos(cos(gamma)*cos(delta)*cos(alpha) + sin(alpha)*sin(gamma))

    References:
    J.M. Bloch, J. Appl. Cryst. 18, 33-36 (1985).
    D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016), eq. 17.
    """
    VERTICAL = basis["vertical"]
    LATERAL = basis["lateral"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("alpha", +VERTICAL, parent=None, role="sample"),
        Stage("Z", +LONGITUDINAL, parent="alpha", role="sample"),
        Stage("delta", -LATERAL, parent="alpha", role="detector"),
        Stage("gamma", +VERTICAL, parent="delta", role="detector"),
    ]
    # zaxis: 4 DOF, N-3=1 constraint needed per mode.
    # All modes require reference vector n̂ (Issue J / #157).
    modes = {
        "zaxis": ConstraintSet(
            [ReferenceConstraint("alpha_i", 0.0)],
            computed=["Z", "delta", "gamma"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "reflectivity": ConstraintSet(
            [ReferenceConstraint("a_eq_b", True)],
            computed=["Z", "delta", "alpha", "gamma"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            "Z-axis four-circle diffractometer (Bloch 1985; Walko 2016 (S1D2)1). "
            "Surface normal parallel to Z-axis. "
            "Sample and detector share the alpha (base) stage."
        ),
        modes=modes,
    )


@register_geometry
def s2d2(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """
    S2D2 four-circle diffractometer (general-inclination geometry).

    Walko (2016) designation: S2D2.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, lateral=+z.

    Two independent sample axes (mu, Z) and two independent detector axes
    (nu, delta), all mechanically decoupled.  The angle of incidence is the
    mu angle; the surface normal is parallel to Z.

    Sample stack (floor first)::

        mu    : vertical,     right-handed
          --> Z : longitudinal, right-handed

    Detector stack (floor first)::

        nu    : vertical,     right-handed
          --> delta : lateral, left-handed

    mu and nu share the same vertical axis; mechanically independent.

    The total scattering angle is (Walko 2016, eq. 18)::

        ttheta = arccos(cos(nu) * cos(delta))

    References:
    K.W. Evans-Lutterodt & M.-T. Tang, J. Appl. Cryst. 28, 318-326 (1995).
    D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016), eq. 18.
    """
    VERTICAL = basis["vertical"]
    LATERAL = basis["lateral"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("Z", +LONGITUDINAL, parent="mu", role="sample"),
        Stage("nu", +VERTICAL, parent=None, role="detector"),
        Stage("delta", -LATERAL, parent="nu", role="detector"),
    ]
    # s2d2: 4 DOF, N-3=1 constraint needed per mode.
    modes = {
        "mu_fixed": ConstraintSet(
            [SampleConstraint("mu", 0.0)],
            computed=["Z", "nu", "delta"],
        ),
        "reflectivity": ConstraintSet(
            [ReferenceConstraint("a_eq_b", True)],
            computed=["mu", "Z", "nu", "delta"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            "S2D2 four-circle diffractometer (Evans-Lutterodt & Tang 1995; "
            "Walko 2016 S2D2). "
            "Fully decoupled sample (mu, Z) and detector (nu, delta) axes."
        ),
        modes=modes,
    )


@register_geometry
def fivec(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """
    Five-circle diffractometer (fourcv mounted on a vertical base).

    Walko (2016) designation: (S3D1)1.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, lateral=+z.

    A standard Eulerian four-circle (fourcv) is mounted on a fifth vertical
    rotation stage (mu) as a base.  The sample and detector motions are coupled
    through mu.  This provides an additional degree of freedom for accessing
    wider regions of reciprocal space, particularly at synchrotron sources.

    Stack (floor first)::

        mu (shared base): vertical,     right-handed
          --> omega (sample): lateral,      left-handed
                --> chi:      longitudinal, right-handed
                      --> phi: lateral,     left-handed
          --> ttheta (detector): lateral,   left-handed

    References:
    E. Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987).
    D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
    """
    VERTICAL = basis["vertical"]
    LATERAL = basis["lateral"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("omega", -LATERAL, parent="mu", role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -LATERAL, parent="chi", role="sample"),
        Stage("ttheta", -LATERAL, parent="mu", role="detector"),
    ]
    # fivec: 5 DOF, N-3=2 constraints needed per mode.
    # With mu=0 the geometry reduces to fourcv; the bisecting solver
    # handles this case identically to fourcv bisecting.
    # Modes where mu != 0 require a tilted-plane solver (not yet implemented).
    modes = {
        "bisecting_4c": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                BisectConstraint("omega", "ttheta"),
            ],
            computed=["omega", "chi", "phi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("chi", 90.0),
            ],
            computed=["omega", "phi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("phi", 0.0),
            ],
            computed=["omega", "chi", "ttheta"],
        ),
        "fixed_mu": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                BisectConstraint("omega", "ttheta"),
            ],
            computed=["omega", "chi", "phi", "ttheta"],
        ),
        "constant_omega_noncoplanar": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("omega", 0.0),
            ],
            computed=["mu", "chi", "phi", "ttheta"],
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            "Five-circle diffractometer: fourcv on vertical mu base "
            "(Vlieg et al. 1987; Walko 2016 (S3D1)1). "
            "Sample and detector coupled through mu."
        ),
        modes=modes,
        default_mode="bisecting_4c",
    )
