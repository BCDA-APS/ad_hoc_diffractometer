# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
presets.py — Demonstrations of diffractometer geometries.

This module demonstrates some of the most common diffractometer
geometries used in synchrotron and laboratory X-ray / neutron
crystallography, each fully-configured, using
:class:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`.
Together, the demos here serve as worked examples for defining custom
geometries.

Demo geometries
---------------

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

Each demo geometry is decorated with
:func:`~ad_hoc_diffractometer.factories.register_geometry`, which registers
it in the geometry registry.  Use
:func:`~ad_hoc_diffractometer.factories.list_geometries` to retrieve all
registered geometries as a ``{name: callable}`` dict.

For historical context on diffractometer design and usage that
predates the closed-form Busing & Levy (1967) angle-setting equations,
see Arndt & Willis (1966) — listed in the project's :doc:`/references`.

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

- ``v`` — vertical scattering plane (synchrotron): ttheta rotates about the transverse axis
- ``h`` — horizontal scattering plane (laboratory): ttheta rotates about the vertical axis

Where no suffix is given, the detector convention is unambiguous (psic, sixc)
or the geometry uses a non-standard detector arrangement (zaxis, s2d2, fivec).

    fourch  (horizontal scattering plane) matches Busing & Levy (1967) Fig. 1b.
    fourcv  (vertical   scattering plane) is the synchrotron convention.
    kappa4ch and kappa4cv follow the same convention as fourch and fourcv.

Kappa angle (alpha) convention
(Walko 2016 Fig. 3; Wyckoff 1985 Fig. 2(b); Thorkildsen 2006 Table 1;
Enraf-Nonius; ITC Vol. C Sec. 2.2.6):
The kappa axis is inclined by alpha degrees from the omega axis,
lying in the plane that contains both omega and the equivalent-
Eulerian chi axis, and tilted from omega toward that chi direction.
Per preset (with omega and chi-equivalent shown as physical
basis-direction lines, ignoring sign of handedness):

- kappa4cv: omega along transverse, chi-eq along vertical;
  kappa lies in the T-V plane between +T and +V.
- kappa4ch: omega along vertical, chi-eq along longitudinal;
  kappa lies in the V-L plane between +V and +L.
- kappa6c:  same as kappa4cv (mounted on top of mu and nu).

Typical value: alpha = 50 deg.

Handedness convention:
    These presets follow Walko (2016) and encode omega/kappa/phi/2theta as
    *left-handed* about their respective signed-axis vectors (e.g.
    ``-TRANSVERSE`` for kappa4cv komega).  ITC Vol. C Sec. 2.2.6.2 (2006)
    instead specifies the standard signs of omega/chi/phi as right-handed
    (only 2theta is left-handed in Hamilton's choice).  Either convention
    is internally consistent; users who prefer the ITC convention can
    construct their own geometry by negating the relevant Stage axis
    vectors.  See ``factories.py`` for how to build a custom geometry.

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
* ITC Vol. C, Sec. 2.2.6 (2006), p. 36
  DOI:10.1107/97809553602060000577 — kappa goniostat;
  cites Wyckoff (1985, p. 334) for the schematic picture and
  states "the κ axis is inclined at 50° to the ω axis".
* H.W. Wyckoff, Methods in Enzymology 114, 330-386 (1985) —
  kappa diffractometer, Fig. 2(b) on p. 334
* G. Thorkildsen, H.B. Larsen & J.A. Beukes,
  J. Appl. Cryst. 39, 151-157 (2006) — three-circle goniostat
  angle calculations, Table 1 (kappa axes)
* D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016) — kappa, zaxis, s2d2, fivec
"""

from __future__ import annotations

import inspect

from .diffractometer import AdHocDiffractometer
from .factories import BASIS_BL
from .factories import BASIS_YOU
from .factories import KAPPA_ALPHA_DEFAULT
from .factories import register_geometry
from .kappa import KappaPseudoAngleConvention
from .kappa import kappa_axis_from_eulerian
from .mode import REQUIRED
from .mode import BisectConstraint
from .mode import ConstraintSet
from .mode import DetectorConstraint
from .mode import ReferenceConstraint
from .mode import SampleConstraint
from .mode import VirtualBisectConstraint
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
def fourcv(basis: dict = BASIS_BL) -> AdHocDiffractometer:
    """
    Four-circle Eulerian diffractometer, vertical scattering plane (synchrotron).

    Walko (2016) designation: S3D1.

    Synchrotron configuration: omega and ttheta both rotate about the
    transverse axis, so the scattering plane is vertical.  This exploits the
    s-polarisation and tighter vertical collimation of synchrotron radiation
    (Walko 2016).

    Default basis: Busing & Levy (1967) — transverse=+x, longitudinal=+y, vertical=+z.

    Sample stack (floor first):
        omega  : transverse,   left-handed
        chi    : longitudinal, right-handed
        phi    : transverse,   left-handed

    Detector (floor, mechanically independent of sample stack):
        ttheta : transverse,   left-handed

    omega and ttheta share the same transverse axis; mechanically independent.

    References: W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
                D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
    """
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("omega", -TRANSVERSE, parent=None, role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -TRANSVERSE, parent="chi", role="sample"),
        Stage("ttheta", -TRANSVERSE, parent=None, role="detector"),
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
        "fixed_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["chi", "phi", "ttheta"],
        ),
        "fixed_psi": ConstraintSet(
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
            "(vertical scattering plane, transverse ttheta, synchrotron)"
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

    Default basis: Busing & Levy (1967) — transverse=+x, longitudinal=+y, vertical=+z.

    Sample stack (floor first):
        omega  : vertical,   left-handed
        chi    : transverse, right-handed
        phi    : vertical,   left-handed

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
        "fixed_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["chi", "phi", "ttheta"],
        ),
        "fixed_psi": ConstraintSet(
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
def psic(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """
    You (1999) '4S+2D' six-circle diffractometer (psic geometry).

    Walko (2016) designation: S4D2.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, transverse=+z.
    Detector axis: transverse.  Vertical scattering plane.

    Sample stack (floor first):
        mu  : vertical,     right-handed
        eta : transverse,   left-handed
        chi : longitudinal, right-handed
        phi : transverse,   left-handed

    Detector stack (floor first):
        nu    : vertical,   right-handed
        delta : transverse, left-handed

    mu and nu share the same vertical rotation axis; mechanically independent.

    Reference: H. You, J. Appl. Cryst. 32, 614-623 (1999).
               DOI: 10.1107/S0021889899001223
    """
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("eta", -TRANSVERSE, parent="mu", role="sample"),
        Stage("chi", +LONGITUDINAL, parent="eta", role="sample"),
        Stage("phi", -TRANSVERSE, parent="chi", role="sample"),
        Stage("nu", +VERTICAL, parent=None, role="detector"),
        Stage("delta", -TRANSVERSE, parent="nu", role="detector"),
    ]
    # psic: 6 DOF, N-3=3 constraints needed per mode (You 1999, S4D2).
    # Vertical bisect pair: eta(transverse) <-> delta(transverse)  => eta = delta/2
    # Horizontal bisect pair: mu(vertical) <-> nu(vertical)  => mu = nu/2
    modes = {
        # ── Vertical scattering plane ───────────────────────────────────────
        "bisecting_vertical": ConstraintSet(
            [
                BisectConstraint("eta", "delta"),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
        ),
        "fixed_phi_vertical": ConstraintSet(
            [
                SampleConstraint("phi", 0.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "delta"],
        ),
        "fixed_chi_vertical": ConstraintSet(
            [
                SampleConstraint("chi", 90.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "phi", "delta"],
        ),
        "fixed_alpha_i_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
                ReferenceConstraint("alpha_i", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_beta_out_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
                ReferenceConstraint("beta_out", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "alpha_eq_beta_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
                ReferenceConstraint("a_eq_b", True),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_psi_vertical": ConstraintSet(
            [
                BisectConstraint("eta", "delta"),
                SampleConstraint("mu", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
        # ── Horizontal scattering plane ─────────────────────────────────────
        "bisecting_horizontal": ConstraintSet(
            [
                BisectConstraint("mu", "nu"),
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
        ),
        "fixed_phi_horizontal": ConstraintSet(
            [
                SampleConstraint("phi", 0.0),
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "chi", "nu"],
        ),
        "fixed_chi_horizontal": ConstraintSet(
            [
                SampleConstraint("chi", 90.0),
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "phi", "nu"],
        ),
        "fixed_alpha_i_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
                ReferenceConstraint("alpha_i", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_beta_out_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
                ReferenceConstraint("beta_out", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "alpha_eq_beta_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
                ReferenceConstraint("a_eq_b", True),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_psi_horizontal": ConstraintSet(
            [
                BisectConstraint("mu", "nu"),
                SampleConstraint("eta", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
        # ── Lifting detector (out-of-plane) ─────────────────────────────────
        "lifting_detector_phi": ConstraintSet(
            [
                SampleConstraint("phi", 0.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("qaz", 90.0),
            ],
            computed=["phi", "nu", "delta"],
        ),
        "lifting_detector_mu": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("eta", 0.0),
                DetectorConstraint("qaz", 90.0),
            ],
            computed=["mu", "nu", "delta"],
        ),
    }
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=BASIS_YOU,
        description=(
            "You (1999) 4S+2D six-circle diffractometer "
            "(transverse detector, vertical scattering plane, synchrotron)"
        ),
        modes=modes,
        default_mode="bisecting_vertical",
    )


@register_geometry
def sixc(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """
    Lohmeier & Vlieg (1993) six-circle surface diffractometer (sixc geometry).

    Also known as the IUCr six-circle diffractometer.
    Walko (2016) designation: (S3D2)1.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, transverse=+z.

    Sample and detector stacks share the alpha (rotary table) base stage,
    making this a coupled geometry.  Useful for surface diffraction.

    From Fig. 1 and §2.1 of Lohmeier & Vlieg (1993):
    alpha and gamma rotate about the vertical axis (x in LV convention).
    omega, phi, and delta all rotate about the transverse axis (z in LV).
    chi rotates about the longitudinal axis (y in LV).

    Stack (floor first)::

        alpha (shared base): vertical,     right-handed  [rotary table]
          --> omega (sample):  transverse,   left-handed
                --> chi:       longitudinal, right-handed
                      --> phi: transverse,   left-handed
          --> delta (detector): transverse,  left-handed
                --> gamma:      vertical,    right-handed

    Reference: M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
    """
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("alpha", +VERTICAL, parent=None, role="sample"),
        Stage("omega", -TRANSVERSE, parent="alpha", role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -TRANSVERSE, parent="chi", role="sample"),
        Stage("delta", -TRANSVERSE, parent="alpha", role="detector"),
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
    The kappa axis lies in the vertical-transverse plane, tilted alpha degrees
    from the vertical toward the transverse axis.

    komega and ttheta both rotate about the transverse axis; the scattering
    plane is vertical (synchrotron convention).

    Default basis: Busing & Levy (1967) — transverse=+x, longitudinal=+y, vertical=+z.

    Sample stack (floor first):

    - ``komega`` — transverse, left-handed
    - ``kappa`` — tilted in the transverse-vertical plane, between +T
      and +V; α from +T toward +V (per Walko 2016 Fig. 3 and
      Thorkildsen 2006 Table 1).
    - ``kphi`` — transverse, left-handed

    Detector (floor, mechanically independent):

    - ``ttheta`` — transverse, left-handed

    komega and ttheta share the same transverse axis; mechanically independent.

    Handedness note: this preset encodes omega/kphi/2theta as left-handed
    about transverse, following Walko (2016).  ITC Vol. C Sec. 2.2.6.2
    (2006) prefers a right-handed sign convention for omega/chi/phi.
    The two conventions are equivalent up to motor-angle sign flips;
    either yields the same physical orientations.  See the module
    docstring for further discussion.

    Parameters
    ----------
    alpha_deg : float
        Kappa tilt angle in degrees (default 50).  Must be in (0, 90).

    References
    ----------
    * D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016) — Fig. 3.
    * G. Thorkildsen, H.B. Larsen & J.A. Beukes, J. Appl. Cryst. 39,
      151-157 (2006) — Table 1, eqn (3).
    * W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
    * ITC Vol. C, Sec. 2.2.6 (2006), p. 36 — α = 50°; cites
      Wyckoff (1985, p. 334) for the schematic picture.
    """
    TRANSVERSE = basis["transverse"]
    VERTICAL = basis["vertical"]
    # Kappa axis per Walko (2016) Fig. 3 and Thorkildsen et al. (2006)
    # Table 1: the kappa arm lies in the transverse-vertical plane,
    # tilted ``alpha_deg`` from the (unsigned) transverse direction
    # toward the (unsigned) vertical direction.  See issue #252.
    # Note: omega itself is left-handed about transverse (encoded as
    # ``-TRANSVERSE`` in the Stage line below); the kappa arm extends
    # into the +T+V quadrant regardless of the omega-handedness sign.
    kax = kappa_axis_from_eulerian(+TRANSVERSE, +VERTICAL, alpha_deg)
    stages = [
        Stage("komega", -TRANSVERSE, parent=None, role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -TRANSVERSE, parent="kappa", role="sample"),
        Stage("ttheta", -TRANSVERSE, parent=None, role="detector"),
    ]
    # kappa4cv: 4 DOF, N-3=1 constraint needed per mode.
    # Virtual Eulerian angles (omega, chi, phi) are mapped to / from
    # real kappa motors via the geometry-aware decomposition in
    # ad_hoc_diffractometer.kappa.  The ``bisecting`` mode uses
    # ``VirtualBisectConstraint`` to enforce *true* bisecting in the
    # virtual Eulerian frame: ``omega_virtual = ttheta/2``.  See
    # issues #226 and #241.
    convention = KappaPseudoAngleConvention(
        n_komega=-TRANSVERSE,
        n_kappa=kax,
        n_kphi=-TRANSVERSE,
        n_chi_eq=+VERTICAL,
    )
    modes = {
        "bisecting": ConstraintSet(
            [VirtualBisectConstraint("omega", "ttheta")],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_kphi": ConstraintSet(
            [SampleConstraint("kphi", 0.0)],
            computed=["komega", "kappa", "ttheta"],
        ),
        "fixed_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_psi": ConstraintSet(
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
        kappa_pseudo_angle_convention=convention,
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

    Default basis: Busing & Levy (1967) — transverse=+x, longitudinal=+y, vertical=+z.

    Sample stack (floor first):

    - ``komega`` — vertical, left-handed
    - ``kappa`` — tilted in the vertical-longitudinal plane, between +V
      and +L; α from +V toward +L (per Wyckoff 1985 Fig. 2(b)).
    - ``kphi`` — vertical, left-handed

    Detector (floor, mechanically independent):

    - ``ttheta`` — vertical, left-handed

    komega and ttheta share the same vertical axis; mechanically independent.

    Handedness note: this preset encodes omega/kphi/2theta as left-handed
    about vertical, following Walko (2016).  ITC Vol. C Sec. 2.2.6.2
    (2006) prefers a right-handed sign convention for omega/chi/phi.
    The two conventions are equivalent up to motor-angle sign flips;
    either yields the same physical orientations.  See the module
    docstring for further discussion.

    Parameters
    ----------
    alpha_deg : float
        Kappa tilt angle in degrees (default 50).  Must be in (0, 90).

    References
    ----------
    * H.W. Wyckoff, Methods in Enzymology 114, 330-386 (1985) —
      Fig. 2(b) on p. 334 (kappa diffractometer).
    * D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
    * W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
    * ITC Vol. C, Sec. 2.2.6 (2006), p. 36 — α = 50°; cites
      Wyckoff (1985, p. 334) for the schematic picture.
    """
    VERTICAL = basis["vertical"]
    LONGITUDINAL = basis["longitudinal"]
    # Kappa axis per Wyckoff (1985) Fig. 2(b): the kappa arm lies in
    # the vertical-longitudinal plane, tilted ``alpha_deg`` from the
    # (unsigned) vertical direction toward the (unsigned) longitudinal
    # direction (toward the X-ray source/sample).  See issue #252.
    # Note: omega itself is left-handed about vertical (encoded as
    # ``-VERTICAL`` in the Stage line below); the kappa arm extends
    # into the +V+L quadrant regardless of the omega-handedness sign.
    kax = kappa_axis_from_eulerian(+VERTICAL, +LONGITUDINAL, alpha_deg)
    stages = [
        Stage("komega", -VERTICAL, parent=None, role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -VERTICAL, parent="kappa", role="sample"),
        Stage("ttheta", -VERTICAL, parent=None, role="detector"),
    ]
    # kappa4ch: 4 DOF, N-3=1 constraint needed per mode.
    # Same mode set as kappa4cv; ``bisecting`` uses VirtualBisectConstraint
    # to enforce true virtual bisecting (omega_virtual = ttheta/2).
    convention = KappaPseudoAngleConvention(
        n_komega=-VERTICAL,
        n_kappa=kax,
        n_kphi=-VERTICAL,
        n_chi_eq=+LONGITUDINAL,
    )
    modes = {
        "bisecting": ConstraintSet(
            [VirtualBisectConstraint("omega", "ttheta")],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_kphi": ConstraintSet(
            [SampleConstraint("kphi", 0.0)],
            computed=["komega", "kappa", "ttheta"],
        ),
        "fixed_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_psi": ConstraintSet(
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
        kappa_pseudo_angle_convention=convention,
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

    See ESRF for a comprehensive example:
    https://www.esrf.fr/home/UsersAndScience/Experiments/CRG/BM02/equipment/diffractometer.html#c

    Extends the kappa4cv geometry with two additional axes (mu, nu) in
    the style of the psic geometry (You 1999), giving full orientation freedom.
    This is the synchrotron configuration with a transverse detector.

    Default basis: You (1999) — vertical=+x, longitudinal=+y, transverse=+z.

    Sample stack (floor first):

    - ``mu`` — vertical, right-handed (outermost)
    - ``komega`` — transverse, left-handed
    - ``kappa`` — tilted in the transverse-vertical plane, between +T
      and +V; α from +T toward +V (per Walko 2016 Fig. 3 and
      Thorkildsen 2006 Table 1).
    - ``kphi`` — transverse, left-handed

    Detector stack (floor first):

    - ``nu`` — vertical, right-handed
    - ``delta`` — transverse, left-handed

    mu and nu share the same vertical axis; mechanically independent.

    Handedness note: this preset encodes komega/kphi/delta as left-handed
    about transverse and mu/nu as right-handed about vertical, following
    Walko (2016) and You (1999) respectively.  ITC Vol. C Sec. 2.2.6.2
    (2006) prefers a right-handed sign convention for omega/chi/phi.
    The two conventions are equivalent up to motor-angle sign flips;
    either yields the same physical orientations.  See the module
    docstring for further discussion.

    Parameters
    ----------
    alpha_deg : float
        Kappa tilt angle in degrees (default 50).  Must be in (0, 90).

    References
    ----------
    * H.H. Sønsteby, D. Chernyshov, M. Getz, O. Nilsen & H. Fjellvåg,
      J. Synchrotron Rad. 20, 644-647 (2013) — six-axis κ
      diffractometer (KUMA6 at SNBL/ESRF).
    * G. Thorkildsen, H.B. Larsen & J.A. Beukes, J. Appl. Cryst. 39,
      151-157 (2006) — Table 1, eqn (3); extends to additional
      rotation axes (§3 last paragraph).
    * D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016) — Fig. 3.
    * H. You, J. Appl. Cryst. 32, 614-623 (1999) — outer-axis layout
      (psic 4S+2D) and You coordinate basis.
    * ITC Vol. C, Sec. 2.2.6 (2006), p. 36 — α = 50°; cites
      Wyckoff (1985, p. 334) for the schematic picture.
    """
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    # Kappa axis per Walko (2016) Fig. 3 and Thorkildsen et al. (2006)
    # Table 1 (same as kappa4cv): the kappa arm lies in the
    # transverse-vertical plane, tilted ``alpha_deg`` from the
    # (unsigned) transverse direction toward the (unsigned) vertical
    # direction.  kappa6c = kappa4cv sample stack mounted on top of
    # the You (1999) ``mu`` outer axis.  See issue #252.
    kax = kappa_axis_from_eulerian(+TRANSVERSE, +VERTICAL, alpha_deg)
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("komega", -TRANSVERSE, parent="mu", role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -TRANSVERSE, parent="kappa", role="sample"),
        Stage("nu", +VERTICAL, parent=None, role="detector"),
        Stage("delta", -TRANSVERSE, parent="nu", role="detector"),
    ]
    # kappa6c: 6 DOF, N-3=3 constraints needed per mode.
    # Vertical bisect pair (kappa-omega bisect): VirtualBisectConstraint
    #   enforces ``omega_virtual = delta/2`` — the physically correct
    #   bisecting condition for a kappa diffractometer (issues #226 and
    #   #241).
    # Horizontal bisect pair: mu(vertical) <-> nu(vertical) => mu = nu/2
    #   (literal motor bisect; mu is a real outer stage, not a kappa motor).
    convention = KappaPseudoAngleConvention(
        n_komega=-TRANSVERSE,
        n_kappa=kax,
        n_kphi=-TRANSVERSE,
        n_chi_eq=+VERTICAL,
    )
    modes = {
        # ── Implemented (generic solver) ────────────────────────────────────
        "bisecting_vertical": ConstraintSet(
            [
                VirtualBisectConstraint("omega", "delta"),
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
                VirtualBisectConstraint("omega", "delta"),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
        ),
        "fixed_nu": ConstraintSet(
            [
                DetectorConstraint("nu", 0.0),
                VirtualBisectConstraint("omega", "delta"),
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
        "fixed_psi_vertical": ConstraintSet(
            [
                VirtualBisectConstraint("omega", "delta"),
                SampleConstraint("mu", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "fixed_psi_horizontal": ConstraintSet(
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
            f"(transverse detector, synchrotron). "
            f"Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
        kappa_pseudo_angle_convention=convention,
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
    Default basis: You (1999) — vertical=+x, longitudinal=+y, transverse=+z.

    Designed for surface diffraction.  The sample surface normal is aligned
    parallel to the Z-axis, so the angle of incidence equals the alpha angle.
    The detector and sample both rotate about the shared alpha (base) axis.

    Stack (floor first)::

        alpha (shared base): vertical,     right-handed
          --> Z     (sample)  : longitudinal, right-handed
          --> delta (detector): transverse,   left-handed
                --> gamma :     vertical,     right-handed

    The total scattering angle is a compound of gamma, delta, and alpha
    (Walko 2016, eq. 17)::

        ttheta = arccos(cos(gamma)*cos(delta)*cos(alpha) + sin(alpha)*sin(gamma))

    References:
    J.M. Bloch, J. Appl. Cryst. 18, 33-36 (1985).
    D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016), eq. 17.
    """
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("alpha", +VERTICAL, parent=None, role="sample"),
        Stage("Z", +LONGITUDINAL, parent="alpha", role="sample"),
        Stage("delta", -TRANSVERSE, parent="alpha", role="detector"),
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
    Default basis: You (1999) — vertical=+x, longitudinal=+y, transverse=+z.

    Two independent sample axes (mu, Z) and two independent detector axes
    (nu, delta), all mechanically decoupled.  The angle of incidence is the
    mu angle; the surface normal is parallel to Z.

    Sample stack (floor first)::

        mu    : vertical,     right-handed
          --> Z : longitudinal, right-handed

    Detector stack (floor first)::

        nu    : vertical,     right-handed
          --> delta : transverse, left-handed

    mu and nu share the same vertical axis; mechanically independent.

    The total scattering angle is (Walko 2016, eq. 18)::

        ttheta = arccos(cos(nu) * cos(delta))

    References:
    K.W. Evans-Lutterodt & M.-T. Tang, J. Appl. Cryst. 28, 318-326 (1995).
    D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016), eq. 18.
    """
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("Z", +LONGITUDINAL, parent="mu", role="sample"),
        Stage("nu", +VERTICAL, parent=None, role="detector"),
        Stage("delta", -TRANSVERSE, parent="nu", role="detector"),
    ]
    # s2d2: 4 DOF, N-3=1 constraint needed per mode.
    modes = {
        "fixed_mu": ConstraintSet(
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
    Default basis: You (1999) — vertical=+x, longitudinal=+y, transverse=+z.

    A standard Eulerian four-circle (fourcv) is mounted on a fifth vertical
    rotation stage (mu) as a base.  The sample and detector motions are coupled
    through mu.  This provides an additional degree of freedom for accessing
    wider regions of reciprocal space, particularly at synchrotron sources.

    Stack (floor first)::

        mu (shared base): vertical,     right-handed
          --> omega (sample): transverse,   left-handed
                --> chi:      longitudinal, right-handed
                      --> phi: transverse,  left-handed
          --> ttheta (detector): transverse, left-handed

    References:
    E. Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987).
    D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
    """
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("omega", -TRANSVERSE, parent="mu", role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -TRANSVERSE, parent="chi", role="sample"),
        Stage("ttheta", -TRANSVERSE, parent="mu", role="detector"),
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
        "fixed_omega_noncoplanar": ConstraintSet(
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
