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
     - Vertical scattering plane (synchrotron).  Migrated to declarative
       YAML (``geometries/fourcv.yml``) in #267; the function in this
       module is a compatibility shim that delegates to the loader.
       See :doc:`/geometries/fourcv`.
   * - :func:`fourch`
     - Eulerian 4-circle
     - Horizontal scattering plane (laboratory).  Migrated to declarative
       YAML (``geometries/fourch.yml``) in #267; the function in this
       module is a compatibility shim that delegates to the loader.
       See :doc:`/geometries/fourch`.
   * - :func:`fivec`
     - Eulerian 5-circle
     - Vlieg et al. (1987), fourcv on a mu base.  Migrated to declarative
       YAML (``geometries/fivec.yml``) in #267; the function in this module
       is a compatibility shim that delegates to the loader.
       See :doc:`/geometries/fivec`.
   * - :func:`psic`
     - Eulerian 6-circle
     - You (1999) 4S+2D.  Migrated to declarative YAML
       (``geometries/psic.yml``) in #267; the function in this module
       is a compatibility shim that delegates to the loader.
       See :doc:`/geometries/psic`.
   * - :func:`sixc`
     - Eulerian 6-circle
     - Lohmeier & Vlieg (1993) surface geometry.  Migrated to declarative
       YAML (``geometries/sixc.yml``) in #267; the function in this module
       is a compatibility shim that delegates to the loader.
       See :doc:`/geometries/sixc`.
   * - :func:`kappa4cv`
     - Kappa 4-circle
     - Vertical scattering plane (synchrotron).  Migrated to declarative
       YAML (``geometries/kappa4cv.yml``) in #267; the function in this
       module is a compatibility shim that delegates to the loader.
       See :doc:`/geometries/kappa4cv`.
   * - :func:`kappa4ch`
     - Kappa 4-circle
     - Horizontal scattering plane (laboratory).  Migrated to declarative
       YAML (``geometries/kappa4ch.yml``) in #267; the function in this
       module is a compatibility shim that delegates to the loader.
       See :doc:`/geometries/kappa4ch`.
   * - :func:`kappa6c`
     - Kappa 6-circle
     - Psic-style outer axes.  Migrated to declarative YAML
       (``geometries/kappa6c.yml``) in #267; the function in this
       module is a compatibility shim that delegates to the loader.
       See :doc:`/geometries/kappa6c`.
   * - :func:`zaxis`
     - Surface / special
     - Bloch (1985) Z-axis geometry.  Migrated to declarative YAML
       (``geometries/zaxis.yml``) in #267; the function in this module
       is a compatibility shim that delegates to the loader.
       See :doc:`/geometries/zaxis`.
   * - :func:`s2d2`
     - Surface / special
     - Evans-Lutterodt & Tang (1995).  Migrated to declarative YAML
       (``geometries/s2d2.yml``) in #267; the function in this module
       is a compatibility shim that delegates to the loader.
       See :doc:`/geometries/s2d2`.

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

from importlib import resources

from .diffractometer import AdHocDiffractometer
from .factories import BASIS_BL
from .factories import BASIS_YOU
from .factories import KAPPA_ALPHA_DEFAULT

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
# ``fourcv`` has been migrated to a declarative YAML file
# (``geometries/fourcv.yml``) — see issue #267.  A compatibility shim
# named ``fourcv`` remains in this module to keep the legacy import
# path ``from ad_hoc_diffractometer.presets import fourcv`` working
# during the staged migration; it delegates to the loader.

# ---------------------------------------------------------------------------
# Eulerian geometries
# ---------------------------------------------------------------------------


# ``fourcv`` was the first demo geometry migrated from a Python factory
# to a declarative YAML file (``geometries/fourcv.yml``); see issue #267.
# This compatibility shim delegates to the loader so existing imports
# (``from ad_hoc_diffractometer.presets import fourcv``) continue to
# work during the staged migration.  The shim will be removed when
# ``presets.py`` is deleted at the end of the migration.
def fourcv(basis: dict = BASIS_BL) -> AdHocDiffractometer:
    """Return the declarative ``fourcv`` geometry (delegates to the YAML loader).

    See ``src/ad_hoc_diffractometer/geometries/fourcv.yml`` for the
    authoritative definition.
    """
    from .geometry_loader import load_geometry_file

    pkg_path = resources.files("ad_hoc_diffractometer.geometries").joinpath(
        "fourcv.yml"
    )
    # ``importlib.resources`` may return a MultiplexedPath / Traversable
    # rather than a real Path; round-trip via ``as_file`` to obtain a
    # filesystem path the loader can read.
    with resources.as_file(pkg_path) as p:
        return load_geometry_file(p, basis=basis)


# ``fourch`` was migrated to a declarative YAML file
# (``geometries/fourch.yml``) in issue #267.  This compatibility shim
# delegates to the loader so existing imports
# (``from ad_hoc_diffractometer.presets import fourch``) continue to
# work during the staged migration.  The shim will be removed when
# ``presets.py`` is deleted at the end of the migration.
def fourch(basis: dict = BASIS_BL) -> AdHocDiffractometer:
    """Return the declarative ``fourch`` geometry (delegates to the YAML loader).

    See ``src/ad_hoc_diffractometer/geometries/fourch.yml`` for the
    authoritative definition.
    """
    from .geometry_loader import load_geometry_file

    pkg_path = resources.files("ad_hoc_diffractometer.geometries").joinpath(
        "fourch.yml"
    )
    with resources.as_file(pkg_path) as p:
        return load_geometry_file(p, basis=basis)


# ``psic`` was migrated to a declarative YAML file
# (``geometries/psic.yml``) in issue #267.  This compatibility shim
# delegates to the loader so existing imports
# (``from ad_hoc_diffractometer.presets import psic``) continue to
# work during the staged migration.  The shim will be removed when
# ``presets.py`` is deleted at the end of the migration.
def psic(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """Return the declarative ``psic`` geometry (delegates to the YAML loader).

    See ``src/ad_hoc_diffractometer/geometries/psic.yml`` for the
    authoritative definition.
    """
    from .geometry_loader import load_geometry_file

    pkg_path = resources.files("ad_hoc_diffractometer.geometries").joinpath("psic.yml")
    with resources.as_file(pkg_path) as p:
        return load_geometry_file(p, basis=basis)


# ``sixc`` was migrated to a declarative YAML file
# (``geometries/sixc.yml``) in issue #267.  This compatibility shim
# delegates to the loader so existing imports
# (``from ad_hoc_diffractometer.presets import sixc``) continue to
# work during the staged migration.  The shim will be removed when
# ``presets.py`` is deleted at the end of the migration.
def sixc(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """Return the declarative ``sixc`` geometry (delegates to the YAML loader).

    See ``src/ad_hoc_diffractometer/geometries/sixc.yml`` for the
    authoritative definition.
    """
    from .geometry_loader import load_geometry_file

    pkg_path = resources.files("ad_hoc_diffractometer.geometries").joinpath("sixc.yml")
    with resources.as_file(pkg_path) as p:
        return load_geometry_file(p, basis=basis)


# ---------------------------------------------------------------------------
# Kappa geometries
# ---------------------------------------------------------------------------


# ``kappa4cv`` was migrated to a declarative YAML file
# (``geometries/kappa4cv.yml``) in issue #267.  This compatibility shim
# delegates to the loader so existing imports
# (``from ad_hoc_diffractometer.presets import kappa4cv``) continue to
# work during the staged migration.  The shim will be removed when
# ``presets.py`` is deleted at the end of the migration.
def kappa4cv(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
    basis: dict = BASIS_BL,
) -> AdHocDiffractometer:
    """Return the declarative ``kappa4cv`` geometry (delegates to the YAML loader).

    See ``src/ad_hoc_diffractometer/geometries/kappa4cv.yml`` for the
    authoritative definition.
    """
    from .geometry_loader import load_geometry_file

    pkg_path = resources.files("ad_hoc_diffractometer.geometries").joinpath(
        "kappa4cv.yml"
    )
    with resources.as_file(pkg_path) as p:
        return load_geometry_file(p, alpha_deg=alpha_deg, basis=basis)


# ``kappa4ch`` was migrated to a declarative YAML file
# (``geometries/kappa4ch.yml``) in issue #267.  This compatibility shim
# delegates to the loader so existing imports
# (``from ad_hoc_diffractometer.presets import kappa4ch``) continue to
# work during the staged migration.  The shim will be removed when
# ``presets.py`` is deleted at the end of the migration.
def kappa4ch(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
    basis: dict = BASIS_BL,
) -> AdHocDiffractometer:
    """Return the declarative ``kappa4ch`` geometry (delegates to the YAML loader).

    See ``src/ad_hoc_diffractometer/geometries/kappa4ch.yml`` for the
    authoritative definition.
    """
    from .geometry_loader import load_geometry_file

    pkg_path = resources.files("ad_hoc_diffractometer.geometries").joinpath(
        "kappa4ch.yml"
    )
    with resources.as_file(pkg_path) as p:
        return load_geometry_file(p, alpha_deg=alpha_deg, basis=basis)


# ``kappa6c`` was migrated to a declarative YAML file
# (``geometries/kappa6c.yml``) in issue #267.  This compatibility shim
# delegates to the loader so existing imports
# (``from ad_hoc_diffractometer.presets import kappa6c``) continue to
# work during the staged migration.  The shim will be removed when
# ``presets.py`` is deleted at the end of the migration.
def kappa6c(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
    basis: dict = BASIS_YOU,
) -> AdHocDiffractometer:
    """Return the declarative ``kappa6c`` geometry (delegates to the YAML loader).

    See ``src/ad_hoc_diffractometer/geometries/kappa6c.yml`` for the
    authoritative definition.
    """
    from .geometry_loader import load_geometry_file

    pkg_path = resources.files("ad_hoc_diffractometer.geometries").joinpath(
        "kappa6c.yml"
    )
    with resources.as_file(pkg_path) as p:
        return load_geometry_file(p, alpha_deg=alpha_deg, basis=basis)


# ---------------------------------------------------------------------------
# General-inclination geometries (Walko 2016, Sections 4.2 and 5)
# ---------------------------------------------------------------------------


# ``zaxis`` was migrated to a declarative YAML file
# (``geometries/zaxis.yml``) in issue #267.  This compatibility shim
# delegates to the loader so existing imports
# (``from ad_hoc_diffractometer.presets import zaxis``) continue to
# work during the staged migration.  The shim will be removed when
# ``presets.py`` is deleted at the end of the migration.
def zaxis(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """Return the declarative ``zaxis`` geometry (delegates to the YAML loader).

    See ``src/ad_hoc_diffractometer/geometries/zaxis.yml`` for the
    authoritative definition.
    """
    from .geometry_loader import load_geometry_file

    pkg_path = resources.files("ad_hoc_diffractometer.geometries").joinpath("zaxis.yml")
    with resources.as_file(pkg_path) as p:
        return load_geometry_file(p, basis=basis)


# ``s2d2`` was migrated to a declarative YAML file
# (``geometries/s2d2.yml``) in issue #267.  This compatibility shim
# delegates to the loader so existing imports
# (``from ad_hoc_diffractometer.presets import s2d2``) continue to
# work during the staged migration.  The shim will be removed when
# ``presets.py`` is deleted at the end of the migration.
def s2d2(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """Return the declarative ``s2d2`` geometry (delegates to the YAML loader).

    See ``src/ad_hoc_diffractometer/geometries/s2d2.yml`` for the
    authoritative definition.
    """
    from .geometry_loader import load_geometry_file

    pkg_path = resources.files("ad_hoc_diffractometer.geometries").joinpath("s2d2.yml")
    with resources.as_file(pkg_path) as p:
        return load_geometry_file(p, basis=basis)


# ``fivec`` was migrated to a declarative YAML file
# (``geometries/fivec.yml``) in issue #267.  This compatibility shim
# delegates to the loader so existing imports
# (``from ad_hoc_diffractometer.presets import fivec``) continue to
# work during the staged migration.  The shim will be removed when
# ``presets.py`` is deleted at the end of the migration.
def fivec(basis: dict = BASIS_YOU) -> AdHocDiffractometer:
    """Return the declarative ``fivec`` geometry (delegates to the YAML loader).

    See ``src/ad_hoc_diffractometer/geometries/fivec.yml`` for the
    authoritative definition.
    """
    from .geometry_loader import load_geometry_file

    pkg_path = resources.files("ad_hoc_diffractometer.geometries").joinpath("fivec.yml")
    with resources.as_file(pkg_path) as p:
        return load_geometry_file(p, basis=basis)
