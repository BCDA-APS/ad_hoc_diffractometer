"""
factories.py — predefined diffractometer geometry factory functions.

Each factory is decorated with @register_geometry, which registers it in
_GEOMETRY_REGISTRY under its function name.  Use list_geometries() to
retrieve all registered factories as a {name: callable} dict.

Extension via entry points
--------------------------
Third-party packages can contribute additional geometry factories by
declaring an entry point in the ``"ad_hoc_diffractometer.geometries"``
group in their ``pyproject.toml``::

    [project.entry-points."ad_hoc_diffractometer.geometries"]
    my_geom = "my_package.module:my_factory_function"

The factory function must accept no required arguments (it may accept
keyword arguments) and return an ``AdHocDiffractometer`` instance.
Entry-point geometries are discovered and loaded automatically when
``list_geometries()`` or ``get_geometry()`` is first called; they do NOT
need to call ``@register_geometry`` themselves.

The built-in geometries (psic, fourcv, fourch, sixc, kappa4cv, kappa4ch,
kappa6c, zaxis, s2d2, fivec) are also declared as entry points in the
package's own ``pyproject.toml`` under the same group, so they are
discoverable by any code that inspects installed entry points.

Naming convention:
    Eulerian geometries:     fourcv, fourch, sixc, psic
    Kappa geometries:        kappa4cv, kappa4ch, kappa6c
    Inclined geometries:     zaxis, s2d2, fivec

Scattering-plane suffix convention (v / h):
    v  ->  vertical   scattering plane  (synchrotron) — ttheta rotates about the lateral axis
    h  ->  horizontal scattering plane  (laboratory)  — ttheta rotates about the vertical axis

    Walko (2016): "at synchrotron sources the scattering plane is usually
    vertical, to take advantage of the (typically) s-polarization of the
    radiation and the higher degree of collimation in the vertical plane."

    Where no suffix is given, the geometry has no single-axis 2theta detector
    (inclined geometries) or the detector convention is unambiguous (psic, sixc).

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

References (chronological):
    W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967)               fourcv / fourch
    J.M. Bloch, J. Appl. Cryst. 18, 33-36 (1985)                          zaxis
    E. Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987)                   fivec
    M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993)            sixc
    K.W. Evans-Lutterodt & M.-T. Tang, J. Appl. Cryst. 28, 318-326 (1995) s2d2
    H. You, J. Appl. Cryst. 32, 614-623 (1999) DOI:10.1107/S0021889899001223   psic
    ITC Vol. C, Sec. 2.2.6 (2006) DOI:10.1107/97809553602060000577         kappa
    D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016)                 kappa, zaxis, s2d2, fivec
"""

from __future__ import annotations

import inspect
import logging
from importlib.metadata import entry_points

import numpy as np

from .axes import kappa_axis
from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT
from .geometry import AdHocDiffractometer
from .stage import Stage

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

#: Maps factory function name -> factory callable.
#: Populated first by @register_geometry at import time, then supplemented
#: by installed third-party entry points the first time list_geometries()
#: or get_geometry() is called.
_GEOMETRY_REGISTRY: dict[str, type] = {}

#: Set to True once entry-point discovery has run, so it only runs once.
_EP_LOADED: bool = False

#: Entry-point group name for geometry plugins.
GEOMETRY_ENTRY_POINT_GROUP = "ad_hoc_diffractometer.geometries"


def register_geometry(func):
    """
    Decorator that registers a geometry factory in _GEOMETRY_REGISTRY.

    The function is stored under its own ``__name__``, so the registry
    key is always identical to the callable's name.  The function is
    returned unchanged; this decorator has no runtime effect on the
    factory itself.

    Third-party packages do **not** need to use this decorator — they
    can instead declare an entry point in the
    ``"ad_hoc_diffractometer.geometries"`` group in their
    ``pyproject.toml`` and the factory will be discovered automatically.

    Example
    -------
    @register_geometry
    def psic() -> AdHocDiffractometer:
        ...
    """
    _GEOMETRY_REGISTRY[func.__name__] = func
    return func


def _load_entry_point_geometries() -> None:
    """
    Discover and load geometry factories from installed entry points.

    Scans the ``"ad_hoc_diffractometer.geometries"`` entry-point group
    for all installed packages (including this package itself) and adds
    any factories not already present in ``_GEOMETRY_REGISTRY``.

    This function is called automatically — and only once — by
    ``list_geometries()`` and ``get_geometry()``.  It is idempotent:
    repeated calls after the first are no-ops.

    Notes
    -----
    Built-in factories are registered via ``@register_geometry`` at
    import time, so they are always present even if entry-point discovery
    fails.  Entry-point discovery supplements the registry with any
    third-party plugins that are installed but were not decorated with
    ``@register_geometry``.

    If loading a particular entry point raises an exception (e.g. the
    plugin package is broken), that entry point is silently skipped so
    that the rest of the registry is unaffected.
    """
    global _EP_LOADED  # noqa: PLW0603
    if _EP_LOADED:
        return
    _EP_LOADED = True

    try:
        eps = entry_points(group=GEOMETRY_ENTRY_POINT_GROUP)
        for ep in eps:
            if ep.name not in _GEOMETRY_REGISTRY:
                try:
                    factory = ep.load()
                    _GEOMETRY_REGISTRY[ep.name] = factory
                except Exception as exc:  # noqa: BLE001
                    logger.debug("Skipping broken entry point %r: %s", ep.name, exc)
    except Exception as exc:  # noqa: BLE001
        logger.debug("entry_points() failed; no plugins loaded: %s", exc)


def list_geometries() -> dict[str, type]:
    """
    Return a copy of the geometry registry as {name: factory_callable}.

    Includes all built-in geometries (registered via ``@register_geometry``
    at import time) plus any third-party geometry plugins installed as
    entry points in the ``"ad_hoc_diffractometer.geometries"`` group.

    Entry-point discovery runs automatically the first time this function
    is called; subsequent calls use the already-populated registry.

    Returns
    -------
    dict
        Keys are factory names (e.g. ``'psic'``, ``'fourcv'``).
        Values are the callable factory functions.

    Examples
    --------
    >>> from ad_hoc_diffractometer import list_geometries
    >>> sorted(list_geometries())
    ['fivec', 'fourch', 'fourcv', 'kappa4ch', 'kappa4cv', 'kappa6c',
     'psic', 's2d2', 'sixc', 'zaxis']
    >>> list_geometries()['psic']()   # instantiate by name
    AdHocDiffractometer(name='psic', ...)
    """
    _load_entry_point_geometries()
    return dict(_GEOMETRY_REGISTRY)


def get_geometry(name: str):
    """
    Return the registered factory function for the named geometry.

    This is the primitive lookup — it returns the callable factory, not an
    instance.  Use make_geometry() if you want an AdHocDiffractometer
    instance directly.

    Parameters
    ----------
    name : str
        Name of the geometry, as registered by @register_geometry
        (e.g. 'psic', 'fourcv', 'kappa4cv').

    Returns
    -------
    callable
        The factory function for the named geometry.

    Raises
    ------
    ValueError
        If no geometry with that name is registered, with a message listing
        the available names.

    Examples
    --------
    >>> from ad_hoc_diffractometer import get_geometry
    >>> factory = get_geometry("psic")
    >>> factory()
    AdHocDiffractometer(name='psic', ...)
    >>> get_geometry("kappa4cv")(alpha_deg=50)
    AdHocDiffractometer(name='kappa4cv', ...)
    """
    _load_entry_point_geometries()
    if name not in _GEOMETRY_REGISTRY:
        available = sorted(_GEOMETRY_REGISTRY.keys())
        raise ValueError(
            f"No geometry named {name!r} is registered. "
            f"Available geometries: {available}."
        )
    return _GEOMETRY_REGISTRY[name]


def make_geometry(name: str, **kwargs) -> AdHocDiffractometer:
    """
    Instantiate a geometry by name, passing keyword arguments to its factory.

    Looks up the factory for the named geometry via get_geometry() and calls
    it with the supplied kwargs.  This is the most convenient entry point for
    config-driven or programmatic geometry selection.

    Parameters
    ----------
    name : str
        Name of the geometry (e.g. 'psic', 'fourcv', 'kappa4cv').
    **kwargs
        Keyword arguments forwarded to the factory function.  Most factories
        take no arguments; kappa factories accept alpha_deg.

    Returns
    -------
    AdHocDiffractometer
        A fully configured diffractometer geometry instance.

    Raises
    ------
    ValueError
        If no geometry with that name is registered.

    Examples
    --------
    >>> from ad_hoc_diffractometer import make_geometry
    >>> make_geometry("psic")
    AdHocDiffractometer(name='psic', ...)
    >>> make_geometry("kappa4cv", alpha_deg=50)
    AdHocDiffractometer(name='kappa4cv', ...)
    >>> make_geometry("kappa6c", alpha_deg=55)
    AdHocDiffractometer(name='kappa6c', ...)
    """
    return get_geometry(name)(**kwargs)


# ---------------------------------------------------------------------------
# Shared basis definitions
# ---------------------------------------------------------------------------

# You (1999) convention: xHat=vertical, yHat=longitudinal, zHat=lateral
_BASIS_YOU = {
    "vertical": XHAT,
    "longitudinal": YHAT,
    "lateral": ZHAT,
}

# Busing & Levy (1967) convention: xHat=lateral, yHat=longitudinal, zHat=vertical
_ZHAT_BL = np.array([0.0, 0.0, 1.0])  # vertical
_XHAT_BL = np.array([1.0, 0.0, 0.0])  # lateral
_YHAT_BL = np.array([0.0, 1.0, 0.0])  # longitudinal
_BASIS_BL = {
    "lateral": _XHAT_BL,
    "longitudinal": _YHAT_BL,
    "vertical": _ZHAT_BL,
}


# ---------------------------------------------------------------------------
# Eulerian geometries
# ---------------------------------------------------------------------------


@register_geometry
def psic() -> AdHocDiffractometer:
    """
    You (1999) '4S+2D' six-circle diffractometer (psic geometry).

    Walko (2016) designation: S4D2.
    Basis: xHat=vertical, yHat=longitudinal, zHat=lateral.
    Detector axis: lateral (-z).  Synchrotron / vertical scattering plane.

    Sample stack (floor first):
        mu  : vertical,     +x, right-handed
        eta : lateral,      -z, left-handed
        chi : longitudinal, +y, right-handed
        phi : lateral,      -z, left-handed

    Detector stack (floor first):
        nu    : vertical, +x, right-handed
        delta : lateral,  -z, left-handed

    mu and nu share the same vertical rotation axis; mechanically independent.

    Reference: H. You, J. Appl. Cryst. 32, 614-623 (1999).
               DOI: 10.1107/S0021889899001223
    """
    stages = [
        Stage("mu", +XHAT, parent=None, role="sample"),
        Stage("eta", -ZHAT, parent="mu", role="sample"),
        Stage("chi", +YHAT, parent="eta", role="sample"),
        Stage("phi", -ZHAT, parent="chi", role="sample"),
        Stage("nu", +XHAT, parent=None, role="detector"),
        Stage("delta", -ZHAT, parent="nu", role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=_BASIS_YOU,
        description=(
            "You (1999) 4S+2D six-circle diffractometer "
            "(lateral detector, vertical scattering plane, synchrotron)"
        ),
    )


@register_geometry
def fourcv() -> AdHocDiffractometer:
    """
    Four-circle Eulerian diffractometer, vertical scattering plane (synchrotron).

    Walko (2016) designation: S3D1.

    Synchrotron configuration: omega and ttheta both rotate about the
    lateral axis, so the scattering plane is vertical.  This exploits the
    s-polarisation and tighter vertical collimation of synchrotron radiation
    (Walko 2016).

    Basis (Busing & Levy convention):
        x = lateral, y = longitudinal, z = vertical.

    Sample stack (floor first):
        omega     : lateral,       -x, left-handed
        chi       : longitudinal,  +y, right-handed
        phi       : lateral,       -x, left-handed

    Detector (floor, mechanically independent of sample stack):
        ttheta : lateral,       -x, left-handed

    omega and ttheta share the same lateral axis; mechanically independent.

    References: W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
                D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
    """
    stages = [
        Stage("omega", -_BASIS_BL["lateral"], parent=None, role="sample"),
        Stage("chi", +_BASIS_BL["longitudinal"], parent="omega", role="sample"),
        Stage("phi", -_BASIS_BL["lateral"], parent="chi", role="sample"),
        Stage("ttheta", -_BASIS_BL["lateral"], parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=_BASIS_BL,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(vertical scattering plane, lateral ttheta, synchrotron)"
        ),
    )


@register_geometry
def fourch() -> AdHocDiffractometer:
    """
    Four-circle Eulerian diffractometer, horizontal scattering plane (laboratory).

    Walko (2016) designation: S3D1.

    Laboratory / default configuration: omega and ttheta both rotate about
    the vertical axis, so the scattering plane is horizontal.  This is the
    geometry described in Busing & Levy (1967), Fig. 1b.

    Basis (Busing & Levy convention):
        x = lateral, y = longitudinal, z = vertical.

    Sample stack (floor first):
        omega     : vertical, -z, left-handed
        chi       : lateral,  +x, right-handed
        phi       : vertical, -z, left-handed

    Detector (floor, mechanically independent):
        ttheta : vertical, -z, left-handed

    omega and ttheta share the same vertical axis; mechanically independent.

    Reference: W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
               D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
    """
    stages = [
        Stage("omega", -_BASIS_BL["vertical"], parent=None, role="sample"),
        Stage("chi", +_BASIS_BL["lateral"], parent="omega", role="sample"),
        Stage("phi", -_BASIS_BL["vertical"], parent="chi", role="sample"),
        Stage("ttheta", -_BASIS_BL["vertical"], parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=_BASIS_BL,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(horizontal scattering plane, vertical ttheta, laboratory)"
        ),
    )


@register_geometry
def sixc() -> AdHocDiffractometer:
    """
    Lohmeier & Vlieg (1993) six-circle surface diffractometer (sixc geometry).

    Walko (2016) designation: (S3D2)1.
    Basis: xHat=vertical, yHat=longitudinal (beam), zHat=lateral.

    Sample and detector stacks share the alpha (rotary table) base stage,
    making this a coupled geometry.  Useful for surface diffraction.

    Stack (floor first):
        alpha (shared base): vertical, +x, right-handed  [rotary table]
          --> omega (sample):     longitudinal, +y, right-handed
                --> chi:          longitudinal, +y, right-handed
                      --> phi:    longitudinal, +y, right-handed
          --> delta (detector):   lateral,      -z, left-handed
                --> gamma:        vertical,     +x, right-handed

    Reference: M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
    """
    stages = [
        Stage("alpha", +XHAT, parent=None, role="sample"),
        Stage("omega", +YHAT, parent="alpha", role="sample"),
        Stage("chi", +YHAT, parent="omega", role="sample"),
        Stage("phi", +YHAT, parent="chi", role="sample"),
        Stage("delta", -ZHAT, parent="alpha", role="detector"),
        Stage("gamma", +XHAT, parent="delta", role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=_BASIS_YOU,
        description=(
            "Lohmeier & Vlieg (1993) six-circle surface diffractometer "
            "(Walko S(3D2)1). "
            "Sample and detector share the alpha (rotary table) base stage."
        ),
    )


# ---------------------------------------------------------------------------
# Kappa geometries
# ---------------------------------------------------------------------------

#: Default kappa tilt angle in degrees (Walko 2016; Enraf-Nonius; ITC Vol. C).
KAPPA_ALPHA_DEFAULT = 50.0


@register_geometry
def kappa4cv(alpha_deg: float = KAPPA_ALPHA_DEFAULT) -> AdHocDiffractometer:
    """
    Four-circle kappa diffractometer, vertical scattering plane (synchrotron).

    Walko (2016) designation: S3D1.

    The chi circle of a standard Eulerian fourcv is replaced by a kappa arm.
    The kappa axis lies in the vertical-lateral plane, tilted alpha degrees
    from the vertical toward the lateral axis.

    komega and ttheta both rotate about the lateral axis; the scattering
    plane is vertical (synchrotron convention).

    Basis (Busing & Levy convention): x=lateral, y=longitudinal, z=vertical.

    Sample stack (floor first):
        komega    : lateral,  -x, left-handed
        kappa     : tilted,   kappa_axis(alpha), right-handed
        kphi      : lateral,  -x, left-handed

    Detector (floor, mechanically independent):
        ttheta : lateral,  -x, left-handed

    komega and ttheta share the same lateral axis; mechanically independent.

    Parameters
    ----------
    alpha_deg : float
        Kappa tilt angle in degrees (default 50).  Must be in (0, 90).

    References: D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
                W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
                ITC Vol. C, Sec. 2.2.6 (2006).
    """
    kax = kappa_axis(alpha_deg, basis=_BASIS_BL)
    stages = [
        Stage("komega", -_BASIS_BL["lateral"], parent=None, role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -_BASIS_BL["lateral"], parent="kappa", role="sample"),
        Stage("ttheta", -_BASIS_BL["lateral"], parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=_BASIS_BL,
        description=(
            f"Four-circle kappa diffractometer, vertical scattering plane "
            f"(synchrotron). Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
    )


@register_geometry
def kappa4ch(alpha_deg: float = KAPPA_ALPHA_DEFAULT) -> AdHocDiffractometer:
    """
    Four-circle kappa diffractometer, horizontal scattering plane (laboratory).

    Walko (2016) designation: S3D1.

    Identical to kappa4cv() but komega and ttheta rotate about the vertical
    axis, giving a horizontal scattering plane (laboratory convention).

    Basis (Busing & Levy convention): x=lateral, y=longitudinal, z=vertical.

    Sample stack (floor first):
        komega    : vertical, -z, left-handed
        kappa     : tilted,   kappa_axis(alpha), right-handed
        kphi      : vertical, -z, left-handed

    Detector (floor, mechanically independent):
        ttheta : vertical, -z, left-handed

    komega and ttheta share the same vertical axis; mechanically independent.

    Parameters
    ----------
    alpha_deg : float
        Kappa tilt angle in degrees (default 50).  Must be in (0, 90).

    References: D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
                W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
                ITC Vol. C, Sec. 2.2.6 (2006).
    """
    kax = kappa_axis(alpha_deg, basis=_BASIS_BL)
    stages = [
        Stage("komega", -_BASIS_BL["vertical"], parent=None, role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -_BASIS_BL["vertical"], parent="kappa", role="sample"),
        Stage("ttheta", -_BASIS_BL["vertical"], parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=_BASIS_BL,
        description=(
            f"Four-circle kappa diffractometer, horizontal scattering plane "
            f"(laboratory). Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
    )


@register_geometry
def kappa6c(alpha_deg: float = KAPPA_ALPHA_DEFAULT) -> AdHocDiffractometer:
    """
    Six-circle kappa diffractometer (psic-style outer axes, kappa inner sample).

    Walko (2016) designation: S4D2.

    Extends the kappa4cv geometry with two additional axes (mu, nu) in
    the style of the psic geometry (You 1999), giving full orientation freedom.
    This is the synchrotron configuration with a lateral detector.

    Basis: xHat=vertical, yHat=longitudinal, zHat=lateral (You 1999 convention).

    Sample stack (floor first):
        mu     : vertical,     +x, right-handed   [outermost]
        komega : lateral,      -z, left-handed
        kappa  : tilted,       kappa_axis(alpha), right-handed
        kphi   : lateral,      -z, left-handed

    Detector stack (floor first):
        nu     : vertical,     +x, right-handed
        delta  : lateral,      -z, left-handed

    mu and nu share the same vertical axis; mechanically independent.

    Parameters
    ----------
    alpha_deg : float
        Kappa tilt angle in degrees (default 50).  Must be in (0, 90).

    References: D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
                H. You, J. Appl. Cryst. 32, 614-623 (1999).
                ITC Vol. C, Sec. 2.2.6 (2006).
    """
    kax = kappa_axis(alpha_deg, basis=_BASIS_YOU)
    stages = [
        Stage("mu", +XHAT, parent=None, role="sample"),
        Stage("komega", -ZHAT, parent="mu", role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -ZHAT, parent="kappa", role="sample"),
        Stage("nu", +XHAT, parent=None, role="detector"),
        Stage("delta", -ZHAT, parent="nu", role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=_BASIS_YOU,
        description=(
            f"Six-circle kappa diffractometer, psic-style outer axes "
            f"(lateral detector, synchrotron). "
            f"Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
    )


# ---------------------------------------------------------------------------
# General-inclination geometries (Walko 2016, Sections 4.2 and 5)
# ---------------------------------------------------------------------------


@register_geometry
def zaxis() -> AdHocDiffractometer:
    """
    Z-axis four-circle diffractometer (general-inclination geometry).

    Walko (2016) designation: (S1D2)1.
    Basis: xHat=vertical, yHat=longitudinal, zHat=lateral.

    Designed for surface diffraction.  The sample surface normal is aligned
    parallel to the Z-axis, so the angle of incidence equals the alpha angle.
    The detector and sample both rotate about the shared alpha (base) axis.

    Stack (floor first):
        alpha (shared base): vertical, +x, right-handed
          --> Z     (sample)  : longitudinal, +y, right-handed
          --> delta (detector): lateral,      -z, left-handed
                --> gamma :     vertical,     +x, right-handed

    The total scattering angle is a compound of gamma, delta, and alpha
    (Walko 2016, eq. 17):
        2theta = arccos(cos(gamma)*cos(delta)*cos(alpha) + sin(alpha)*sin(gamma))

    References: J.M. Bloch, J. Appl. Cryst. 18, 33-36 (1985).
                D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016), eq. 17.
    """
    stages = [
        Stage("alpha", +XHAT, parent=None, role="sample"),
        Stage("Z", +YHAT, parent="alpha", role="sample"),
        Stage("delta", -ZHAT, parent="alpha", role="detector"),
        Stage("gamma", +XHAT, parent="delta", role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=_BASIS_YOU,
        description=(
            "Z-axis four-circle diffractometer (Bloch 1985; Walko 2016 (S1D2)1). "
            "Surface normal parallel to Z-axis. "
            "Sample and detector share the alpha (base) stage."
        ),
    )


@register_geometry
def s2d2() -> AdHocDiffractometer:
    """
    S2D2 four-circle diffractometer (general-inclination geometry).

    Walko (2016) designation: S2D2.
    Basis: xHat=vertical, yHat=longitudinal, zHat=lateral.

    Two independent sample axes (mu, Z) and two independent detector axes
    (nu, delta), all mechanically decoupled.  The angle of incidence is the
    mu angle; the surface normal is parallel to Z.

    Sample stack (floor first):
        mu    : vertical,     +x, right-handed
          --> Z : longitudinal, +y, right-handed

    Detector stack (floor first):
        nu    : vertical,     +x, right-handed
          --> delta : lateral, -z, left-handed

    mu and nu share the same vertical axis; mechanically independent.

    The total scattering angle is (Walko 2016, eq. 18):
        2theta = arccos(cos(nu) * cos(delta))

    References: K.W. Evans-Lutterodt & M.-T. Tang, J. Appl. Cryst.
                    28, 318-326 (1995).
                D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016), eq. 18.
    """
    stages = [
        Stage("mu", +XHAT, parent=None, role="sample"),
        Stage("Z", +YHAT, parent="mu", role="sample"),
        Stage("nu", +XHAT, parent=None, role="detector"),
        Stage("delta", -ZHAT, parent="nu", role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=_BASIS_YOU,
        description=(
            "S2D2 four-circle diffractometer (Evans-Lutterodt & Tang 1995; "
            "Walko 2016 S2D2). "
            "Fully decoupled sample (mu, Z) and detector (nu, delta) axes."
        ),
    )


@register_geometry
def fivec() -> AdHocDiffractometer:
    """
    Five-circle diffractometer (fourcv mounted on a vertical base).

    Walko (2016) designation: (S3D1)1.
    Basis: xHat=vertical, yHat=longitudinal, zHat=lateral.

    A standard Eulerian four-circle (fourcv) is mounted on a fifth vertical
    rotation stage (mu) as a base.  The sample and detector motions are coupled
    through mu.  This provides an additional degree of freedom for accessing
    wider regions of reciprocal space, particularly at synchrotron sources.

    Stack (floor first):
        mu (shared base): vertical, +x, right-handed
          --> omega (sample): lateral,  -z, left-handed
                --> chi:      lateral,  -z, left-handed (chi axis)
                      --> phi: lateral, -z, left-handed
          --> ttheta (detector): lateral, -z, left-handed

    Note: the inner fourc uses the You (1999) basis (xHat=vertical) rather
    than the Busing & Levy basis, since mu is the outermost vertical axis.

    References: E. Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987).
                D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
    """
    stages = [
        # Shared base
        Stage("mu", +XHAT, parent=None, role="sample"),
        # Sample stack on mu
        Stage("omega", -ZHAT, parent="mu", role="sample"),
        Stage("chi", +YHAT, parent="omega", role="sample"),
        Stage("phi", -ZHAT, parent="chi", role="sample"),
        # Detector stack on mu
        Stage("ttheta", -ZHAT, parent="mu", role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=_BASIS_YOU,
        description=(
            "Five-circle diffractometer: fourcv on vertical mu base "
            "(Vlieg et al. 1987; Walko 2016 (S3D1)1). "
            "Sample and detector coupled through mu."
        ),
    )
