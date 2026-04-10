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

Writing a plugin factory
------------------------
Each factory accepts an optional ``basis`` keyword argument (defaulting to
the canonical convention for that geometry).  Inside the factory, resolve
physical-direction aliases locally::

    from ad_hoc_diffractometer import (
        AdHocDiffractometer, Stage, register_geometry, BASIS_YOU
    )

    @register_geometry
    def my_geometry(basis=BASIS_YOU):
        VERTICAL     = basis["vertical"]
        LATERAL      = basis["lateral"]
        LONGITUDINAL = basis["longitudinal"]
        stages = [
            Stage("omega", -LATERAL,      role="sample"),
            Stage("chi",   +LONGITUDINAL, parent="omega", role="sample"),
            Stage("phi",   -LATERAL,      parent="chi",   role="sample"),
            Stage("ttheta", -LATERAL,     role="detector"),
        ]
        return AdHocDiffractometer(
            name="my_geometry", stages=stages, basis=basis,
        )

The two public basis dicts ``BASIS_YOU`` and ``BASIS_BL`` are exported from
``ad_hoc_diffractometer``.  Pass ``basis=BASIS_BL`` for Busing & Levy
convention geometries.

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
BASIS_YOU = {
    "vertical": XHAT,
    "longitudinal": YHAT,
    "lateral": ZHAT,
}
#: Alias for backward compatibility.
_BASIS_YOU = BASIS_YOU

# Busing & Levy (1967) convention: xHat=lateral, yHat=longitudinal, zHat=vertical
BASIS_BL = {
    "lateral": np.array([1.0, 0.0, 0.0]),
    "longitudinal": np.array([0.0, 1.0, 0.0]),
    "vertical": np.array([0.0, 0.0, 1.0]),
}
#: Alias for backward compatibility.
_BASIS_BL = BASIS_BL


# ---------------------------------------------------------------------------
# Eulerian geometries
# ---------------------------------------------------------------------------


@register_geometry
def psic(basis: dict = _BASIS_YOU) -> AdHocDiffractometer:
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
def fourcv(basis: dict = _BASIS_BL) -> AdHocDiffractometer:
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
def fourch(basis: dict = _BASIS_BL) -> AdHocDiffractometer:
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
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(horizontal scattering plane, vertical ttheta, laboratory)"
        ),
    )


@register_geometry
def sixc(basis: dict = _BASIS_YOU) -> AdHocDiffractometer:
    """
    Lohmeier & Vlieg (1993) six-circle surface diffractometer (sixc geometry).

    Walko (2016) designation: (S3D2)1.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, lateral=+z.

    Sample and detector stacks share the alpha (rotary table) base stage,
    making this a coupled geometry.  Useful for surface diffraction.

    Stack (floor first):
        alpha (shared base): vertical,     right-handed  [rotary table]
          --> omega (sample):     longitudinal, right-handed
                --> chi:          longitudinal, right-handed
                      --> phi:    longitudinal, right-handed
          --> delta (detector):   lateral,      left-handed
                --> gamma:        vertical,     right-handed

    Reference: M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
    """
    VERTICAL = basis["vertical"]
    LATERAL = basis["lateral"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("alpha", +VERTICAL, parent=None, role="sample"),
        Stage("omega", +LONGITUDINAL, parent="alpha", role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", +LONGITUDINAL, parent="chi", role="sample"),
        Stage("delta", -LATERAL, parent="alpha", role="detector"),
        Stage("gamma", +VERTICAL, parent="delta", role="detector"),
    ]
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
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
def kappa4cv(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
    basis: dict = _BASIS_BL,
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
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            f"Four-circle kappa diffractometer, vertical scattering plane "
            f"(synchrotron). Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
    )


@register_geometry
def kappa4ch(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
    basis: dict = _BASIS_BL,
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
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            f"Four-circle kappa diffractometer, horizontal scattering plane "
            f"(laboratory). Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
    )


@register_geometry
def kappa6c(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
    basis: dict = _BASIS_YOU,
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
    )


# ---------------------------------------------------------------------------
# General-inclination geometries (Walko 2016, Sections 4.2 and 5)
# ---------------------------------------------------------------------------


@register_geometry
def zaxis(basis: dict = _BASIS_YOU) -> AdHocDiffractometer:
    """
    Z-axis four-circle diffractometer (general-inclination geometry).

    Walko (2016) designation: (S1D2)1.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, lateral=+z.

    Designed for surface diffraction.  The sample surface normal is aligned
    parallel to the Z-axis, so the angle of incidence equals the alpha angle.
    The detector and sample both rotate about the shared alpha (base) axis.

    Stack (floor first):
        alpha (shared base): vertical,     right-handed
          --> Z     (sample)  : longitudinal, right-handed
          --> delta (detector): lateral,      left-handed
                --> gamma :     vertical,     right-handed

    The total scattering angle is a compound of gamma, delta, and alpha
    (Walko 2016, eq. 17):
        ttheta = arccos(cos(gamma)*cos(delta)*cos(alpha) + sin(alpha)*sin(gamma))

    References: J.M. Bloch, J. Appl. Cryst. 18, 33-36 (1985).
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
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            "Z-axis four-circle diffractometer (Bloch 1985; Walko 2016 (S1D2)1). "
            "Surface normal parallel to Z-axis. "
            "Sample and detector share the alpha (base) stage."
        ),
    )


@register_geometry
def s2d2(basis: dict = _BASIS_YOU) -> AdHocDiffractometer:
    """
    S2D2 four-circle diffractometer (general-inclination geometry).

    Walko (2016) designation: S2D2.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, lateral=+z.

    Two independent sample axes (mu, Z) and two independent detector axes
    (nu, delta), all mechanically decoupled.  The angle of incidence is the
    mu angle; the surface normal is parallel to Z.

    Sample stack (floor first):
        mu    : vertical,     right-handed
          --> Z : longitudinal, right-handed

    Detector stack (floor first):
        nu    : vertical,     right-handed
          --> delta : lateral, left-handed

    mu and nu share the same vertical axis; mechanically independent.

    The total scattering angle is (Walko 2016, eq. 18):
        ttheta = arccos(cos(nu) * cos(delta))

    References: K.W. Evans-Lutterodt & M.-T. Tang, J. Appl. Cryst.
                    28, 318-326 (1995).
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
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            "S2D2 four-circle diffractometer (Evans-Lutterodt & Tang 1995; "
            "Walko 2016 S2D2). "
            "Fully decoupled sample (mu, Z) and detector (nu, delta) axes."
        ),
    )


@register_geometry
def fivec(basis: dict = _BASIS_YOU) -> AdHocDiffractometer:
    """
    Five-circle diffractometer (fourcv mounted on a vertical base).

    Walko (2016) designation: (S3D1)1.
    Default basis: You (1999) — vertical=+x, longitudinal=+y, lateral=+z.

    A standard Eulerian four-circle (fourcv) is mounted on a fifth vertical
    rotation stage (mu) as a base.  The sample and detector motions are coupled
    through mu.  This provides an additional degree of freedom for accessing
    wider regions of reciprocal space, particularly at synchrotron sources.

    Stack (floor first):
        mu (shared base): vertical,     right-handed
          --> omega (sample): lateral,      left-handed
                --> chi:      longitudinal, right-handed
                      --> phi: lateral,     left-handed
          --> ttheta (detector): lateral,   left-handed

    References: E. Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987).
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
    return AdHocDiffractometer(
        name=inspect.currentframe().f_code.co_name,
        stages=stages,
        basis=basis,
        description=(
            "Five-circle diffractometer: fourcv on vertical mu base "
            "(Vlieg et al. 1987; Walko 2016 (S3D1)1). "
            "Sample and detector coupled through mu."
        ),
    )
