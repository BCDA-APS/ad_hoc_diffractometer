"""
factories.py — predefined diffractometer geometry factory functions.

Each factory is decorated with @register_geometry, which registers it in
_GEOMETRY_REGISTRY under its function name.  Use list_geometries() to
retrieve all registered factories as a {name: callable} dict.

Naming convention:
    Eulerian geometries:     fourc_v, fourc_h, sixc, psic
    Kappa geometries:        kappa4c, kappa4c_h, kappa6c
    Inclined geometries:     zaxis, s2d2, fivec

Detector axis suffix convention (_v / _h):
    _v  ->  vertical detector axis   (laboratory / horizontal scattering plane)
    _h  ->  lateral   detector axis  (synchrotron / vertical scattering plane)

    Walko (2016): "at synchrotron sources the scattering plane is usually
    vertical, to take advantage of the (typically) s-polarization of the
    radiation and the higher degree of collimation in the vertical plane."

    Where no suffix is given, the geometry has no single-axis 2theta detector
    (inclined geometries) or the detector convention is unambiguous (psic, sixc).

    fourc_v (vertical detector) is the default/laboratory convention.
    kappa4c (vertical detector) is the default/laboratory convention.

Kappa angle (alpha) convention (Walko 2016; Enraf-Nonius; ITC Vol. C Sec. 2.2.6):
    The kappa axis lies in the vertical-lateral plane, tilted alpha degrees
    from the vertical axis toward the lateral axis.  Typical value: 50 deg.

Walko (2016) designations:
    S3D1      fourc_v, fourc_h, kappa4c, kappa4c_h
    S4D2      psic, kappa6c
    (S3D2)1   sixc   (sample and detector share base stage)
    (S3D1)1   fivec  (fourc_v mounted on vertical base)
    (S1D2)1   zaxis  (sample and detector share alpha base stage)
    S2D2      s2d2   (fully decoupled sample/detector pairs)

References (chronological):
    W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967)               fourc_v / fourc_h
    J.M. Bloch, J. Appl. Cryst. 18, 33-36 (1985)                          zaxis
    E. Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987)                   fivec
    M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993)            sixc
    K.W. Evans-Lutterodt & M.-T. Tang, J. Appl. Cryst. 28, 318-326 (1995) s2d2
    H. You, J. Appl. Cryst. 32, 614-623 (1999) DOI:10.1107/S0021889899001223   psic
    ITC Vol. C, Sec. 2.2.6 (2006) DOI:10.1107/97809553602060000577         kappa
    D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016)                 kappa, zaxis, s2d2, fivec
"""

from __future__ import annotations

import numpy as np

from .axes import kappa_axis
from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT
from .geometry import AdHocDiffractometer
from .stage import Stage

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

#: Maps factory function name -> factory callable.
#: Populated automatically by @register_geometry.
_GEOMETRY_REGISTRY: dict[str, type] = {}


def register_geometry(func):
    """
    Decorator that registers a geometry factory in _GEOMETRY_REGISTRY.

    The function is stored under its own __name__, so the registry key
    is always identical to the callable's name.  The function is returned
    unchanged; this decorator has no runtime effect on the factory itself.

    Example
    -------
    @register_geometry
    def psic() -> AdHocDiffractometer:
        ...
    """
    _GEOMETRY_REGISTRY[func.__name__] = func
    return func


def list_geometries() -> dict[str, type]:
    """
    Return a copy of the geometry registry as {name: factory_callable}.

    All entries were registered via @register_geometry at import time.

    Returns
    -------
    dict
        Keys are factory function names (e.g. 'psic', 'fourc_v', 'kappa4c').
        Values are the callable factory functions.

    Examples
    --------
    >>> from ad_hoc_diffractometer import list_geometries
    >>> list_geometries()
    {'psic': <function psic ...>, 'fourc_v': <function fourc_v ...>, ...}
    >>> list_geometries()['psic']()   # instantiate by name
    AdHocDiffractometer(name='psic', ...)
    """
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
        (e.g. 'psic', 'fourc_v', 'kappa4c').

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
    >>> get_geometry("kappa4c")(alpha_deg=50)
    AdHocDiffractometer(name='kappa4c', ...)
    """
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
        Name of the geometry (e.g. 'psic', 'fourc_v', 'kappa4c').
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
    >>> make_geometry("kappa4c", alpha_deg=50)
    AdHocDiffractometer(name='kappa4c', ...)
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
        name="psic",
        stages=stages,
        basis=_BASIS_YOU,
        description=(
            "You (1999) 4S+2D six-circle diffractometer "
            "(lateral detector, vertical scattering plane, synchrotron)"
        ),
    )


@register_geometry
def fourc_v() -> AdHocDiffractometer:
    """
    Busing & Levy (1967) four-circle Eulerian diffractometer, vertical detector.

    Walko (2016) designation: S3D1.

    This is the default / laboratory configuration: the detector arm swings
    in a vertical plane (horizontal scattering plane), which minimises
    gravitational distortions.  This is the geometry described in Busing &
    Levy (1967).

    Basis (Busing & Levy convention):
        x = lateral  (scattering vector at zero angles)
        y = longitudinal (along the beam)
        z = vertical
    Right-handed: lateral x longitudinal = vertical.

    Sample stack (floor first):
        omega     : vertical, -z, left-handed
        chi       : lateral,  +x, right-handed
        phi       : vertical, -z, left-handed

    Detector (floor, mechanically independent of sample stack):
        two_theta : vertical, -z, left-handed

    omega and two_theta share the same vertical axis; mechanically independent.

    Reference: W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
    """
    stages = [
        Stage("omega", -_ZHAT_BL, parent=None, role="sample"),
        Stage("chi", +_XHAT_BL, parent="omega", role="sample"),
        Stage("phi", -_ZHAT_BL, parent="chi", role="sample"),
        Stage("two_theta", -_ZHAT_BL, parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name="fourc_v",
        stages=stages,
        basis=_BASIS_BL,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(vertical detector, horizontal scattering plane, laboratory)"
        ),
    )


@register_geometry
def fourc_h() -> AdHocDiffractometer:
    """
    Busing & Levy (1967) four-circle Eulerian diffractometer, lateral detector.

    Walko (2016) designation: S3D1.

    Synchrotron configuration: the scattering plane is vertical (horizontal
    detector axis), exploiting the s-polarisation and tighter vertical
    collimation of synchrotron radiation (Walko 2016).

    Basis (Busing & Levy convention):
        x = lateral, y = longitudinal, z = vertical.

    Sample stack (floor first):
        omega     : vertical, -z, left-handed
        chi       : lateral,  +x, right-handed
        phi       : vertical, -z, left-handed

    Detector (floor, mechanically independent):
        two_theta : lateral,  -x, left-handed

    omega and two_theta share the same vertical axis; mechanically independent.

    Reference: W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
               D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016).
    """
    stages = [
        Stage("omega", -_ZHAT_BL, parent=None, role="sample"),
        Stage("chi", +_XHAT_BL, parent="omega", role="sample"),
        Stage("phi", -_ZHAT_BL, parent="chi", role="sample"),
        Stage("two_theta", -_XHAT_BL, parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name="fourc_h",
        stages=stages,
        basis=_BASIS_BL,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(lateral detector, vertical scattering plane, synchrotron)"
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
        name="sixc",
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
def kappa4c(alpha_deg: float = KAPPA_ALPHA_DEFAULT) -> AdHocDiffractometer:
    """
    Four-circle kappa diffractometer, vertical detector (laboratory default).

    Walko (2016) designation: S3D1.

    The chi circle of a standard Eulerian fourc_v is replaced by a kappa arm.
    The kappa axis lies in the vertical-lateral plane, tilted alpha degrees
    from the vertical toward the lateral axis.

    Basis (Busing & Levy convention): x=lateral, y=longitudinal, z=vertical.

    Sample stack (floor first):
        komega    : vertical, -z, left-handed
        kappa     : tilted,   kappa_axis(alpha), right-handed
        kphi      : vertical, -z, left-handed

    Detector (floor, mechanically independent):
        two_theta : vertical, -z, left-handed

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
        Stage("komega", -_ZHAT_BL, parent=None, role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -_ZHAT_BL, parent="kappa", role="sample"),
        Stage("two_theta", -_ZHAT_BL, parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name="kappa4c",
        stages=stages,
        basis=_BASIS_BL,
        description=(
            f"Four-circle kappa diffractometer, vertical detector (laboratory). "
            f"Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
    )


@register_geometry
def kappa4c_h(alpha_deg: float = KAPPA_ALPHA_DEFAULT) -> AdHocDiffractometer:
    """
    Four-circle kappa diffractometer, lateral detector (synchrotron).

    Walko (2016) designation: S3D1.

    Identical to kappa4c() but the two_theta (detector) axis is lateral,
    corresponding to the synchrotron configuration with a vertical
    scattering plane.

    Basis (Busing & Levy convention): x=lateral, y=longitudinal, z=vertical.

    Sample stack (floor first):
        komega    : vertical, -z, left-handed
        kappa     : tilted,   kappa_axis(alpha), right-handed
        kphi      : vertical, -z, left-handed

    Detector (floor, mechanically independent):
        two_theta : lateral,  -x, left-handed

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
        Stage("komega", -_ZHAT_BL, parent=None, role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -_ZHAT_BL, parent="kappa", role="sample"),
        Stage("two_theta", -_XHAT_BL, parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name="kappa4c_h",
        stages=stages,
        basis=_BASIS_BL,
        description=(
            f"Four-circle kappa diffractometer, lateral detector (synchrotron). "
            f"Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
    )


@register_geometry
def kappa6c(alpha_deg: float = KAPPA_ALPHA_DEFAULT) -> AdHocDiffractometer:
    """
    Six-circle kappa diffractometer (psic-style outer axes, kappa inner sample).

    Walko (2016) designation: S4D2.

    Extends the kappa4c geometry with two additional axes (mu, nu) in
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
        name="kappa6c",
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
        name="zaxis",
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
        name="s2d2",
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
    Five-circle diffractometer (fourc_v mounted on a vertical base).

    Walko (2016) designation: (S3D1)1.
    Basis: xHat=vertical, yHat=longitudinal, zHat=lateral.

    A standard Eulerian four-circle (fourc_v) is mounted on a fifth vertical
    rotation stage (mu) as a base.  The sample and detector motions are coupled
    through mu.  This provides an additional degree of freedom for accessing
    wider regions of reciprocal space, particularly at synchrotron sources.

    Stack (floor first):
        mu (shared base): vertical, +x, right-handed
          --> omega (sample): lateral,  -z, left-handed
                --> chi:      lateral,  -z, left-handed (chi axis)
                      --> phi: lateral, -z, left-handed
          --> two_theta (detector): lateral, -z, left-handed

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
        Stage("two_theta", -ZHAT, parent="mu", role="detector"),
    ]
    return AdHocDiffractometer(
        name="fivec",
        stages=stages,
        basis=_BASIS_YOU,
        description=(
            "Five-circle diffractometer: fourc_v on vertical mu base "
            "(Vlieg et al. 1987; Walko 2016 (S3D1)1). "
            "Sample and detector coupled through mu."
        ),
    )
