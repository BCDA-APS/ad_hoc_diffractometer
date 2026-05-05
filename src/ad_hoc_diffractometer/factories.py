# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
factories.py — geometry registry and shared definitions.

This module provides the **geometry registry** infrastructure: the
``@register_geometry`` decorator, discovery of third-party geometry plugins
via entry points, and the ``list_geometries()`` / ``get_geometry()`` /
``make_geometry()`` lookup functions.

It also defines the shared **basis dictionaries** (``BASIS_YOU``,
``BASIS_BL``) and the default **kappa tilt angle** (``KAPPA_ALPHA_DEFAULT``)
used by the demo geometries in
:mod:`ad_hoc_diffractometer.presets`.

Demo geometries
---------------

The 10 demo geometries (``psic``, ``fourcv``, ``fourch``, ``sixc``,
``kappa4cv``, ``kappa4ch``, ``kappa6c``, ``zaxis``, ``s2d2``,
``fivec``) live in :mod:`ad_hoc_diffractometer.presets`.  Access them
as::

    import ad_hoc_diffractometer as ahd

    g = ahd.presets.fourcv()

Extension via entry points
--------------------------
Third-party packages can contribute additional geometries by
declaring an entry point in the ``"ad_hoc_diffractometer.geometries"``
group in their ``pyproject.toml``::

    [project.entry-points."ad_hoc_diffractometer.geometries"]
    my_geom = "my_package.module:my_geometry_function"

Each geometry function must accept no required arguments (it may accept
keyword arguments) and return an ``AdHocDiffractometer`` instance.
Entry-point geometries are discovered and loaded automatically when
``list_geometries()`` or ``get_geometry()`` is first called; they do NOT
need to call ``@register_geometry`` themselves.

**Geometry names must be globally unique.**  If an entry-point name
duplicates an already-registered name (whether a built-in or a
previously loaded plugin), a ``ValueError`` is raised at discovery time.
This prevents silent shadowing: an external package cannot overwrite
``fourcv``, ``psic``, or any other registered geometry.

Writing a custom geometry
-------------------------
Each geometry function accepts an optional ``basis`` keyword argument
(defaulting to the canonical convention for that geometry).  Inside the
function, resolve physical-direction aliases locally::

    from ad_hoc_diffractometer import AdHocDiffractometer, register_geometry
    from ad_hoc_diffractometer.factories import BASIS_YOU
    from ad_hoc_diffractometer.stage import Stage

    @register_geometry
    def my_geometry(basis=BASIS_YOU):
        VERTICAL     = basis["vertical"]
        TRANSVERSE   = basis["transverse"]
        LONGITUDINAL = basis["longitudinal"]
        stages = [
            Stage("omega", -TRANSVERSE,      role="sample"),
            Stage("chi",   +LONGITUDINAL, parent="omega", role="sample"),
            Stage("phi",   -TRANSVERSE,      parent="chi",   role="sample"),
            Stage("ttheta", -TRANSVERSE,     role="detector"),
        ]
        return AdHocDiffractometer(
            name="my_geometry", stages=stages, basis=basis,
        )

The two public basis dicts ``BASIS_YOU`` and ``BASIS_BL`` are available
from this module.  Pass ``basis=BASIS_BL`` for Busing & Levy convention
geometries.

References
----------
* W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967)
* H. You, J. Appl. Cryst. 32, 614-623 (1999)
  DOI:10.1107/S0021889899001223
* ITC Vol. C, Sec. 2.2.6 (2006)
  DOI:10.1107/97809553602060000577
* D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016)
"""

from __future__ import annotations

import logging
from importlib.metadata import entry_points

from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_GEOMETRY_REGISTRY: dict[str, type] = {}
"""Maps geometry name to a callable that returns a configured
``AdHocDiffractometer``.

Populated first by ``@register_geometry`` at import time, then supplemented
by installed third-party entry points the first time :func:`list_geometries`
or :func:`get_geometry` is called.
"""

_EP_LOADED: bool = False
"""``True`` once entry-point discovery has run, so it only runs once."""

GEOMETRY_ENTRY_POINT_GROUP = "ad_hoc_diffractometer.geometries"
"""Entry-point group name for geometry plugins."""


def get_geometry(name: str):
    """
    Return the registered geometry callable by name.

    This is the primitive lookup — it returns the callable that
    constructs the geometry, not an instance.  Use make_geometry()
    if you want an AdHocDiffractometer instance directly.

    Parameters
    ----------
    name : str
        Name of the geometry, as registered by @register_geometry
        (e.g. 'psic', 'fourcv', 'kappa4cv').

    Returns
    -------
    callable
        A callable that returns a configured AdHocDiffractometer
        for the named geometry.

    Raises
    ------
    ValueError
        If no geometry with that name is registered, with a message listing
        the available names.

    Examples
    --------
    >>> from ad_hoc_diffractometer import get_geometry
    >>> make_psic = get_geometry("psic")
    >>> make_psic()
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


def list_geometries() -> dict[str, type]:
    """
    Return a copy of the geometry registry as {name: callable}.

    Includes all built-in geometries (registered via ``@register_geometry``
    at import time) plus any third-party geometry plugins installed as
    entry points in the ``"ad_hoc_diffractometer.geometries"`` group.

    Entry-point discovery runs automatically the first time this function
    is called; subsequent calls use the already-populated registry.

    Returns
    -------
    dict
        Keys are geometry names (e.g. ``'psic'``, ``'fourcv'``).
        Values are callables that return a configured
        ``AdHocDiffractometer``.

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


def _load_entry_point_geometries() -> None:
    """
    Discover and load geometry callables from installed entry points.

    Scans the ``"ad_hoc_diffractometer.geometries"`` entry-point group
    for all installed packages (including this package itself) and adds
    any geometry callables not already present in ``_GEOMETRY_REGISTRY``.

    This function is called automatically — and only once — by
    ``list_geometries()`` and ``get_geometry()``.  It is idempotent:
    repeated calls after the first are no-ops.

    Notes
    -----
    Built-in demo geometries are registered via ``@register_geometry``
    at import time, so they are always present even if entry-point
    discovery fails.  Entry-point discovery supplements the registry
    with any third-party plugins that are installed but were not
    decorated with ``@register_geometry``.

    Each geometry name must be unique across all installed packages.  If
    an entry-point name collides with an already-registered name (whether
    a built-in or a previously loaded plugin), a ``ValueError`` is raised
    identifying the conflicting name and its source.  This prevents silent
    shadowing of built-in geometries and ambiguous duplicate registrations.

    If loading a particular entry point raises an exception *other than* a
    name collision (e.g. the plugin package is broken or missing), that
    entry point is silently skipped so that the rest of the registry is
    unaffected.

    Raises
    ------
    ValueError
        If an entry-point name duplicates an already-registered geometry
        name.
    """
    global _EP_LOADED  # noqa: PLW0603
    if _EP_LOADED:
        return
    _EP_LOADED = True

    # Register the packaged declarative-YAML geometries (issue #267)
    # before consulting third-party entry points, so plugin name
    # collisions surface against the canonical names.
    try:
        from .geometry_loader import _register_packaged_geometries

        _register_packaged_geometries()
    except Exception as exc:  # noqa: BLE001 — non-fatal at import time
        logger.debug("packaged geometry loader failed: %s", exc)

    try:
        eps = entry_points(group=GEOMETRY_ENTRY_POINT_GROUP)
        for ep in eps:
            try:
                factory = ep.load()
            except Exception as exc:  # noqa: BLE001
                logger.debug("Skipping broken entry point %r: %s", ep.name, exc)
                continue
            if ep.name in _GEOMETRY_REGISTRY:
                existing = _GEOMETRY_REGISTRY[ep.name]
                if existing is factory:
                    # Same callable re-declared via entry point (e.g. a built-in
                    # that is both @register_geometry'd and listed in pyproject.toml).
                    # This is not a conflict — skip silently.
                    continue
                raise ValueError(
                    f"Geometry name {ep.name!r} is already registered. "
                    f"Each geometry name must be unique across all installed "
                    f"packages. The entry point {ep.name!r} from {ep.value!r} "
                    f"conflicts with the existing registration "
                    f"{existing!r}. "
                    f"Rename the geometry in your package to resolve this conflict."
                )
            _GEOMETRY_REGISTRY[ep.name] = factory
    except ValueError:
        raise
    except Exception as exc:  # noqa: BLE001
        logger.debug("entry_points() failed; no plugins loaded: %s", exc)


def make_geometry(name: str, **kwargs):
    """
    Instantiate a geometry by name, passing keyword arguments to its
    constructor callable.

    Looks up the geometry by name via get_geometry() and calls it with
    the supplied kwargs.  This is the most convenient entry point for
    config-driven or programmatic geometry selection.

    Parameters
    ----------
    name : str
        Name of the geometry (e.g. 'psic', 'fourcv', 'kappa4cv').
    **kwargs
        Keyword arguments forwarded to the geometry constructor.  Most
        demo geometries take no arguments; kappa demo geometries accept
        alpha_deg.

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


def register_geometry(func):
    """
    Decorator that registers a geometry callable in _GEOMETRY_REGISTRY.

    The function is stored under its own ``__name__``, so the registry
    key is always identical to the callable's name.  The function is
    returned unchanged; this decorator has no runtime effect on the
    callable itself.

    Third-party packages do **not** need to use this decorator — they
    can instead declare an entry point in the
    ``"ad_hoc_diffractometer.geometries"`` group in their
    ``pyproject.toml`` and the geometry will be discovered automatically.

    Example
    -------
    ::

        @register_geometry
        def psic() -> AdHocDiffractometer:
            ...
    """
    _GEOMETRY_REGISTRY[func.__name__] = func
    return func


# ---------------------------------------------------------------------------
# Shared basis definitions
# ---------------------------------------------------------------------------

BASIS_YOU = {
    "vertical": XHAT,
    "longitudinal": YHAT,
    "transverse": ZHAT,
}
"""Basis vector dictionary for the You (1999) coordinate convention.

Maps physical direction names to Cartesian unit vectors:

- ``"vertical"`` → ``XHAT`` (+x, opposite to gravitational acceleration)
- ``"longitudinal"`` → ``YHAT`` (+y, along the beam)
- ``"transverse"`` → ``ZHAT`` (+z, completes the right-handed system: vertical × longitudinal)

Default basis used by ``psic``, ``sixc``, ``kappa6c``,
``zaxis``, ``s2d2``, and ``fivec`` (in :mod:`ad_hoc_diffractometer.presets`).
"""

BASIS_BL = {
    "vertical": ZHAT,
    "longitudinal": YHAT,
    "transverse": XHAT,
}
"""Basis vector dictionary for the Busing & Levy (1967) coordinate convention.

Maps physical direction names to Cartesian unit vectors:

- ``"transverse"`` → +x
- ``"longitudinal"`` → +y (along the beam)
- ``"vertical"`` → +z (opposite to gravitational acceleration)

Used by ``fourcv``, ``fourch``, ``kappa4cv``, and ``kappa4ch``
(in :mod:`ad_hoc_diffractometer.presets`).
"""

BASIS_DEFAULT = {
    "vertical": YHAT,
    "longitudinal": ZHAT,
    "transverse": XHAT,
}
"""Neutral basis used by the declarative geometry loader when a YAML file
omits the ``basis:`` key (issue #267).

Maps physical direction names to Cartesian unit vectors:

- ``"vertical"`` → ``YHAT`` (+y)
- ``"longitudinal"`` → ``ZHAT`` (+z)
- ``"transverse"`` → ``XHAT`` (+x)

This basis is **deliberately distinct** from both :data:`BASIS_YOU` (You
1999) and :data:`BASIS_BL` (Busing & Levy 1967) so that the package
does not appear to espouse either literature convention as a "project
default."  YAML files that should match a literature convention must
declare the basis explicitly (``basis: BL``, ``basis: YOU``, or an
explicit mapping); files that omit the key opt into this neutral
fallback.

The Cartesian basis ``(XHAT, YHAT, ZHAT)`` is right-handed; the
mapping above is one of three possible cyclic permutations of physical
names onto those axes.  Compare:

- :data:`BASIS_YOU` uses ``(vertical, longitudinal, transverse) = (X, Y, Z)``.
- :data:`BASIS_BL` uses ``(transverse, longitudinal, vertical) = (X, Y, Z)``.
- ``BASIS_DEFAULT`` uses ``(transverse, vertical, longitudinal) = (X, Y, Z)``,
  i.e. ``(vertical, longitudinal, transverse) = (Y, Z, X)`` — the
  remaining cyclic permutation.
"""

# ---------------------------------------------------------------------------
# Kappa angle default
# ---------------------------------------------------------------------------

KAPPA_ALPHA_DEFAULT = 50.0
"""Default kappa tilt angle in degrees (Walko 2016; Enraf-Nonius; ITC Vol. C)."""
