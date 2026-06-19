# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
geometry_loader.py — declarative-YAML geometry loader (issue #267).

Parses YAML files (or in-memory YAML text) describing a diffractometer
geometry and returns a fully-configured
:class:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`
instance.

The schema is documented authoritatively in
``docs/source/reference/declarative_geometry_schema.md``.  This module
implements the schema; the how-to
``docs/source/howto/custom_geometry.md`` walks through writing a file
from scratch.

Public API
----------

- :func:`load_geometry_file` — construct a geometry from a path or YAML
  text (no registry mutation).
- :func:`register_geometry_file` — parse a YAML file from any path and
  add it to the geometry registry under its declared name (or an
  explicit ``name=`` override).
- :func:`register_geometry_yaml` — parse an in-memory YAML string and
  add it to the geometry registry under a caller-supplied name (the
  in-memory companion to :func:`register_geometry_file`).
- :data:`KIND_KEY` — the top-level marker key
  (``"ad_hoc_diffractometer_geometry"``).
- :data:`SUPPORTED_REVISIONS` — the schema revisions this loader
  understands.
- :data:`CURRENT_REVISION` — the schema revision newly-authored YAML
  files SHOULD declare.

Schema marker
-------------
Every declarative geometry file must declare itself at the top level::

    ad_hoc_diffractometer_geometry:
        schema_revision: 1

The presence of this key (with a well-formed mapping value containing a
supported integer ``schema_revision``) is the signal that the document
is a geometry declaration.  ``schema_revision`` is a fixed property of
the schema this file conforms to — not a per-file edit counter; users
should treat the value verbatim and only change it when migrating the
file to a different declarative-geometry schema revision.  Files
lacking the marker are treated by the polymorphic-string dispatch in
:func:`load_geometry_file` as paths, not as YAML text.

Polymorphic source dispatch
---------------------------
:func:`load_geometry_file` accepts either a filesystem path (as
``str``, :class:`pathlib.Path`, or any :class:`os.PathLike`) **or** a
``str`` containing YAML text.  The dispatch rule is:

1. ``Path`` / ``os.PathLike`` → always a path; ``FileNotFoundError`` if absent.
2. ``str`` parses as YAML *and* declares the
   :data:`KIND_KEY` marker → treat as YAML text.
3. ``str`` either fails to parse as YAML *or* parses but lacks the marker
   → treat as a filesystem path.  ``FileNotFoundError`` if the file does
   not exist; the error message names both attempts.

Strict unknown-key policy
-------------------------
Unknown top-level or nested keys are rejected with a ``ValueError`` that
names the offending key, the containing context, and the accepted-key
set for that context.  Future schema revisions may relax this policy.
"""

from __future__ import annotations

import logging
import os
import warnings
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from .axes import parse_axis
from .diffractometer import AdHocDiffractometer
from .factories import BASIS_BL
from .factories import BASIS_DEFAULT
from .factories import BASIS_YOU
from .factories import KAPPA_ALPHA_DEFAULT
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

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public schema constants
# ---------------------------------------------------------------------------

KIND_KEY: str = "ad_hoc_diffractometer_geometry"
"""Top-level marker key that identifies a YAML document as a declarative
diffractometer geometry."""

SUPPORTED_REVISIONS: frozenset[int] = frozenset({1})
"""Schema revisions this loader can parse."""

CURRENT_REVISION: int = 1
"""Schema revision that newly-authored YAML files SHOULD declare."""


def get_schema_text() -> str:
    """Return the JSON Schema (revision 1) as a string.

    The schema is shipped as ``ad_hoc_diffractometer/geometries/schema.json``
    and is intended for editor tooling, documentation, and any external
    consumer that wants a machine-readable description of the
    declarative geometry format.  The loader itself does **not** consume
    this schema at runtime — its hand-written validator enforces the
    same rules without taking on a ``jsonschema`` runtime dependency.
    """
    return (
        resources.files("ad_hoc_diffractometer.geometries")
        .joinpath("schema.json")
        .read_text()
    )


def get_schema() -> dict:
    """Return the parsed JSON Schema as a Python ``dict``."""
    import json

    return json.loads(get_schema_text())


_BASIS_SHORTHANDS: dict[str, dict[str, np.ndarray]] = {
    "BL": BASIS_BL,
    "YOU": BASIS_YOU,
    "DEFAULT": BASIS_DEFAULT,
}
"""Shorthand strings accepted as values of the top-level ``basis:`` key."""

# Accepted-key sets per schema context.  Used by the strict unknown-key
# policy to produce informative error messages.
_TOP_LEVEL_KEYS: frozenset[str] = frozenset(
    {
        KIND_KEY,
        "name",
        "documentation",
        "basis",
        "parameters",
        "kappa_chi_eq",
        "kappa_eulerian_chi",
        "stages",
        "modes",
    }
)

_MARKER_KEYS: frozenset[str] = frozenset({"schema_revision"})
_PARAMETERS_KEYS: frozenset[str] = frozenset({"alpha_deg"})
_STAGE_KEYS: frozenset[str] = frozenset({"name", "axis", "parent", "role"})
_MODE_KEYS: frozenset[str] = frozenset(
    {"default", "constraints", "computed", "extras", "cut_points"}
)

_CONSTRAINT_TYPES: frozenset[str] = frozenset(
    {"bisect", "virtual_bisect", "sample", "detector", "reference"}
)

_CONSTRAINT_KEYS: dict[str, frozenset[str]] = {
    "bisect": frozenset({"type", "stage1", "stage2"}),
    "virtual_bisect": frozenset({"type", "stage1", "stage2"}),
    "sample": frozenset({"type", "stage", "value"}),
    "detector": frozenset({"type", "stage", "value"}),
    "reference": frozenset({"type", "name", "value"}),
}


# ---------------------------------------------------------------------------
# Internal validation helpers
# ---------------------------------------------------------------------------


class GeometrySchemaError(ValueError):
    """Raised when a YAML document violates the declarative geometry schema."""


def _check_unknown_keys(
    mapping: dict[str, Any],
    accepted: frozenset[str],
    context: str,
) -> None:
    """Enforce the strict unknown-key policy for one mapping context.

    Parameters
    ----------
    mapping : dict
        The mapping to inspect.
    accepted : frozenset of str
        The keys allowed in ``mapping``.
    context : str
        Human-readable description of the location (e.g. ``"top-level"``,
        ``"inside ad_hoc_diffractometer_geometry"``,
        ``"inside stages[2]"``).
    """
    extras = set(mapping) - accepted
    if extras:
        raise GeometrySchemaError(
            f"unrecognised key(s) {sorted(extras)!r} {context}; "
            f"accepted keys: {sorted(accepted)!r}."
        )


def _validate_marker(doc: dict[str, Any]) -> int:
    """Validate the ``ad_hoc_diffractometer_geometry`` top-level marker.

    Returns
    -------
    int
        The schema revision declared by the document.

    Raises
    ------
    GeometrySchemaError
        If the marker is missing, malformed, or declares an unsupported
        revision.
    """
    if KIND_KEY not in doc:
        raise GeometrySchemaError(
            f"missing top-level marker key {KIND_KEY!r}; every declarative "
            f"geometry file must declare its kind and schema revision, e.g.\n"
            f"    {KIND_KEY}:\n"
            f"        schema_revision: {CURRENT_REVISION}"
        )
    marker = doc[KIND_KEY]
    if not isinstance(marker, dict):
        raise GeometrySchemaError(
            f"{KIND_KEY!r} value must be a mapping with a 'schema_revision' "
            f"integer; got {type(marker).__name__!r} ({marker!r})."
        )
    _check_unknown_keys(marker, _MARKER_KEYS, f"inside {KIND_KEY!r}")
    if "schema_revision" not in marker:
        raise GeometrySchemaError(
            f"{KIND_KEY!r} mapping must contain a 'schema_revision' integer key."
        )
    revision = marker["schema_revision"]
    if not isinstance(revision, int) or isinstance(revision, bool):
        raise GeometrySchemaError(
            f"{KIND_KEY!r}.schema_revision must be an integer; "
            f"got {type(revision).__name__!r} ({revision!r})."
        )
    if revision not in SUPPORTED_REVISIONS:
        raise GeometrySchemaError(
            f"unsupported geometry schema revision {revision}; this build "
            f"of ad_hoc_diffractometer supports {sorted(SUPPORTED_REVISIONS)!r}. "
            f"Upgrade the package to read newer revisions, or downgrade the "
            f"YAML to a supported revision."
        )
    return revision


def _resolve_basis(value: Any) -> dict[str, np.ndarray]:
    """Resolve a YAML ``basis`` value to a ``dict[str, ndarray]``.

    Accepts:

    - one of the shorthand strings ``"BL"``, ``"YOU"``, ``"DEFAULT"``
    - a mapping with the three physical-direction keys
      (``vertical``, ``longitudinal``, ``transverse``) whose values are
      either signed-axis strings (``"+x"``, ``"-vertical"``, ...) or
      length-3 numeric sequences.

    The returned dict is validated for unit length, mutual orthogonality,
    and right-handedness (``vertical × longitudinal ≈ transverse``).
    """
    if isinstance(value, str):
        key = value.strip().upper()
        if key not in _BASIS_SHORTHANDS:
            raise GeometrySchemaError(
                f"unknown basis shorthand {value!r}; "
                f"accepted: {sorted(_BASIS_SHORTHANDS)!r}, or supply an "
                f"explicit mapping with 'vertical', 'longitudinal', "
                f"'transverse' keys."
            )
        return {
            k: np.asarray(v, dtype=float).copy()
            for k, v in _BASIS_SHORTHANDS[key].items()
        }

    if not isinstance(value, dict):
        raise GeometrySchemaError(
            f"basis must be a shorthand string or a mapping; got "
            f"{type(value).__name__!r} ({value!r})."
        )
    required = {"vertical", "longitudinal", "transverse"}
    missing = required - set(value)
    if missing:
        raise GeometrySchemaError(
            f"basis mapping is missing required key(s) {sorted(missing)!r}; "
            f"need exactly {sorted(required)!r}."
        )
    extras = set(value) - required
    if extras:
        raise GeometrySchemaError(
            f"basis mapping has unrecognised key(s) {sorted(extras)!r}; "
            f"accepted: {sorted(required)!r}."
        )
    out: dict[str, np.ndarray] = {}
    for key in ("vertical", "longitudinal", "transverse"):
        item = value[key]
        if isinstance(item, str):
            try:
                vec = parse_axis(item)
            except ValueError as exc:
                raise GeometrySchemaError(f"basis['{key}']: {exc}") from None
        elif isinstance(item, list | tuple) and len(item) == 3:
            try:
                vec = np.asarray([float(c) for c in item], dtype=float)
            except (TypeError, ValueError) as exc:
                raise GeometrySchemaError(
                    f"basis['{key}']: numeric vector entries must be numbers; "
                    f"got {item!r} ({exc})."
                ) from None
        else:
            raise GeometrySchemaError(
                f"basis['{key}']: expected a signed-axis string "
                f"(e.g. '+x', '-vertical') or a length-3 numeric vector; "
                f"got {type(item).__name__!r} ({item!r})."
            )
        out[key] = vec
    _validate_basis(out)
    return out


def _validate_basis(basis: dict[str, np.ndarray]) -> None:
    """Verify a basis is orthonormal within tolerance.

    The three vectors must each have unit length and be mutually
    orthogonal.  Both literature conventions accepted by this package
    (Busing & Levy 1967 and You 1999) are orthonormal, but they assign
    physical names to Cartesian axes in mirror-image patterns:

    - You 1999: ``(vertical, longitudinal, transverse) = (X, Y, Z)``
      so ``vertical × longitudinal = transverse``.
    - Busing & Levy 1967:
      ``(transverse, longitudinal, vertical) = (X, Y, Z)`` so
      ``transverse × longitudinal = vertical`` and equivalently
      ``vertical × longitudinal = -transverse``.

    Both are physically meaningful and right-handed in their own
    canonical ordering; this check therefore enforces only orthonormality
    and not any particular chirality of the (vertical, longitudinal,
    transverse) ordering.
    """
    tol = 1e-9
    v = basis["vertical"]
    lon = basis["longitudinal"]
    t = basis["transverse"]
    for name, vec in (("vertical", v), ("longitudinal", lon), ("transverse", t)):
        norm = float(np.linalg.norm(vec))
        if abs(norm - 1.0) > tol:
            raise GeometrySchemaError(
                f"basis['{name}'] is not a unit vector (|v| = {norm:.6g}); "
                f"every basis vector must have unit length."
            )
    pairs = (
        ("vertical", "longitudinal", v, lon),
        ("longitudinal", "transverse", lon, t),
        ("vertical", "transverse", v, t),
    )
    for n1, n2, a, b in pairs:
        dot = float(np.dot(a, b))
        if abs(dot) > tol:
            raise GeometrySchemaError(
                f"basis['{n1}'] and basis['{n2}'] are not orthogonal "
                f"(dot product {dot:.6g})."
            )
    # Also require det = ±1 (degenerate / co-planar bases would be
    # caught by orthonormality; this is a defensive cross-check).
    det = float(np.linalg.det(np.column_stack([v, lon, t])))
    if abs(abs(det) - 1.0) > tol:
        raise GeometrySchemaError(
            f"basis is degenerate or non-orthonormal "
            f"(|det([v, l, t])| = {abs(det):.6g}, expected 1)."
        )


def _resolve_axis(
    spec: Any,
    basis: dict[str, np.ndarray],
    *,
    alpha_deg: float | None,
    kappa_chi_eq: np.ndarray | None,
    stage_name: str,
) -> np.ndarray:
    """Resolve one stage's ``axis:`` field to a numpy vector.

    Supported forms:

    - signed-axis string: ``'+transverse'``, ``'-vertical'``, ``'+x'``, etc.
    - length-3 numeric sequence: ``[0.0, 0.0, 1.0]``.
    - kappa-tilt mapping: ``{kappa_eulerian: [<ref>]}`` resolves to
      :func:`kappa_axis_from_eulerian` using the geometry's
      ``parameters.alpha_deg`` and the top-level ``kappa_chi_eq:``.  The
      single argument is the outer (komega) axis.
    """
    if isinstance(spec, str):
        try:
            return parse_axis(spec, basis=basis)
        except ValueError as exc:
            raise GeometrySchemaError(f"stage {stage_name!r}: {exc}") from None
    if isinstance(spec, list | tuple) and len(spec) == 3:
        try:
            return np.asarray([float(c) for c in spec], dtype=float)
        except (TypeError, ValueError) as exc:
            raise GeometrySchemaError(
                f"stage {stage_name!r}: numeric axis entries must be "
                f"numbers; got {spec!r} ({exc})."
            ) from None
    if isinstance(spec, dict) and "kappa_eulerian" in spec:
        if alpha_deg is None or kappa_chi_eq is None:
            raise GeometrySchemaError(
                f"stage {stage_name!r}: 'kappa_eulerian' axis form requires "
                f"the top-level 'parameters.alpha_deg' and 'kappa_chi_eq' "
                f"fields; one or both are missing."
            )
        extra = set(spec) - {"kappa_eulerian"}
        if extra:
            raise GeometrySchemaError(
                f"stage {stage_name!r}: kappa_eulerian axis spec accepts "
                f"only one key; got extra(s) {sorted(extra)!r}."
            )
        komega_spec = spec["kappa_eulerian"]
        if isinstance(komega_spec, str):
            n_komega = parse_axis(komega_spec, basis=basis)
        elif isinstance(komega_spec, list | tuple) and len(komega_spec) == 3:
            n_komega = np.asarray([float(c) for c in komega_spec], dtype=float)
        else:
            raise GeometrySchemaError(
                f"stage {stage_name!r}: 'kappa_eulerian' value must be a "
                f"signed-axis string or a length-3 numeric vector; got "
                f"{komega_spec!r}."
            )
        return kappa_axis_from_eulerian(n_komega, kappa_chi_eq, alpha_deg)
    raise GeometrySchemaError(
        f"stage {stage_name!r}: axis must be a signed-axis string, a "
        f"length-3 numeric vector, or a 'kappa_eulerian' mapping; got "
        f"{spec!r}."
    )


def _derive_kappa_eulerian_chi(
    n_komega: np.ndarray,
    basis: dict[str, np.ndarray],
    source_label: str,
) -> np.ndarray:
    """Derive the equivalent-Eulerian chi pseudo-angle axis when the
    YAML does not declare ``kappa_eulerian_chi`` explicitly (issue
    #284).

    Picks the first basis direction perpendicular to ``n_komega`` in
    the conventional order ``(+longitudinal, +vertical, +transverse)``.
    Every standard Eulerian preset shipped with this package
    (fourcv, fourch, psic, sixc, fivec) puts its ``chi`` rotation
    about ``+longitudinal``; honoring that order makes the kappa
    preset's equivalent-Eulerian decomposition match its sister
    Eulerian preset's chi pseudo-angle exactly.

    Parameters
    ----------
    n_komega : numpy.ndarray, shape (3,)
        Outer kappa stage axis.  Need not be normalized.
    basis : dict[str, numpy.ndarray]
        The geometry's basis dictionary with keys ``vertical``,
        ``longitudinal``, ``transverse``.
    source_label : str
        Label naming the source YAML file or string; used in error
        messages.

    Returns
    -------
    n_chi_eq : numpy.ndarray, shape (3,)
        A unit-magnitude basis direction perpendicular to
        ``n_komega``.

    Raises
    ------
    GeometrySchemaError
        If no basis direction is perpendicular to ``n_komega`` within
        tolerance ``1e-9`` — i.e. ``n_komega`` is not aligned to any
        single basis axis.  In that pathological case the YAML must
        declare ``kappa_eulerian_chi`` explicitly.
    """
    n_om = np.asarray(n_komega, dtype=float)
    n_om = n_om / np.linalg.norm(n_om)
    for label in ("longitudinal", "vertical", "transverse"):
        candidate = np.asarray(basis[label], dtype=float)
        candidate = candidate / np.linalg.norm(candidate)
        if abs(float(np.dot(n_om, candidate))) < 1e-9:
            return candidate
    raise GeometrySchemaError(
        f"{source_label}: cannot derive the equivalent-Eulerian chi "
        f"axis automatically because the outer kappa axis (komega = "
        f"{n_komega.tolist()!r}) is not perpendicular to any single "
        f"basis direction.  Declare 'kappa_eulerian_chi' explicitly "
        f"in the top-level YAML."
    )


def _resolve_extras(extras: dict[str, Any]) -> dict[str, Any]:
    """Map the literal string ``'REQUIRED'`` to the
    :data:`ad_hoc_diffractometer.mode.REQUIRED` sentinel.

    All other values pass through unchanged.  YAML ``null`` (Python
    ``None``) is preserved.
    """
    out: dict[str, Any] = {}
    for k, v in extras.items():
        if isinstance(v, str) and v == "REQUIRED":
            out[k] = REQUIRED
        else:
            out[k] = v
    return out


def _build_constraint(spec: Any, *, mode_name: str, index: int):
    """Translate one constraint spec dict into a constraint instance."""
    if not isinstance(spec, dict):
        raise GeometrySchemaError(
            f"mode {mode_name!r} constraints[{index}]: each constraint must "
            f"be a mapping; got {type(spec).__name__!r}."
        )
    ctype = spec.get("type")
    if ctype is None:
        raise GeometrySchemaError(
            f"mode {mode_name!r} constraints[{index}]: missing required "
            f"'type' key; one of {sorted(_CONSTRAINT_TYPES)!r}."
        )
    if ctype not in _CONSTRAINT_TYPES:
        raise GeometrySchemaError(
            f"mode {mode_name!r} constraints[{index}]: unknown constraint "
            f"type {ctype!r}; accepted: {sorted(_CONSTRAINT_TYPES)!r}."
        )
    accepted = _CONSTRAINT_KEYS[ctype]
    _check_unknown_keys(
        spec,
        accepted,
        f"in mode {mode_name!r} constraints[{index}] (type={ctype!r})",
    )
    if ctype == "bisect":
        for key in ("stage1", "stage2"):
            if key not in spec:
                raise GeometrySchemaError(
                    f"mode {mode_name!r} constraints[{index}]: bisect requires {key!r}."
                )
        return BisectConstraint(spec["stage1"], spec["stage2"])
    if ctype == "virtual_bisect":
        for key in ("stage1", "stage2"):
            if key not in spec:
                raise GeometrySchemaError(
                    f"mode {mode_name!r} constraints[{index}]: virtual_bisect "
                    f"requires {key!r}."
                )
        return VirtualBisectConstraint(spec["stage1"], spec["stage2"])
    if ctype == "sample":
        for key in ("stage", "value"):
            if key not in spec:
                raise GeometrySchemaError(
                    f"mode {mode_name!r} constraints[{index}]: sample requires {key!r}."
                )
        return SampleConstraint(spec["stage"], spec["value"])
    if ctype == "detector":
        for key in ("stage", "value"):
            if key not in spec:
                raise GeometrySchemaError(
                    f"mode {mode_name!r} constraints[{index}]: detector "
                    f"requires {key!r}."
                )
        return DetectorConstraint(spec["stage"], spec["value"])
    # ctype == "reference"
    for key in ("name", "value"):
        if key not in spec:
            raise GeometrySchemaError(
                f"mode {mode_name!r} constraints[{index}]: reference requires {key!r}."
            )
    return ReferenceConstraint(spec["name"], spec["value"])


def _build_constraint_set(name: str, spec: dict[str, Any]) -> ConstraintSet:
    """Translate one mode-spec mapping into a :class:`ConstraintSet`."""
    if not isinstance(spec, dict):
        raise GeometrySchemaError(
            f"mode {name!r}: must be a mapping; got {type(spec).__name__!r}."
        )
    _check_unknown_keys(spec, _MODE_KEYS, f"in mode {name!r}")
    constraints_spec = spec.get("constraints", [])
    if not isinstance(constraints_spec, list):
        raise GeometrySchemaError(
            f"mode {name!r}: 'constraints' must be a list; got "
            f"{type(constraints_spec).__name__!r}."
        )
    constraints = [
        _build_constraint(c, mode_name=name, index=i)
        for i, c in enumerate(constraints_spec)
    ]
    computed = spec.get("computed")
    if computed is not None and not isinstance(computed, list):
        raise GeometrySchemaError(
            f"mode {name!r}: 'computed' must be a list of stage names; "
            f"got {type(computed).__name__!r}."
        )
    extras_raw = spec.get("extras") or {}
    if not isinstance(extras_raw, dict):
        raise GeometrySchemaError(
            f"mode {name!r}: 'extras' must be a mapping; got "
            f"{type(extras_raw).__name__!r}."
        )
    cut_points = spec.get("cut_points") or {}
    if not isinstance(cut_points, dict):
        raise GeometrySchemaError(
            f"mode {name!r}: 'cut_points' must be a mapping; got "
            f"{type(cut_points).__name__!r}."
        )
    return ConstraintSet(
        constraints=constraints,
        computed=computed,
        extras=_resolve_extras(extras_raw),
        cut_points=dict(cut_points),
    )


# ---------------------------------------------------------------------------
# Polymorphic-source dispatch
# ---------------------------------------------------------------------------


def _is_geometry_yaml_text(text: str) -> tuple[bool, Any]:
    """Try to parse ``text`` as a YAML geometry declaration.

    Returns
    -------
    (is_geometry, parsed)
        ``is_geometry`` is ``True`` when the text parses successfully,
        produces a mapping, and contains the :data:`KIND_KEY` marker as
        a top-level key (regardless of the marker value's well-formedness).
        ``parsed`` is the parsed document (or ``None`` if parsing failed).

    Notes
    -----
    The marker value's well-formedness (mapping vs scalar, supported
    revision, etc.) does **not** affect this helper.  Once the
    :data:`KIND_KEY` key is present, the document has *declared* its
    intent; full schema validation then happens in
    :func:`_construct_from_doc`, which raises a :class:`GeometrySchemaError`
    on any malformed marker.  This means a user who typed an incorrect
    marker value (e.g. ``ad_hoc_diffractometer_geometry: 1``) gets a
    clear schema error rather than a misleading
    :class:`FileNotFoundError`.
    """
    try:
        parsed = yaml.safe_load(text)
    except yaml.YAMLError:
        return False, None
    if not isinstance(parsed, dict):
        return False, parsed
    if KIND_KEY not in parsed:
        return False, parsed
    return True, parsed


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


def _construct_from_doc(
    doc: dict[str, Any],
    *,
    source_label: str,
    overrides: dict[str, Any],
) -> AdHocDiffractometer:
    """Build an :class:`AdHocDiffractometer` from a parsed document.

    Parameters
    ----------
    doc : dict
        The YAML document (already parsed).  Must declare the marker.
    source_label : str
        Human-readable label used in error messages (typically the file
        path or ``"<in-memory YAML text>"``).
    overrides : dict
        Caller-supplied keyword overrides (currently ``basis`` and
        ``alpha_deg``).
    """
    if not isinstance(doc, dict):
        raise GeometrySchemaError(
            f"{source_label}: top-level YAML document must be a mapping; "
            f"got {type(doc).__name__!r}."
        )
    _validate_marker(doc)
    _check_unknown_keys(doc, _TOP_LEVEL_KEYS, "at top level")

    if "name" not in doc:
        raise GeometrySchemaError(
            f"{source_label}: missing required top-level key 'name'."
        )
    name = doc["name"]
    if not isinstance(name, str) or not name:
        raise GeometrySchemaError(
            f"{source_label}: 'name' must be a non-empty string; got {name!r}."
        )

    documentation = doc.get("documentation", "")
    if not isinstance(documentation, str):
        raise GeometrySchemaError(
            f"{source_label}: 'documentation' must be a string; got "
            f"{type(documentation).__name__!r}."
        )

    # parameters block (currently only alpha_deg)
    parameters = doc.get("parameters") or {}
    if not isinstance(parameters, dict):
        raise GeometrySchemaError(
            f"{source_label}: 'parameters' must be a mapping; got "
            f"{type(parameters).__name__!r}."
        )
    _check_unknown_keys(parameters, _PARAMETERS_KEYS, "in 'parameters'")

    alpha_deg: float | None = None
    if "alpha_deg" in parameters:
        alpha_deg = float(parameters["alpha_deg"])
    if "alpha_deg" in overrides and overrides["alpha_deg"] is not None:
        alpha_deg = float(overrides["alpha_deg"])
    elif alpha_deg is None and ("kappa_chi_eq" in doc or "kappa_eulerian_chi" in doc):
        # File declares a kappa pseudo-angle equivalent but no alpha; default.
        alpha_deg = KAPPA_ALPHA_DEFAULT

    # Basis resolution (precedence: caller kwarg > YAML 'basis:' > BASIS_DEFAULT)
    basis: dict[str, np.ndarray]
    if overrides.get("basis") is not None:
        ov = overrides["basis"]
        if isinstance(ov, dict):
            # caller-supplied dict; copy and validate
            basis = {k: np.asarray(v, dtype=float).copy() for k, v in ov.items()}
            required = {"vertical", "longitudinal", "transverse"}
            missing = required - set(basis)
            if missing:
                raise GeometrySchemaError(
                    f"{source_label}: caller-supplied basis is missing "
                    f"key(s) {sorted(missing)!r}."
                )
            _validate_basis(basis)
        else:
            basis = _resolve_basis(ov)
    elif "basis" in doc:
        basis = _resolve_basis(doc["basis"])
    else:
        warnings.warn(
            f"{source_label} does not declare a basis; using BASIS_DEFAULT "
            f"(vertical=YHAT, longitudinal=ZHAT, transverse=XHAT). Real "
            f"instruments should declare an explicit basis (e.g. "
            f"`basis: BL` or `basis: YOU`).",
            UserWarning,
            stacklevel=3,
        )
        basis = {k: np.asarray(v, dtype=float).copy() for k, v in BASIS_DEFAULT.items()}

    # Optional kappa_chi_eq for kappa-arm tilt-direction (Walko's
    # geometric formula in ``_resolve_axis``).  This direction lies in
    # the plane spanned by ``n_komega`` and ``n_kappa`` (it is the
    # in-plane perpendicular of ``n_komega``).  It is NOT the
    # equivalent-Eulerian chi-pseudo-angle axis (which is generally a
    # different direction; see ``kappa_eulerian_chi`` below and issue
    # #284).
    kappa_chi_eq: np.ndarray | None = None
    if "kappa_chi_eq" in doc:
        kappa_chi_eq_spec = doc["kappa_chi_eq"]
        if isinstance(kappa_chi_eq_spec, str):
            kappa_chi_eq = parse_axis(kappa_chi_eq_spec, basis=basis)
        elif (
            isinstance(kappa_chi_eq_spec, list | tuple) and len(kappa_chi_eq_spec) == 3
        ):
            kappa_chi_eq = np.asarray(
                [float(c) for c in kappa_chi_eq_spec], dtype=float
            )
        else:
            raise GeometrySchemaError(
                f"{source_label}: 'kappa_chi_eq' must be a signed-axis "
                f"string or a length-3 numeric vector; got "
                f"{kappa_chi_eq_spec!r}."
            )

    # Optional kappa_eulerian_chi for the equivalent-Eulerian chi
    # pseudo-angle axis (issue #284).  Distinct from ``kappa_chi_eq``:
    # this is the axis the kappa→Eulerian decomposition rotates about
    # for the virtual ``chi`` angle, and should match the corresponding
    # non-kappa Eulerian preset's ``chi`` axis (fourcv/fourch/psic chi
    # is conventionally ``+longitudinal``).  When omitted, the loader
    # derives it from the first basis direction perpendicular to
    # ``n_komega`` in the conventional order
    # (``+longitudinal``, ``+vertical``, ``+transverse``).
    kappa_eulerian_chi: np.ndarray | None = None
    if "kappa_eulerian_chi" in doc:
        spec = doc["kappa_eulerian_chi"]
        if isinstance(spec, str):
            kappa_eulerian_chi = parse_axis(spec, basis=basis)
        elif isinstance(spec, list | tuple) and len(spec) == 3:
            kappa_eulerian_chi = np.asarray([float(c) for c in spec], dtype=float)
        else:
            raise GeometrySchemaError(
                f"{source_label}: 'kappa_eulerian_chi' must be a "
                f"signed-axis string or a length-3 numeric vector; "
                f"got {spec!r}."
            )

    # Stages
    if "stages" not in doc:
        raise GeometrySchemaError(
            f"{source_label}: missing required top-level key 'stages'."
        )
    stages_spec = doc["stages"]
    if not isinstance(stages_spec, list) or not stages_spec:
        raise GeometrySchemaError(f"{source_label}: 'stages' must be a non-empty list.")
    stages: list[Stage] = []
    for i, stage_spec in enumerate(stages_spec):
        if not isinstance(stage_spec, dict):
            raise GeometrySchemaError(
                f"{source_label}: stages[{i}] must be a mapping; got "
                f"{type(stage_spec).__name__!r}."
            )
        _check_unknown_keys(stage_spec, _STAGE_KEYS, f"in stages[{i}]")
        for key in ("name", "axis", "parent", "role"):
            if key not in stage_spec:
                raise GeometrySchemaError(
                    f"{source_label}: stages[{i}] missing required key {key!r}."
                )
        sname = stage_spec["name"]
        axis = _resolve_axis(
            stage_spec["axis"],
            basis,
            alpha_deg=alpha_deg,
            kappa_chi_eq=kappa_chi_eq,
            stage_name=sname,
        )
        stages.append(
            Stage(
                name=sname,
                axis=axis,
                parent=stage_spec["parent"],
                role=stage_spec["role"],
            )
        )

    # Modes
    if "modes" not in doc:
        raise GeometrySchemaError(
            f"{source_label}: missing required top-level key 'modes'."
        )
    modes_spec = doc["modes"]
    if not isinstance(modes_spec, dict) or not modes_spec:
        raise GeometrySchemaError(
            f"{source_label}: 'modes' must be a non-empty mapping."
        )
    modes: dict[str, ConstraintSet] = {}
    default_mode: str | None = None
    for mname, mspec in modes_spec.items():
        if not isinstance(mname, str):
            raise GeometrySchemaError(
                f"{source_label}: mode names must be strings; got "
                f"{mname!r} ({type(mname).__name__!r})."
            )
        if isinstance(mspec, dict) and mspec.get("default") is True:
            if default_mode is not None:
                raise GeometrySchemaError(
                    f"{source_label}: more than one mode is marked "
                    f"'default: true' ({default_mode!r}, {mname!r}); "
                    f"exactly one default is allowed."
                )
            default_mode = mname
        modes[mname] = _build_constraint_set(mname, mspec)

    # Build the kappa pseudo-angle convention if applicable.
    kappa_convention: KappaPseudoAngleConvention | None = None
    if alpha_deg is not None and (
        kappa_chi_eq is not None or kappa_eulerian_chi is not None
    ):
        try:
            n_komega = next(s.axis for s in stages if s.name == "komega")
            n_kappa = next(s.axis for s in stages if s.name == "kappa")
            n_kphi = next(s.axis for s in stages if s.name == "kphi")
        except StopIteration:
            raise GeometrySchemaError(
                f"{source_label}: kappa parameters declared but the stage "
                f"list does not contain the canonical kappa stage names "
                f"'komega', 'kappa', 'kphi'.  The declarative loader "
                f"synthesizes the KappaPseudoAngleConvention from these "
                f"names; rename your stages accordingly or omit the "
                f"'kappa_chi_eq' / 'kappa_eulerian_chi' / "
                f"'parameters.alpha_deg' fields."
            ) from None

        # Resolve the equivalent-Eulerian chi axis (issue #284).
        # Precedence:
        #   1. explicit ``kappa_eulerian_chi`` field (caller override);
        #   2. derived: first basis direction perpendicular to
        #      ``n_komega`` in the conventional order
        #      ``(+longitudinal, +vertical, +transverse)``.
        # The conventional choice is ``+longitudinal``: every standard
        # 4-/6-circle Eulerian preset shipped with this package
        # (fourcv, fourch, psic, sixc, fivec) puts its ``chi`` rotation
        # about ``+longitudinal``.  Aligning the kappa equivalent-
        # Eulerian decomposition with that choice makes
        # ``forward()`` reachability of a kappa preset match its
        # non-kappa sister (fourcv↔kappa4cv, fourch↔kappa4ch,
        # psic↔kappa6c).
        if kappa_eulerian_chi is not None:
            n_chi_eq = kappa_eulerian_chi
        else:
            n_chi_eq = _derive_kappa_eulerian_chi(n_komega, basis, source_label)

        kappa_convention = KappaPseudoAngleConvention(
            n_komega=n_komega,
            n_kappa=n_kappa,
            n_kphi=n_kphi,
            n_chi_eq=n_chi_eq,
        )

    return AdHocDiffractometer(
        name=name,
        stages=stages,
        basis=basis,
        description=documentation,
        modes=modes,
        default_mode=default_mode,
        kappa_alpha_deg=alpha_deg,
        kappa_pseudo_angle_convention=kappa_convention,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_geometry_file(
    source: str | os.PathLike[str] | Path,
    **kwargs: Any,
) -> AdHocDiffractometer:
    """Construct a geometry from a path or YAML text.

    Parameters
    ----------
    source : str | os.PathLike | pathlib.Path
        Either a filesystem path to a YAML file, or a ``str`` containing
        YAML text.  See the dispatch rule in the module docstring.
    **kwargs
        Override values for ``basis`` and ``alpha_deg``.  Caller-supplied
        ``basis`` always wins over whatever the YAML declares.

    Returns
    -------
    AdHocDiffractometer

    Raises
    ------
    FileNotFoundError
        If the source resolves to a path that does not exist.
    GeometrySchemaError
        If the YAML document declares the geometry marker but is otherwise
        malformed.
    """
    accepted = {"basis", "alpha_deg"}
    extras = set(kwargs) - accepted
    if extras:
        raise TypeError(
            f"load_geometry_file: unrecognised keyword argument(s) "
            f"{sorted(extras)!r}; accepted: {sorted(accepted)!r}."
        )

    # Branch 1: explicit Path / os.PathLike → always a path.
    if isinstance(source, Path | os.PathLike) and not isinstance(source, str):
        path = Path(source)
        if not path.is_file():
            raise FileNotFoundError(
                f"declarative geometry file not found: {str(path)!r}."
            )
        return _construct_from_doc(
            yaml.safe_load(path.read_text()),
            source_label=str(path),
            overrides=kwargs,
        )

    # Branch 2: str — try as YAML first.
    if not isinstance(source, str):
        raise TypeError(
            f"load_geometry_file: 'source' must be a path or a str of YAML "
            f"text; got {type(source).__name__!r}."
        )

    is_geometry, parsed = _is_geometry_yaml_text(source)
    if is_geometry:
        return _construct_from_doc(
            parsed,
            source_label="<in-memory YAML text>",
            overrides=kwargs,
        )

    # Branch 3: fall back to filesystem path.
    path = Path(source)
    if not path.is_file():
        raise FileNotFoundError(
            f"declarative geometry source {source!r} not found.  "
            f"It was first attempted as YAML text and was not recognized "
            f"as a geometry declaration (no {KIND_KEY!r} marker), then "
            f"attempted as a filesystem path and the file does not exist."
        )
    return _construct_from_doc(
        yaml.safe_load(path.read_text()),
        source_label=str(path),
        overrides=kwargs,
    )


def register_geometry_file(
    path: str | os.PathLike[str] | Path,
    *,
    name: str | None = None,
) -> str:
    """Parse a YAML file and add it to the geometry registry.

    Parameters
    ----------
    path : str | os.PathLike | pathlib.Path
        Filesystem path to the YAML file.  Unlike
        :func:`load_geometry_file`, this function does **not** accept
        in-memory YAML text — registry entries need a stable name
        source, which a filesystem path provides naturally via the
        filename stem.
    name : str, optional
        Override for the registry key.  Defaults to the filename stem.

    Returns
    -------
    str
        The registry name under which the geometry was registered.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the registry already contains an entry under the same name.
    """
    # Local import to avoid a circular dependency at module load time.
    from . import factories as _factories

    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"register_geometry_file: file not found: {str(p)!r}.")
    doc = yaml.safe_load(p.read_text())
    # Validate eagerly so problems surface at registration time.
    geom = _construct_from_doc(doc, source_label=str(p), overrides={})
    final_name = name if name is not None else p.stem
    registry = _factories._GEOMETRY_REGISTRY  # noqa: SLF001
    if final_name in registry:
        existing = registry[final_name]
        raise ValueError(
            f"register_geometry_file: a geometry named {final_name!r} is "
            f"already registered ({existing!r}); pass an explicit name= "
            f"to register this file under a different name."
        )

    def _factory(**kwargs: Any) -> AdHocDiffractometer:
        # Re-read the file on each call so on-disk edits are picked up.
        fresh = yaml.safe_load(p.read_text())
        return _construct_from_doc(fresh, source_label=str(p), overrides=kwargs)

    _factory.__name__ = final_name
    _factory.__doc__ = geom.description.splitlines()[0] if geom.description else ""
    registry[final_name] = _factory
    return final_name


def register_geometry_yaml(
    yaml_text: str,
    *,
    name: str,
) -> str:
    """Parse an in-memory YAML string and add it to the geometry registry.

    The in-memory companion to :func:`register_geometry_file`.  Use it
    when the YAML definition of a geometry is already available as a
    Python ``str`` (for example, persisted inside another configuration
    document) and you want it discoverable via :func:`list_geometries`
    and :func:`make_geometry` without round-tripping through the
    filesystem.

    Parameters
    ----------
    yaml_text : str
        The YAML document declaring the geometry.  Must contain the
        :data:`KIND_KEY` marker with a supported
        :data:`SUPPORTED_REVISIONS` value.  Parsed eagerly so schema
        errors surface at registration time rather than at first
        :func:`make_geometry` call (matching
        :func:`register_geometry_file`).
    name : str
        Required registry key under which the geometry will be
        installed.  Unlike :func:`register_geometry_file` there is no
        filesystem path from which to derive a default, so this
        argument is mandatory.

    Returns
    -------
    str
        The registry name (the value passed in as ``name``).

    Raises
    ------
    ValueError
        If the registry already contains an entry under ``name``.
    GeometrySchemaError
        If ``yaml_text`` is not a well-formed declarative geometry
        document.

    Notes
    -----
    The registered factory re-parses ``yaml_text`` on every
    :func:`make_geometry` call so that each invocation returns a fresh
    :class:`AdHocDiffractometer` instance, mirroring the re-read
    semantics of :func:`register_geometry_file`.
    """
    # Local import to avoid a circular dependency at module load time.
    from . import factories as _factories

    source_label = f"<in-memory:{name}>"
    doc = yaml.safe_load(yaml_text)
    # Validate eagerly so problems surface at registration time.
    geom = _construct_from_doc(doc, source_label=source_label, overrides={})
    registry = _factories._GEOMETRY_REGISTRY  # noqa: SLF001
    if name in registry:
        existing = registry[name]
        raise ValueError(
            f"register_geometry_yaml: a geometry named {name!r} is "
            f"already registered ({existing!r}); pass a different "
            f"name= to register this YAML under another name."
        )

    # Capture yaml_text in the closure so each call re-parses the same
    # source.  This mirrors register_geometry_file's re-read semantics
    # (every make_geometry() returns a fresh AdHocDiffractometer).
    def _factory(**kwargs: Any) -> AdHocDiffractometer:
        fresh = yaml.safe_load(yaml_text)
        return _construct_from_doc(fresh, source_label=source_label, overrides=kwargs)

    _factory.__name__ = name
    _factory.__doc__ = geom.description.splitlines()[0] if geom.description else ""
    registry[name] = _factory
    return name


def _register_packaged_geometries() -> None:
    """Register every ``*.yml`` file shipped under
    :mod:`ad_hoc_diffractometer.geometries` with the geometry registry.

    Called once at package import time from
    :mod:`ad_hoc_diffractometer.factories`.  Subsequent calls are no-ops
    in the sense that already-registered names are skipped (with the
    same uniqueness guard as :func:`register_geometry_file`).
    """
    from . import factories as _factories

    pkg = "ad_hoc_diffractometer.geometries"
    try:
        files = resources.files(pkg)
    except (ImportError, ModuleNotFoundError):  # pragma: no cover
        logger.debug("packaged geometries directory %r not found", pkg)
        return
    registry = _factories._GEOMETRY_REGISTRY  # noqa: SLF001
    for entry in sorted(files.iterdir(), key=lambda e: e.name):
        if not entry.name.endswith(".yml"):
            continue
        stem = entry.name[: -len(".yml")]

        # Bind ``entry`` and ``stem`` into the closure explicitly so each
        # iteration captures its own values.
        def _make_factory(_entry=entry, _stem=stem):
            def _factory(**kwargs: Any) -> AdHocDiffractometer:
                text = _entry.read_text()
                fresh = yaml.safe_load(text)
                return _construct_from_doc(
                    fresh,
                    source_label=f"{pkg}/{_entry.name}",
                    overrides=kwargs,
                )

            _factory.__name__ = _stem
            return _factory

        factory = _make_factory()
        # Eagerly construct once to catch schema errors at import time.
        try:
            geom = factory()
        except Exception:  # pragma: no cover
            logger.exception(
                "failed to load packaged geometry %r; skipping", entry.name
            )
            continue
        factory.__doc__ = geom.description.splitlines()[0] if geom.description else ""
        if stem in registry:
            logger.debug(
                "declarative YAML overrides existing registration for %r",
                stem,
            )
        registry[stem] = factory
