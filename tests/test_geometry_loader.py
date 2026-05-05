# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for :mod:`ad_hoc_diffractometer.geometry_loader`.

Covers the declarative-YAML schema (#267):

- the ``ad_hoc_diffractometer_geometry`` marker (presence, well-formedness,
  supported revision)
- the polymorphic ``source`` argument of :func:`load_geometry_file`
  (``Path``, ``os.PathLike``, ``str`` of YAML text, ``str`` of a path)
- the ``basis`` precedence (caller kwarg > YAML ``basis:`` > BASIS_DEFAULT)
- shorthand and explicit basis forms (``"BL"``, ``"YOU"``, mapping of
  signed-axis strings, mapping of numeric vectors)
- strict unknown-key policy at every nesting level
- :data:`REQUIRED` sentinel mapping in ``extras``
- constraint-type dispatch (``bisect``, ``virtual_bisect``, ``sample``,
  ``detector``, ``reference``)
- :func:`register_geometry_file` (registry insertion, name override,
  duplicate-name guard)
- module-level constants ``KIND_KEY``, ``SUPPORTED_REVISIONS``,
  ``CURRENT_REVISION``
"""

from __future__ import annotations

import os
import re

import numpy as np
import pytest

from ad_hoc_diffractometer import load_geometry_file
from ad_hoc_diffractometer import register_geometry_file
from ad_hoc_diffractometer.factories import _GEOMETRY_REGISTRY
from ad_hoc_diffractometer.factories import BASIS_BL
from ad_hoc_diffractometer.factories import BASIS_DEFAULT
from ad_hoc_diffractometer.factories import BASIS_YOU
from ad_hoc_diffractometer.geometry_loader import CURRENT_REVISION
from ad_hoc_diffractometer.geometry_loader import KIND_KEY
from ad_hoc_diffractometer.geometry_loader import SUPPORTED_REVISIONS
from ad_hoc_diffractometer.geometry_loader import GeometrySchemaError
from ad_hoc_diffractometer.geometry_loader import get_schema
from ad_hoc_diffractometer.geometry_loader import get_schema_text
from ad_hoc_diffractometer.mode import REQUIRED
from ad_hoc_diffractometer.mode import BisectConstraint
from ad_hoc_diffractometer.mode import DetectorConstraint
from ad_hoc_diffractometer.mode import ReferenceConstraint
from ad_hoc_diffractometer.mode import SampleConstraint

# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------

# Minimal complete YAML for a 1-stage sample + 1-stage detector geometry
# (passes ConstraintSet validation and the demo loader without needing
# the full fourcv complexity).
_MINIMAL_YAML = """\
ad_hoc_diffractometer_geometry:
    revision: 1
name: tiny
documentation: A minimal demo geometry used in loader unit tests.
basis: BL
stages:
    - {name: omega,  axis: -transverse, parent: null,  role: sample}
    - {name: ttheta, axis: -transverse, parent: null,  role: detector}
modes:
    bisecting:
        default: true
        constraints:
            - {type: bisect, stage1: omega, stage2: ttheta}
        computed: [omega, ttheta]
"""


def _yaml_with(extra: str, base: str = _MINIMAL_YAML) -> str:
    """Return ``base`` with ``extra`` appended at the end (for adding keys)."""
    return base + extra


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------


def test_kind_key_value():
    assert KIND_KEY == "ad_hoc_diffractometer_geometry"


def test_current_revision_in_supported_set():
    assert CURRENT_REVISION in SUPPORTED_REVISIONS


def test_supported_revisions_is_frozenset_of_int():
    assert isinstance(SUPPORTED_REVISIONS, frozenset)
    for r in SUPPORTED_REVISIONS:
        assert isinstance(r, int) and not isinstance(r, bool)


# ---------------------------------------------------------------------------
# Polymorphic ``source`` dispatch
# ---------------------------------------------------------------------------


def test_load_from_pathlib_path(tmp_path):
    p = tmp_path / "tiny.yml"
    p.write_text(_MINIMAL_YAML)
    g = load_geometry_file(p)
    assert g.name == "tiny"


def test_load_from_str_path(tmp_path):
    p = tmp_path / "tiny.yml"
    p.write_text(_MINIMAL_YAML)
    g = load_geometry_file(str(p))
    assert g.name == "tiny"


def test_load_from_oslike_path(tmp_path):
    class Wrap(os.PathLike):
        def __init__(self, p):
            self._p = p

        def __fspath__(self):
            return str(self._p)

    p = tmp_path / "tiny.yml"
    p.write_text(_MINIMAL_YAML)
    g = load_geometry_file(Wrap(p))
    assert g.name == "tiny"


def test_load_from_yaml_text():
    g = load_geometry_file(_MINIMAL_YAML)
    assert g.name == "tiny"


def test_str_without_marker_falls_back_to_path(tmp_path):
    """A str that doesn't validate as a geometry is tried as a path."""
    not_a_file = str(tmp_path / "nonexistent.yml")
    with pytest.raises(FileNotFoundError, match=re.escape(not_a_file)):
        load_geometry_file(not_a_file)


def test_str_yaml_without_marker_falls_back_to_path(tmp_path):
    """Valid YAML *without* the geometry marker is tried as a path.

    The combined diagnostic in the FileNotFoundError message names both
    attempts.
    """
    yaml_no_marker = "name: oops\nstages: []\n"
    with pytest.raises(FileNotFoundError, match=KIND_KEY):
        load_geometry_file(yaml_no_marker)


def test_pathlib_path_missing_raises_filenotfound(tmp_path):
    p = tmp_path / "missing.yml"
    with pytest.raises(FileNotFoundError, match="missing.yml"):
        load_geometry_file(p)


def test_unsupported_source_type_raises_typeerror():
    with pytest.raises(TypeError, match="must be a path or a str"):
        load_geometry_file(12345)


# ---------------------------------------------------------------------------
# Marker validation
# ---------------------------------------------------------------------------


def test_missing_marker_in_loaded_text_falls_back_to_path():
    """A document parsed but missing the marker is treated as a path."""
    text = "name: foo\nstages: []\nmodes: {}\n"
    with pytest.raises(FileNotFoundError):
        load_geometry_file(text)


@pytest.mark.parametrize(
    "marker_value, match_pattern",
    [
        pytest.param(
            "scalar",
            r"value must be a mapping",
            id="marker-not-a-mapping",
        ),
        pytest.param(
            {"future_key": 1},
            r"unrecognised key.*future_key",
            id="unknown-marker-subkey",
        ),
        pytest.param(
            {"revision": "1"},
            r"revision must be an integer",
            id="revision-string",
        ),
        pytest.param(
            {"revision": True},
            r"revision must be an integer",
            id="revision-bool",
        ),
        pytest.param(
            {"revision": 999},
            r"unsupported geometry schema revision 999",
            id="revision-unsupported",
        ),
        pytest.param(
            {},
            r"must contain a 'revision' integer key",
            id="revision-missing",
        ),
    ],
)
def test_malformed_marker_raises(marker_value, match_pattern):
    """Malformed markers raise *without* falling back to a path attempt."""
    import yaml as _yaml

    doc = {
        KIND_KEY: marker_value,
        "name": "x",
        "stages": [],
        "modes": {},
    }
    text = _yaml.safe_dump(doc)
    with pytest.raises(GeometrySchemaError, match=match_pattern):
        load_geometry_file(text)


# ---------------------------------------------------------------------------
# Strict unknown-key policy
# ---------------------------------------------------------------------------


def test_unknown_top_level_key_rejected():
    bad = _MINIMAL_YAML + "documantation: typo!\n"  # misspelt deliberately
    with pytest.raises(GeometrySchemaError, match="documantation"):
        load_geometry_file(bad)


def test_unknown_stage_key_rejected():
    bad = _MINIMAL_YAML.replace(
        "{name: omega,  axis: -transverse, parent: null,  role: sample}",
        "{name: omega, axis: -transverse, parent: null, role: sample, foo: bar}",
    )
    with pytest.raises(GeometrySchemaError, match=r"foo"):
        load_geometry_file(bad)


def test_unknown_mode_key_rejected():
    bad = _MINIMAL_YAML.replace(
        "        computed: [omega, ttheta]",
        "        computed: [omega, ttheta]\n        flavor: spicy",
    )
    with pytest.raises(GeometrySchemaError, match="flavor"):
        load_geometry_file(bad)


def test_unknown_constraint_key_rejected():
    bad = _MINIMAL_YAML.replace(
        "{type: bisect, stage1: omega, stage2: ttheta}",
        "{type: bisect, stage1: omega, stage2: ttheta, oops: 1}",
    )
    with pytest.raises(GeometrySchemaError, match="oops"):
        load_geometry_file(bad)


# ---------------------------------------------------------------------------
# Basis resolution and precedence
# ---------------------------------------------------------------------------


def test_basis_shorthand_BL():
    g = load_geometry_file(_MINIMAL_YAML)
    np.testing.assert_array_equal(g.basis["vertical"], BASIS_BL["vertical"])
    np.testing.assert_array_equal(g.basis["transverse"], BASIS_BL["transverse"])


def test_basis_shorthand_YOU():
    text = _MINIMAL_YAML.replace("basis: BL", "basis: YOU")
    # YOU has vertical=XHAT, so the omega axis (-transverse) becomes (-ZHAT).
    g = load_geometry_file(text)
    np.testing.assert_array_equal(g.basis["vertical"], BASIS_YOU["vertical"])


def test_basis_omitted_warns_and_uses_default(recwarn):
    """A YAML without ``basis:`` warns and uses BASIS_DEFAULT."""
    text = _MINIMAL_YAML.replace("basis: BL\n", "")
    with pytest.warns(UserWarning, match="BASIS_DEFAULT"):
        g = load_geometry_file(text)
    np.testing.assert_array_equal(g.basis["vertical"], BASIS_DEFAULT["vertical"])


def test_basis_caller_kwarg_overrides_yaml():
    """Caller-supplied basis= wins over the YAML declaration."""
    g = load_geometry_file(_MINIMAL_YAML, basis=BASIS_YOU)
    np.testing.assert_array_equal(g.basis["vertical"], BASIS_YOU["vertical"])


def test_basis_explicit_mapping_strings():
    text = _MINIMAL_YAML.replace(
        "basis: BL",
        "basis: {vertical: '+z', longitudinal: '+y', transverse: '+x'}",
    )
    g = load_geometry_file(text)
    np.testing.assert_array_equal(g.basis["vertical"], np.array([0.0, 0.0, 1.0]))


def test_basis_explicit_mapping_numeric():
    text = _MINIMAL_YAML.replace(
        "basis: BL",
        "basis:\n"
        "    vertical: [0.0, 0.0, 1.0]\n"
        "    longitudinal: [0.0, 1.0, 0.0]\n"
        "    transverse: [1.0, 0.0, 0.0]\n",
    )
    g = load_geometry_file(text)
    np.testing.assert_array_equal(g.basis["vertical"], np.array([0.0, 0.0, 1.0]))


def test_basis_explicit_mapping_numeric_rotated():
    """Non-axis-aligned numeric basis is accepted when orthonormal.

    Verifies the schema and loader both support an arbitrary
    orthonormal triple, e.g. a 45° rotation in the X-Y plane.
    """
    s = 0.7071067811865476  # 1/sqrt(2)
    text = _MINIMAL_YAML.replace(
        "basis: BL",
        "basis:\n"
        f"    vertical: [{s}, {s}, 0.0]\n"
        f"    longitudinal: [-{s}, {s}, 0.0]\n"
        "    transverse: [0.0, 0.0, 1.0]\n",
    )
    g = load_geometry_file(text)
    np.testing.assert_allclose(g.basis["vertical"], [s, s, 0.0])
    np.testing.assert_allclose(g.basis["longitudinal"], [-s, s, 0.0])
    np.testing.assert_allclose(g.basis["transverse"], [0.0, 0.0, 1.0])


@pytest.mark.parametrize(
    "yaml_basis, match_pattern",
    [
        pytest.param("basis: ZZZ", r"unknown basis shorthand", id="bad-shorthand"),
        pytest.param(
            "basis:\n    vertical: [0.5, 0.5, 0.5]\n"
            "    longitudinal: [0.0, 1.0, 0.0]\n"
            "    transverse: [1.0, 0.0, 0.0]\n",
            r"not a unit vector",
            id="non-unit",
        ),
        pytest.param(
            "basis:\n    vertical: [1.0, 0.0, 0.0]\n"
            "    longitudinal: [1.0, 0.0, 0.0]\n"
            "    transverse: [1.0, 0.0, 0.0]\n",
            r"not orthogonal",
            id="non-orthogonal",
        ),
        pytest.param(
            "basis: 42",
            r"shorthand string or a mapping",
            id="basis-not-mapping",
        ),
    ],
)
def test_basis_invalid(yaml_basis, match_pattern):
    text = _MINIMAL_YAML.replace("basis: BL", yaml_basis)
    with pytest.raises(GeometrySchemaError, match=match_pattern):
        load_geometry_file(text)


# ---------------------------------------------------------------------------
# Constraint-type dispatch
# ---------------------------------------------------------------------------


def test_bisect_constraint_built():
    g = load_geometry_file(_MINIMAL_YAML)
    cs = g.modes["bisecting"]
    assert isinstance(cs.constraints[0], BisectConstraint)
    assert cs.constraints[0].sample_stage == "omega"
    assert cs.constraints[0].detector_stage == "ttheta"


def test_sample_and_detector_constraints_built():
    text = _MINIMAL_YAML.replace(
        "    bisecting:\n"
        "        default: true\n"
        "        constraints:\n"
        "            - {type: bisect, stage1: omega, stage2: ttheta}\n"
        "        computed: [omega, ttheta]\n",
        "    bisecting:\n"
        "        default: true\n"
        "        constraints:\n"
        "            - {type: bisect, stage1: omega, stage2: ttheta}\n"
        "        computed: [omega, ttheta]\n"
        "    pinned:\n"
        "        constraints:\n"
        "            - {type: sample, stage: omega, value: 7.5}\n"
        "            - {type: detector, stage: ttheta, value: 15.0}\n"
        "        computed: [omega, ttheta]\n",
    )
    g = load_geometry_file(text)
    pinned = g.modes["pinned"]
    sc, dc = pinned.constraints
    assert isinstance(sc, SampleConstraint) and sc.value == 7.5
    assert isinstance(dc, DetectorConstraint) and dc.value == 15.0


def test_reference_constraint_built():
    text = _MINIMAL_YAML.replace(
        "        computed: [omega, ttheta]\n",
        "        computed: [omega, ttheta]\n"
        "        extras:\n"
        "            n_hat: REQUIRED\n"
        "            psi: null\n",
    )
    text = text.replace(
        "{type: bisect, stage1: omega, stage2: ttheta}",
        "{type: reference, name: psi, value: 0.0}",
    )
    g = load_geometry_file(text)
    cs = g.modes["bisecting"]
    assert isinstance(cs.constraints[0], ReferenceConstraint)
    assert cs.extras["n_hat"] is REQUIRED
    assert cs.extras["psi"] is None


def test_unknown_constraint_type_raises():
    bad = _MINIMAL_YAML.replace(
        "{type: bisect, stage1: omega, stage2: ttheta}",
        "{type: SOMETHING_ELSE, stage1: omega, stage2: ttheta}",
    )
    with pytest.raises(GeometrySchemaError, match="unknown constraint type"):
        load_geometry_file(bad)


def test_two_modes_marked_default_raises():
    bad = _MINIMAL_YAML.replace(
        "        computed: [omega, ttheta]\n",
        "        computed: [omega, ttheta]\n"
        "    second_default:\n"
        "        default: true\n"
        "        constraints:\n"
        "            - {type: sample, stage: omega, value: 0.0}\n"
        "        computed: [ttheta]\n",
    )
    with pytest.raises(GeometrySchemaError, match="more than one mode"):
        load_geometry_file(bad)


# ---------------------------------------------------------------------------
# register_geometry_file()
# ---------------------------------------------------------------------------


@pytest.fixture
def restore_registry():
    """Snapshot and restore _GEOMETRY_REGISTRY around a test."""
    original = dict(_GEOMETRY_REGISTRY)
    try:
        yield
    finally:
        _GEOMETRY_REGISTRY.clear()
        _GEOMETRY_REGISTRY.update(original)


def test_register_geometry_file_uses_filename_stem(tmp_path, restore_registry):
    p = tmp_path / "demo_loader_alpha.yml"
    p.write_text(_MINIMAL_YAML.replace("name: tiny", "name: demo_loader_alpha"))
    name = register_geometry_file(p)
    assert name == "demo_loader_alpha"
    assert "demo_loader_alpha" in _GEOMETRY_REGISTRY


def test_register_geometry_file_explicit_name(tmp_path, restore_registry):
    p = tmp_path / "x.yml"
    p.write_text(_MINIMAL_YAML.replace("name: tiny", "name: ignore_me"))
    name = register_geometry_file(p, name="my_lab_diffr")
    assert name == "my_lab_diffr"
    assert "my_lab_diffr" in _GEOMETRY_REGISTRY


def test_register_geometry_file_duplicate_raises(tmp_path, restore_registry):
    p = tmp_path / "demo_dup.yml"
    p.write_text(_MINIMAL_YAML.replace("name: tiny", "name: demo_dup"))
    register_geometry_file(p)
    with pytest.raises(ValueError, match="already registered"):
        register_geometry_file(p)


def test_register_geometry_file_missing(tmp_path, restore_registry):
    with pytest.raises(FileNotFoundError):
        register_geometry_file(tmp_path / "does_not_exist.yml")


# ---------------------------------------------------------------------------
# load_geometry_file rejects unknown kwargs
# ---------------------------------------------------------------------------


def test_load_geometry_file_rejects_unknown_kwargs(tmp_path):
    p = tmp_path / "x.yml"
    p.write_text(_MINIMAL_YAML)
    with pytest.raises(TypeError, match="unrecognised keyword"):
        load_geometry_file(p, made_up_kwarg=1)


# ---------------------------------------------------------------------------
# Shipped JSON Schema
# ---------------------------------------------------------------------------


def test_get_schema_text_returns_str():
    text = get_schema_text()
    assert isinstance(text, str)
    assert text.lstrip().startswith("{")


def test_get_schema_returns_dict_with_expected_top_level():
    schema = get_schema()
    assert isinstance(schema, dict)
    assert schema["$schema"].startswith("https://json-schema.org/")
    assert KIND_KEY in schema["properties"]
    # required list matches loader behavior
    assert set(schema["required"]) == {
        KIND_KEY,
        "name",
        "stages",
        "modes",
    }


def test_schema_marker_revision_enum_matches_supported_revisions():
    """The schema's revision enum must equal SUPPORTED_REVISIONS."""
    schema = get_schema()
    enum = schema["properties"][KIND_KEY]["properties"]["revision"]["enum"]
    assert set(enum) == set(SUPPORTED_REVISIONS)


def test_schema_constraint_types_match_loader_types():
    """The schema's constraint type enum must equal the loader's accepted set."""
    schema = get_schema()
    types = set()
    for variant in schema["$defs"]["constraint"]["oneOf"]:
        types.add(variant["properties"]["type"]["const"])
    expected = {"bisect", "virtual_bisect", "sample", "detector", "reference"}
    assert types == expected
