# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
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
- :func:`register_geometry_yaml` (in-memory companion: registry
  insertion, required ``name=``, duplicate-name guard, eager schema
  validation, re-parse on every call)
- module-level constants ``KIND_KEY``, ``SUPPORTED_REVISIONS``,
  ``CURRENT_REVISION``
"""

from __future__ import annotations

import os
import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

from ad_hoc_diffractometer import load_geometry_file
from ad_hoc_diffractometer import register_geometry_file
from ad_hoc_diffractometer import register_geometry_yaml
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
    schema_revision: 1
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
            {"schema_revision": "1"},
            r"schema_revision must be an integer",
            id="schema_revision-string",
        ),
        pytest.param(
            {"schema_revision": True},
            r"schema_revision must be an integer",
            id="schema_revision-bool",
        ),
        pytest.param(
            {"schema_revision": 999},
            r"unsupported geometry schema revision 999",
            id="schema_revision-unsupported",
        ),
        pytest.param(
            {},
            r"must contain a 'schema_revision' integer key",
            id="schema_revision-missing",
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
# register_geometry_yaml()  (issue #288 — in-memory companion)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, expected_name, context",
    [
        pytest.param(
            {"name": "inmem_alpha"},
            "inmem_alpha",
            does_not_raise(),
            id="registers-under-supplied-name",
        ),
        pytest.param(
            {},
            None,
            pytest.raises(TypeError, match="name"),
            id="name-is-required",
        ),
    ],
)
def test_register_geometry_yaml_name_argument(
    kwargs, expected_name, context, restore_registry
):
    """``name`` is a required keyword-only argument."""
    with context:
        returned = register_geometry_yaml(_MINIMAL_YAML, **kwargs)
        assert returned == expected_name
        assert expected_name in _GEOMETRY_REGISTRY


def test_register_geometry_yaml_make_geometry_round_trip(restore_registry):
    """make_geometry() on a registered in-memory geometry produces an
    AdHocDiffractometer consistent with load_geometry_file(yaml_text)."""
    import ad_hoc_diffractometer as ahd

    register_geometry_yaml(_MINIMAL_YAML, name="inmem_round_trip")
    g_registry = ahd.make_geometry("inmem_round_trip")
    g_direct = load_geometry_file(_MINIMAL_YAML)
    assert g_registry.name == g_direct.name == "tiny"
    assert list(g_registry._stages) == list(g_direct._stages)
    assert list(g_registry.modes) == list(g_direct.modes)


def test_register_geometry_yaml_duplicate_raises(restore_registry):
    """Second registration under the same name raises ValueError
    (same contract as register_geometry_file)."""
    register_geometry_yaml(_MINIMAL_YAML, name="inmem_dup")
    with pytest.raises(ValueError, match="already registered"):
        register_geometry_yaml(_MINIMAL_YAML, name="inmem_dup")


def test_register_geometry_yaml_eager_schema_validation(restore_registry):
    """Schema errors surface at registration time, not at first
    make_geometry() call (matches register_geometry_file behavior)."""
    # Drop the required top-level 'stages' key.
    bad_yaml = (
        "ad_hoc_diffractometer_geometry:\n"
        "    schema_revision: 1\n"
        "name: broken\n"
        "basis: BL\n"
        "modes:\n"
        "    m:\n"
        "        default: true\n"
        "        constraints: []\n"
        "        computed: []\n"
    )
    with pytest.raises(GeometrySchemaError, match="stages"):
        register_geometry_yaml(bad_yaml, name="inmem_broken")
    # Registration was rejected — nothing landed in the registry.
    assert "inmem_broken" not in _GEOMETRY_REGISTRY


def test_register_geometry_yaml_factory_reparses_text(restore_registry):
    """Each make_geometry() call re-parses the captured YAML text, so
    every call returns a fresh AdHocDiffractometer instance (mirroring
    register_geometry_file's per-call re-read semantics)."""
    import ad_hoc_diffractometer as ahd

    register_geometry_yaml(_MINIMAL_YAML, name="inmem_fresh")
    g1 = ahd.make_geometry("inmem_fresh")
    g2 = ahd.make_geometry("inmem_fresh")
    assert g1 is not g2
    assert g1.name == g2.name == "tiny"


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
    """The schema's schema_revision enum must equal SUPPORTED_REVISIONS."""
    schema = get_schema()
    enum = schema["properties"][KIND_KEY]["properties"]["schema_revision"]["enum"]
    assert set(enum) == set(SUPPORTED_REVISIONS)


def test_schema_constraint_types_match_loader_types():
    """The schema's constraint type enum must equal the loader's accepted set."""
    schema = get_schema()
    types = set()
    for variant in schema["$defs"]["constraint"]["oneOf"]:
        types.add(variant["properties"]["type"]["const"])
    expected = {"bisect", "virtual_bisect", "sample", "detector", "reference"}
    assert types == expected


# ---------------------------------------------------------------------------
# Error-path coverage for `_construct_from_doc` and helpers
# ---------------------------------------------------------------------------
#
# Each parametrised case here exercises one error branch in the loader.
# These tests are intentionally tightly scoped to lift coverage on
# ``geometry_loader.py`` to 100 %; they are not meant to express new
# product behavior beyond what the user-facing tests above already cover.


def _yaml_doc_to_text(doc: dict) -> str:
    """Render a Python dict as a YAML document the loader can parse."""
    import yaml as _yaml

    return _yaml.safe_dump(doc, sort_keys=False)


def _minimal_doc(**overrides):
    """Return a minimal valid geometry document as a Python dict."""
    base = {
        KIND_KEY: {"schema_revision": 1},
        "name": "tiny",
        "documentation": "smoke",
        "basis": "BL",
        "stages": [
            {"name": "omega", "axis": "-transverse", "parent": None, "role": "sample"},
            {
                "name": "ttheta",
                "axis": "-transverse",
                "parent": None,
                "role": "detector",
            },
        ],
        "modes": {
            "bisecting": {
                "default": True,
                "constraints": [
                    {"type": "bisect", "stage1": "omega", "stage2": "ttheta"}
                ],
                "computed": ["omega", "ttheta"],
            }
        },
    }
    base.update(overrides)
    return base


@pytest.mark.parametrize(
    "doc, match_pattern",
    [
        pytest.param(
            _minimal_doc(name=""),
            r"'name' must be a non-empty string",
            id="empty-name",
        ),
        pytest.param(
            _minimal_doc(name=42),
            r"'name' must be a non-empty string",
            id="non-string-name",
        ),
        pytest.param(
            _minimal_doc(documentation=42),
            r"'documentation' must be a string",
            id="non-string-documentation",
        ),
        pytest.param(
            _minimal_doc(parameters="not-a-mapping"),
            r"'parameters' must be a mapping",
            id="parameters-not-mapping",
        ),
        pytest.param(
            _minimal_doc(parameters={"unknown_key": 1}),
            r"unrecognised key.*unknown_key",
            id="parameters-unknown-key",
        ),
    ],
)
def test_top_level_field_errors(doc, match_pattern):
    """Error branches in `_construct_from_doc` for top-level fields."""
    with pytest.raises(GeometrySchemaError, match=match_pattern):
        load_geometry_file(_yaml_doc_to_text(doc))


@pytest.mark.parametrize(
    "doc, match_pattern",
    [
        pytest.param(
            _minimal_doc(stages="oops"),
            r"'stages' must be a non-empty list",
            id="stages-not-list",
        ),
        pytest.param(
            _minimal_doc(stages=[]),
            r"'stages' must be a non-empty list",
            id="stages-empty",
        ),
        pytest.param(
            _minimal_doc(stages=["not-a-mapping"]),
            r"stages\[0\] must be a mapping",
            id="stage-not-mapping",
        ),
        pytest.param(
            _minimal_doc(
                stages=[{"axis": "+x", "parent": None, "role": "sample"}],
            ),
            r"stages\[0\] missing required key 'name'",
            id="stage-missing-name",
        ),
    ],
)
def test_stage_field_errors(doc, match_pattern):
    with pytest.raises(GeometrySchemaError, match=match_pattern):
        load_geometry_file(_yaml_doc_to_text(doc))


@pytest.mark.parametrize(
    "modes_value, match_pattern",
    [
        pytest.param(
            "oops", r"'modes' must be a non-empty mapping", id="modes-not-mapping"
        ),
        pytest.param({}, r"'modes' must be a non-empty mapping", id="modes-empty"),
    ],
)
def test_modes_field_errors(modes_value, match_pattern):
    doc = _minimal_doc()
    doc["modes"] = modes_value
    with pytest.raises(GeometrySchemaError, match=match_pattern):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_mode_name_must_be_string():
    """Non-string mode name raises GeometrySchemaError."""
    doc = _minimal_doc()
    doc["modes"] = {
        42: {
            "constraints": [{"type": "bisect", "stage1": "omega", "stage2": "ttheta"}],
            "computed": ["omega", "ttheta"],
        }
    }
    with pytest.raises(GeometrySchemaError, match="mode names must be strings"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_mode_spec_not_mapping():
    """Mode spec must itself be a mapping."""
    doc = _minimal_doc()
    doc["modes"] = {"weird": "scalar-not-mapping"}
    with pytest.raises(GeometrySchemaError, match="must be a mapping"):
        load_geometry_file(_yaml_doc_to_text(doc))


@pytest.mark.parametrize(
    "constraints_value, match_pattern",
    [
        pytest.param(
            "oops", r"'constraints' must be a list", id="constraints-not-list"
        ),
    ],
)
def test_mode_constraints_field_errors(constraints_value, match_pattern):
    doc = _minimal_doc()
    doc["modes"] = {
        "wonky": {"constraints": constraints_value, "computed": []},
    }
    with pytest.raises(GeometrySchemaError, match=match_pattern):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_mode_computed_must_be_list():
    doc = _minimal_doc()
    doc["modes"] = {
        "wonky": {"constraints": [], "computed": "not-a-list"},
    }
    with pytest.raises(GeometrySchemaError, match="'computed' must be a list"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_mode_extras_must_be_mapping():
    doc = _minimal_doc()
    doc["modes"] = {
        "wonky": {"constraints": [], "computed": [], "extras": "not-a-mapping"},
    }
    with pytest.raises(GeometrySchemaError, match="'extras' must be a mapping"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_mode_cut_points_must_be_mapping():
    doc = _minimal_doc()
    doc["modes"] = {
        "wonky": {"constraints": [], "computed": [], "cut_points": "not-a-mapping"},
    }
    with pytest.raises(GeometrySchemaError, match="'cut_points' must be a mapping"):
        load_geometry_file(_yaml_doc_to_text(doc))


@pytest.mark.parametrize(
    "constraint_spec, match_pattern",
    [
        pytest.param(
            "scalar",
            r"each constraint must be a mapping",
            id="constraint-not-mapping",
        ),
        pytest.param(
            {"stage1": "omega", "stage2": "ttheta"},
            r"missing required 'type' key",
            id="constraint-missing-type",
        ),
        pytest.param(
            {"type": "bisect", "stage1": "omega"},
            r"bisect requires 'stage2'",
            id="bisect-missing-stage2",
        ),
        pytest.param(
            {"type": "virtual_bisect", "stage1": "omega"},
            r"virtual_bisect requires 'stage2'",
            id="virtual_bisect-missing-stage2",
        ),
        pytest.param(
            {"type": "sample", "stage": "omega"},
            r"sample requires 'value'",
            id="sample-missing-value",
        ),
        pytest.param(
            {"type": "detector", "value": 0.0},
            r"detector requires 'stage'",
            id="detector-missing-stage",
        ),
        pytest.param(
            {"type": "reference", "value": 0.0},
            r"reference requires 'name'",
            id="reference-missing-name",
        ),
    ],
)
def test_constraint_dispatch_errors(constraint_spec, match_pattern):
    doc = _minimal_doc()
    doc["modes"] = {
        "wonky": {
            "constraints": [constraint_spec],
            "computed": [],
        },
    }
    with pytest.raises(GeometrySchemaError, match=match_pattern):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_basis_omitted_uses_default_warning_value():
    """When `basis:` is omitted, the loader emits a UserWarning AND uses
    BASIS_DEFAULT (covers the warnings.warn branch)."""
    import numpy as np

    doc = _minimal_doc()
    del doc["basis"]
    with pytest.warns(UserWarning, match="BASIS_DEFAULT"):
        g = load_geometry_file(_yaml_doc_to_text(doc))
    np.testing.assert_array_equal(g.basis["vertical"], BASIS_DEFAULT["vertical"])


def test_basis_caller_dict_missing_keys_raises():
    """Caller-supplied dict missing required basis keys raises."""
    doc = _minimal_doc()
    with pytest.raises(GeometrySchemaError, match="caller-supplied basis is missing"):
        load_geometry_file(
            _yaml_doc_to_text(doc),
            basis={"vertical": [1.0, 0.0, 0.0]},  # missing two keys
        )


def test_basis_caller_string_shorthand():
    """Caller may pass a shorthand string instead of a dict."""
    import numpy as np

    doc = _minimal_doc()
    g = load_geometry_file(_yaml_doc_to_text(doc), basis="YOU")
    # YOU has vertical=XHAT
    np.testing.assert_array_equal(g.basis["vertical"], [1.0, 0.0, 0.0])


def test_basis_numeric_vector_non_numeric_entries_rejected():
    """Numeric basis vector with non-numeric entries raises."""
    doc = _minimal_doc()
    doc["basis"] = {
        "vertical": ["a", "b", "c"],
        "longitudinal": [0.0, 1.0, 0.0],
        "transverse": [1.0, 0.0, 0.0],
    }
    with pytest.raises(GeometrySchemaError, match="numeric vector entries"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_basis_value_invalid_type_in_mapping():
    """Basis-mapping value of a wrong type raises."""
    doc = _minimal_doc()
    doc["basis"] = {
        "vertical": 42,  # not a string, not a vector
        "longitudinal": [0.0, 1.0, 0.0],
        "transverse": [1.0, 0.0, 0.0],
    }
    with pytest.raises(GeometrySchemaError, match="signed-axis string"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_basis_string_lookup_failure_in_mapping():
    """A bad signed-axis string inside a basis mapping raises."""
    doc = _minimal_doc()
    doc["basis"] = {
        "vertical": "+nonsense",
        "longitudinal": "+y",
        "transverse": "+x",
    }
    with pytest.raises(GeometrySchemaError, match=r"basis\['vertical'\]"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_kappa_chi_eq_invalid_form():
    """kappa_chi_eq must be a string or length-3 numeric vector."""
    doc = _minimal_doc()
    doc["parameters"] = {"alpha_deg": 50.0}
    doc["kappa_chi_eq"] = 42  # neither
    with pytest.raises(GeometrySchemaError, match="'kappa_chi_eq' must be"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_kappa_chi_eq_numeric_vector_form():
    """kappa_chi_eq may be a length-3 numeric vector.

    With issue #284 the ``kappa_chi_eq`` field controls only the
    kappa-arm tilt direction (input to Walko's formula in
    ``_resolve_axis``).  This test verifies the numeric-vector form
    parses without error and that the kappa stage axis is built
    correctly from it; the synthesized convention's ``n_chi_eq`` is
    a separate concept (auto-derived from the basis, here
    ``+longitudinal``) and is verified by other tests.
    """
    import math

    import numpy as np

    doc = {
        KIND_KEY: {"schema_revision": 1},
        "name": "kchi_numeric",
        "documentation": "test",
        "basis": "BL",
        "parameters": {"alpha_deg": 50.0},
        "kappa_chi_eq": [0.0, 0.0, 1.0],
        "stages": [
            {
                "name": "komega",
                "axis": "-transverse",
                "parent": None,
                "role": "sample",
            },
            {
                "name": "kappa",
                "axis": {"kappa_eulerian": "+transverse"},
                "parent": "komega",
                "role": "sample",
            },
            {
                "name": "kphi",
                "axis": "-transverse",
                "parent": "kappa",
                "role": "sample",
            },
            {
                "name": "ttheta",
                "axis": "-transverse",
                "parent": None,
                "role": "detector",
            },
        ],
        "modes": {
            "bisecting": {
                "default": True,
                "constraints": [
                    {"type": "virtual_bisect", "stage1": "omega", "stage2": "ttheta"}
                ],
                "computed": ["komega", "kappa", "kphi", "ttheta"],
            }
        },
    }
    g = load_geometry_file(_yaml_doc_to_text(doc))
    # The kappa stage axis is built from the unsigned outer (+transverse)
    # tilted alpha_deg toward kappa_chi_eq (+vertical):
    expected_kappa = np.cos(math.radians(50.0)) * np.array([1.0, 0.0, 0.0]) + np.sin(
        math.radians(50.0)
    ) * np.array([0.0, 0.0, 1.0])
    np.testing.assert_allclose(g.stage("kappa").axis, expected_kappa, atol=1e-12)
    # The synthesized convention's n_chi_eq is auto-derived as the
    # first basis direction perpendicular to n_komega; for BL with
    # n_komega = -transverse that yields +longitudinal.
    np.testing.assert_allclose(
        g.kappa_pseudo_angle_convention.n_chi_eq, [0.0, 1.0, 0.0], atol=1e-12
    )


def test_kappa_default_alpha_when_kappa_chi_eq_present():
    """alpha_deg defaults to KAPPA_ALPHA_DEFAULT when only kappa_chi_eq is set."""

    doc = {
        KIND_KEY: {"schema_revision": 1},
        "name": "kappa_default_alpha",
        "documentation": "test",
        "basis": "BL",
        "kappa_chi_eq": "+vertical",
        "stages": [
            {
                "name": "komega",
                "axis": "-transverse",
                "parent": None,
                "role": "sample",
            },
            {
                "name": "kappa",
                "axis": {"kappa_eulerian": "+transverse"},
                "parent": "komega",
                "role": "sample",
            },
            {
                "name": "kphi",
                "axis": "-transverse",
                "parent": "kappa",
                "role": "sample",
            },
            {
                "name": "ttheta",
                "axis": "-transverse",
                "parent": None,
                "role": "detector",
            },
        ],
        "modes": {
            "bisecting": {
                "default": True,
                "constraints": [
                    {"type": "virtual_bisect", "stage1": "omega", "stage2": "ttheta"}
                ],
                "computed": ["komega", "kappa", "kphi", "ttheta"],
            }
        },
    }
    g = load_geometry_file(_yaml_doc_to_text(doc))
    # The kappa stage axis was synthesized via kappa_axis_from_eulerian
    # using the numeric komega vector; sanity-check it's a unit vector.
    import numpy as np

    np.testing.assert_allclose(np.linalg.norm(g._stages["kappa"].axis), 1.0, atol=1e-12)  # noqa: SLF001


def test_axis_invalid_form_rejected():
    """Stage axis that is neither string, list, nor kappa_eulerian raises."""
    doc = _minimal_doc()
    doc["stages"][0]["axis"] = 42
    with pytest.raises(GeometrySchemaError, match="must be a signed-axis string"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_axis_bad_signed_string_rejected():
    """A bad signed-axis string in a stage axis surfaces as a stage error."""
    doc = _minimal_doc()
    doc["stages"][0]["axis"] = "+nonsense"
    with pytest.raises(GeometrySchemaError, match=r"stage 'omega'"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_axis_numeric_vector_non_numeric_entries_rejected():
    """Numeric stage axis with non-numeric entries raises."""
    doc = _minimal_doc()
    doc["stages"][0]["axis"] = ["a", "b", "c"]
    with pytest.raises(GeometrySchemaError, match="numeric axis entries"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_axis_numeric_vector_form_accepted():
    """Numeric stage axis is accepted (covers happy path)."""
    import numpy as np

    doc = _minimal_doc()
    doc["stages"][0]["axis"] = [-1.0, 0.0, 0.0]
    g = load_geometry_file(_yaml_doc_to_text(doc))
    np.testing.assert_array_equal(g._stages["omega"].axis, [-1.0, 0.0, 0.0])  # noqa: SLF001


def test_axis_kappa_eulerian_without_alpha_or_chi_eq_rejected():
    """kappa_eulerian axis form requires both alpha_deg and kappa_chi_eq."""
    doc = _minimal_doc()
    doc["stages"][0]["axis"] = {"kappa_eulerian": "+transverse"}
    with pytest.raises(GeometrySchemaError, match=r"requires the top-level"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_axis_kappa_eulerian_extra_keys_rejected():
    """kappa_eulerian axis spec rejects unknown sub-keys."""
    doc = _minimal_doc()
    doc["parameters"] = {"alpha_deg": 50.0}
    doc["kappa_chi_eq"] = "+vertical"
    doc["stages"][0]["axis"] = {"kappa_eulerian": "+transverse", "extra": 1}
    with pytest.raises(GeometrySchemaError, match=r"only one key"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_axis_kappa_eulerian_invalid_inner_value_rejected():
    """kappa_eulerian inner value must be string or length-3 numeric vector."""
    doc = _minimal_doc()
    doc["parameters"] = {"alpha_deg": 50.0}
    doc["kappa_chi_eq"] = "+vertical"
    doc["stages"][0]["axis"] = {"kappa_eulerian": 42}
    with pytest.raises(GeometrySchemaError, match=r"'kappa_eulerian' value must be"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_axis_kappa_eulerian_inner_numeric_vector_form():
    """kappa_eulerian inner value may be a length-3 numeric vector."""
    import numpy as np

    doc = {
        KIND_KEY: {"schema_revision": 1},
        "name": "kev_numeric",
        "documentation": "test",
        "basis": "BL",
        "parameters": {"alpha_deg": 50.0},
        "kappa_chi_eq": "+vertical",
        "stages": [
            {
                "name": "komega",
                "axis": "-transverse",
                "parent": None,
                "role": "sample",
            },
            {
                "name": "kappa",
                "axis": {"kappa_eulerian": [1.0, 0.0, 0.0]},
                "parent": "komega",
                "role": "sample",
            },
            {
                "name": "kphi",
                "axis": "-transverse",
                "parent": "kappa",
                "role": "sample",
            },
            {
                "name": "ttheta",
                "axis": "-transverse",
                "parent": None,
                "role": "detector",
            },
        ],
        "modes": {
            "bisecting": {
                "default": True,
                "constraints": [
                    {"type": "virtual_bisect", "stage1": "omega", "stage2": "ttheta"}
                ],
                "computed": ["komega", "kappa", "kphi", "ttheta"],
            }
        },
    }
    g = load_geometry_file(_yaml_doc_to_text(doc))
    np.testing.assert_allclose(np.linalg.norm(g._stages["kappa"].axis), 1.0, atol=1e-12)  # noqa: SLF001


def test_kappa_synthesis_missing_canonical_stage_names_raises():
    """If parameters.alpha_deg + kappa_chi_eq are set but the kappa stage
    names are not canonical (komega/kappa/kphi), the loader raises."""
    doc = {
        KIND_KEY: {"schema_revision": 1},
        "name": "weird_kappa",
        "documentation": "test",
        "basis": "BL",
        "parameters": {"alpha_deg": 50.0},
        "kappa_chi_eq": "+vertical",
        "stages": [
            # Note: not named komega/kappa/kphi — names are wholly different.
            {"name": "alpha", "axis": "+vertical", "parent": None, "role": "sample"},
        ],
        "modes": {
            "any": {
                "default": True,
                "constraints": [],
                "computed": ["alpha"],
            }
        },
    }
    with pytest.raises(GeometrySchemaError, match="canonical kappa stage names"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_load_yaml_text_with_yaml_parser_error_falls_back_to_path():
    """A string that fails YAML parsing is treated as a path."""
    # Tab + colon in a way YAML rejects as scanner error,
    # AND not a valid file path.
    bad = "\t: not yaml\n@@@"
    with pytest.raises(FileNotFoundError):
        load_geometry_file(bad)


def test_load_geometry_file_path_object_missing_raises():
    """An explicit Path that doesn't exist raises FileNotFoundError."""
    from pathlib import Path

    with pytest.raises(FileNotFoundError, match="not found"):
        load_geometry_file(Path("/no/such/path.yml"))


def test_register_geometry_decorator_round_trip(restore_registry):
    """The @register_geometry decorator (factories.py) registers a callable."""
    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer import register_geometry
    from ad_hoc_diffractometer.factories import _GEOMETRY_REGISTRY  # noqa: SLF001

    @register_geometry
    def _my_demo_geom_267():
        # Build a trivial geometry directly via the loader.
        return ahd.load_geometry_file(_MINIMAL_YAML)

    try:
        assert "_my_demo_geom_267" in _GEOMETRY_REGISTRY
        g = ahd.make_geometry("_my_demo_geom_267")
        assert g.name == "tiny"
    finally:
        _GEOMETRY_REGISTRY.pop("_my_demo_geom_267", None)


def test_construct_from_doc_missing_marker_directly():
    from ad_hoc_diffractometer.geometry_loader import _construct_from_doc

    with pytest.raises(GeometrySchemaError, match="missing top-level marker"):
        _construct_from_doc(
            {"name": "x", "stages": [], "modes": {}},
            source_label="t",
            overrides={},
        )


def test_construct_from_doc_top_level_not_mapping():
    from ad_hoc_diffractometer.geometry_loader import _construct_from_doc

    with pytest.raises(GeometrySchemaError, match="must be a mapping"):
        _construct_from_doc(["not", "a", "mapping"], source_label="t", overrides={})


def test_basis_mapping_extra_keys_rejected():
    doc = _minimal_doc()
    doc["basis"] = {
        "vertical": "+z",
        "longitudinal": "+y",
        "transverse": "+x",
        "extra_key": "+x",
    }
    with pytest.raises(GeometrySchemaError, match="unrecognised key"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_basis_mapping_missing_keys_rejected():
    doc = _minimal_doc()
    doc["basis"] = {"vertical": "+z"}
    with pytest.raises(GeometrySchemaError, match="missing required key"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_basis_value_not_a_mapping_or_string_rejected():
    doc = _minimal_doc()
    doc["basis"] = 42
    with pytest.raises(GeometrySchemaError, match="shorthand string or a mapping"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_axis_resolve_final_fallback_unknown_form():
    """Stage axis spec that is none of the accepted forms hits the
    final fallback raise in `_resolve_axis`."""
    doc = _minimal_doc()
    # An empty dict (no kappa_eulerian key) lands at the last raise in _resolve_axis.
    doc["stages"][0]["axis"] = {"unrelated_key": 1}
    with pytest.raises(GeometrySchemaError, match="must be a signed-axis string"):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_register_geometry_file_factory_rereads_file(tmp_path, restore_registry):
    """The registered factory re-reads the file on each call so on-disk
    edits are picked up without re-registration (covers the re-read path)."""
    p = tmp_path / "rereadable.yml"
    p.write_text(_MINIMAL_YAML.replace("name: tiny", "name: rereadable"))
    name = register_geometry_file(p)
    g1 = load_geometry_file(p)  # not via registry
    import ad_hoc_diffractometer as ahd

    g2 = ahd.make_geometry(name)  # via registry — exercises the re-read path
    assert g1.name == g2.name == "rereadable"


def test_top_level_missing_name():
    doc = _minimal_doc()
    del doc["name"]
    with pytest.raises(
        GeometrySchemaError, match="missing required top-level key 'name'"
    ):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_top_level_missing_stages():
    doc = _minimal_doc()
    del doc["stages"]
    with pytest.raises(
        GeometrySchemaError, match="missing required top-level key 'stages'"
    ):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_top_level_missing_modes():
    doc = _minimal_doc()
    del doc["modes"]
    with pytest.raises(
        GeometrySchemaError, match="missing required top-level key 'modes'"
    ):
        load_geometry_file(_yaml_doc_to_text(doc))


def test_basis_defensive_det_check_via_monkeypatch(monkeypatch):
    """Defensive det-check in `_validate_basis` when somehow det is not 1."""
    import numpy as np

    from ad_hoc_diffractometer.geometry_loader import _validate_basis

    # All three vectors are unit-length and pairwise orthogonal in finite
    # precision, but we monkey-patch np.linalg.det to force a non-unit
    # determinant so the defensive branch fires.
    real_det = np.linalg.det
    monkeypatch.setattr(np.linalg, "det", lambda *a, **kw: 0.0)
    basis = {
        "vertical": np.array([0.0, 0.0, 1.0]),
        "longitudinal": np.array([0.0, 1.0, 0.0]),
        "transverse": np.array([1.0, 0.0, 0.0]),
    }
    with pytest.raises(GeometrySchemaError, match="degenerate or non-orthonormal"):
        _validate_basis(basis)
    monkeypatch.setattr(np.linalg, "det", real_det)
