# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.reflection.

Covers:
  - Reflection dataclass: construction, normalisation, validation, __eq__, __repr__
  - ReflectionList: add, remove, clear, dict-like interface,
    setor1/setor2, orienting_reflections, cross-geometry safety
  - AdHocDiffractometer.add_reflection() convenience wrapper
  - Reflection.to_dict() / from_dict(): hkl, angles, wavelength, name
  - ReflectionList.to_dict() / from_dict(): ordering, or1/or2 preserved
"""

import json
import re
from contextlib import nullcontext as does_not_raise

import pytest

from ad_hoc_diffractometer import Reflection
from ad_hoc_diffractometer import ReflectionList
from ad_hoc_diffractometer import kappa6c
from ad_hoc_diffractometer import psic

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_PSIC_STAGES = {"mu", "eta", "chi", "phi", "nu", "delta"}
_PSIC_ANGLES = {
    "mu": 0.0,
    "eta": 20.0,
    "chi": 90.0,
    "phi": 0.0,
    "nu": 0.0,
    "delta": 40.0,
}


@pytest.fixture
def rl():
    """Fresh ReflectionList for psic geometry."""
    return ReflectionList(geometry_name="psic", valid_stages=_PSIC_STAGES)


# ---------------------------------------------------------------------------
# Reflection dataclass
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "hkl, angles, wavelength, context",
    [
        pytest.param(
            (1, 0, 0),
            {"mu": 0.0},
            1.5406,
            does_not_raise(),
            id="valid-integer-hkl",
        ),
        pytest.param(
            (0.5, 0.0, 0.0),
            {"mu": 0.0},
            1.5406,
            does_not_raise(),
            id="valid-fractional-hkl",
        ),
        pytest.param(
            (1, 0, 0),
            {"mu": 0.0},
            None,
            does_not_raise(),
            id="valid-wavelength-none",
        ),
        pytest.param(
            (1, 0, 0),
            {"mu": 0.0},
            0.0,
            pytest.raises(ValueError, match=re.escape("must be > 0")),
            id="invalid-wavelength-zero",
        ),
        pytest.param(
            (1, 0, 0),
            {"mu": 0.0},
            -1.5,
            pytest.raises(ValueError, match=re.escape("must be > 0")),
            id="invalid-wavelength-negative",
        ),
    ],
)
def test_reflection_construction(hkl, angles, wavelength, context):
    with context:
        r = Reflection(name="r1", hkl=hkl, angles=angles, wavelength=wavelength)
        assert r.hkl == tuple(float(v) for v in hkl)
        assert all(isinstance(v, float) for v in r.angles.values())


def test_reflection_name_stored():
    r = Reflection(name="Si_111", hkl=(1, 1, 1), angles={})
    assert r.name == "Si_111"


def test_reflection_hkl_normalised_to_float():
    r = Reflection(name="r1", hkl=(1, 2, 3), angles={})
    assert r.hkl == (1.0, 2.0, 3.0)


def test_reflection_geometry_name_stored():
    r = Reflection(name="r1", hkl=(1, 0, 0), angles={}, geometry_name="psic")
    assert r.geometry_name == "psic"


def test_reflection_geometry_name_default_none():
    r = Reflection(name="r1", hkl=(1, 0, 0), angles={})
    assert r.geometry_name is None


def test_reflection_repr_contains_name_and_geometry():
    r = Reflection(name="Si_111", hkl=(1, 1, 1), angles={}, geometry_name="psic")
    assert "Si_111" in repr(r)
    assert "psic" in repr(r)


# ---------------------------------------------------------------------------
# Reflection.__eq__
# ---------------------------------------------------------------------------


def test_reflection_eq_identical():
    r1 = Reflection(
        name="r1",
        hkl=(1, 0, 0),
        angles={"mu": 0.0},
        wavelength=1.5406,
        geometry_name="psic",
    )
    r2 = Reflection(
        name="r1",
        hkl=(1, 0, 0),
        angles={"mu": 0.0},
        wavelength=1.5406,
        geometry_name="psic",
    )
    assert r1 == r2


def test_reflection_eq_different_name():
    r1 = Reflection(name="r1", hkl=(1, 0, 0), angles={}, geometry_name="psic")
    r2 = Reflection(name="r2", hkl=(1, 0, 0), angles={}, geometry_name="psic")
    assert r1 != r2


def test_reflection_eq_different_hkl():
    r1 = Reflection(name="r1", hkl=(1, 0, 0), angles={}, geometry_name="psic")
    r2 = Reflection(name="r1", hkl=(0, 1, 0), angles={}, geometry_name="psic")
    assert r1 != r2


def test_reflection_eq_different_angles():
    r1 = Reflection(name="r1", hkl=(1, 0, 0), angles={"mu": 0.0}, geometry_name="psic")
    r2 = Reflection(name="r1", hkl=(1, 0, 0), angles={"mu": 5.0}, geometry_name="psic")
    assert r1 != r2


def test_reflection_eq_different_wavelength():
    r1 = Reflection(
        name="r1", hkl=(1, 0, 0), angles={}, wavelength=1.5406, geometry_name="psic"
    )
    r2 = Reflection(
        name="r1", hkl=(1, 0, 0), angles={}, wavelength=0.7107, geometry_name="psic"
    )
    assert r1 != r2


def test_reflection_eq_different_geometry_name():
    r1 = Reflection(name="r1", hkl=(1, 0, 0), angles={}, geometry_name="psic")
    r2 = Reflection(name="r1", hkl=(1, 0, 0), angles={}, geometry_name="kappa6c")
    assert r1 != r2


def test_reflection_eq_not_implemented_for_non_reflection():
    r = Reflection(name="r1", hkl=(1, 0, 0), angles={})
    assert r.__eq__("not a reflection") is NotImplemented
    assert r != "not a reflection"


def test_reflection_eq_within_default_tolerance():
    """hkl values differing by less than default atol compare equal."""
    r1 = Reflection(
        name="r1", hkl=(1.0, 0.0, 0.0), angles={"mu": 20.0}, geometry_name="psic"
    )
    r2 = Reflection(
        name="r1", hkl=(1.0000003, 0.0, 0.0), angles={"mu": 20.0}, geometry_name="psic"
    )
    assert r1 == r2  # within default atol=5e-7


def test_reflection_eq_outside_default_tolerance():
    r1 = Reflection(
        name="r1", hkl=(1.0, 0.0, 0.0), angles={"mu": 20.0}, geometry_name="psic"
    )
    r2 = Reflection(
        name="r1", hkl=(1.001, 0.0, 0.0), angles={"mu": 20.0}, geometry_name="psic"
    )
    assert r1 != r2


def test_reflection_eq_explicit_atol():
    r1 = Reflection(
        name="r1", hkl=(1.0, 0.0, 0.0), angles={"mu": 20.0}, geometry_name="psic"
    )
    r2 = Reflection(
        name="r1", hkl=(1.005, 0.0, 0.0), angles={"mu": 20.0}, geometry_name="psic"
    )
    assert r1.__eq__(r2, atol=0.01) is True
    assert r1.__eq__(r2, atol=0.001) is False


def test_reflection_eq_angles_tolerance():
    """Angle values also compared with tolerance."""
    r1 = Reflection(name="r1", hkl=(1, 0, 0), angles={"mu": 20.0}, geometry_name="psic")
    r2 = Reflection(
        name="r1", hkl=(1, 0, 0), angles={"mu": 20.0000003}, geometry_name="psic"
    )
    assert r1 == r2  # within default atol


def test_reflection_eq_wavelength_tolerance():
    r1 = Reflection(
        name="r1", hkl=(1, 0, 0), angles={}, wavelength=1.5406000, geometry_name="psic"
    )
    r2 = Reflection(
        name="r1", hkl=(1, 0, 0), angles={}, wavelength=1.5406003, geometry_name="psic"
    )
    assert r1 == r2  # within default atol


# ---------------------------------------------------------------------------
# ReflectionList construction
# ---------------------------------------------------------------------------


def test_reflection_list_empty_by_default(rl):
    assert len(rl) == 0
    assert list(rl) == []


def test_reflection_list_geometry_name(rl):
    assert rl.geometry_name == "psic"


def test_reflection_list_valid_stages(rl):
    assert rl.valid_stages == _PSIC_STAGES


def test_reflection_list_repr(rl):
    assert "psic" in repr(rl)


# ---------------------------------------------------------------------------
# ReflectionList.add()
# ---------------------------------------------------------------------------


def test_add_returns_reflection(rl):
    r = rl.add("r1", hkl=(1, 0, 0), angles={"mu": 0.0})
    assert isinstance(r, Reflection)
    assert r.name == "r1"
    assert r.hkl == (1.0, 0.0, 0.0)
    assert r.geometry_name == "psic"


def test_add_wavelength_stored(rl):
    r = rl.add("r1", hkl=(1, 0, 0), angles={}, wavelength=1.5406)
    assert r.wavelength == pytest.approx(1.5406)


def test_add_duplicate_name_raises(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    with pytest.raises(ValueError, match=re.escape("already exists")):
        rl.add("r1", hkl=(0, 1, 0), angles={})


def test_add_unknown_angle_key_raises(rl):
    with pytest.raises(ValueError, match=re.escape("not stage names")):
        rl.add("r1", hkl=(1, 0, 0), angles={"bogus": 0.0})


def test_add_subset_of_stages_allowed(rl):
    r = rl.add("r1", hkl=(1, 0, 0), angles={"mu": 0.0, "eta": 20.0})
    assert set(r.angles) == {"mu", "eta"}


# ---------------------------------------------------------------------------
# ReflectionList dict-like interface
# ---------------------------------------------------------------------------


def test_getitem(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    assert rl["r1"].name == "r1"


def test_getitem_missing_raises(rl):
    with pytest.raises(KeyError):
        _ = rl["missing"]


def test_contains(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    assert "r1" in rl
    assert "missing" not in rl


def test_len(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.add("r2", hkl=(0, 1, 0), angles={})
    assert len(rl) == 2


def test_iter_ordering(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.add("r2", hkl=(0, 1, 0), angles={})
    rl.add("r3", hkl=(0, 0, 1), angles={})
    assert list(rl) == ["r1", "r2", "r3"]


def test_delitem_clears_or_designation(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.setor1("r1")
    del rl["r1"]
    assert "r1" not in rl
    assert rl.orienting_reflections == []


# ---------------------------------------------------------------------------
# ReflectionList.remove() and .clear()
# ---------------------------------------------------------------------------


def test_remove(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.remove("r1")
    assert "r1" not in rl


def test_remove_missing_raises(rl):
    with pytest.raises(KeyError):
        rl.remove("missing")


def test_remove_then_readd(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.remove("r1")
    rl.add("r1", hkl=(0, 1, 0), angles={})
    assert rl["r1"].hkl == (0.0, 1.0, 0.0)


def test_clear(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.add("r2", hkl=(0, 1, 0), angles={})
    rl.setor1("r1")
    rl.setor2("r2")
    rl.clear()
    assert len(rl) == 0
    assert rl.orienting_reflections == []


# ---------------------------------------------------------------------------
# setor1 / setor2 / orienting_reflections
# ---------------------------------------------------------------------------


def test_orienting_reflections_empty_by_default(rl):
    assert rl.orienting_reflections == []


def test_setor1_by_name(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.setor1("r1")
    ors = rl.orienting_reflections
    assert len(ors) == 1
    assert ors[0].name == "r1"


def test_setor1_by_object(rl):
    r = rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.setor1(r)
    assert rl.orienting_reflections[0].name == "r1"


def test_setor2_by_name(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.add("r2", hkl=(0, 1, 0), angles={})
    rl.setor1("r1")
    rl.setor2("r2")
    ors = rl.orienting_reflections
    assert len(ors) == 2
    assert ors[0].name == "r1"
    assert ors[1].name == "r2"


def test_setor2_by_object(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    r2 = rl.add("r2", hkl=(0, 1, 0), angles={})
    rl.setor1("r1")
    rl.setor2(r2)
    assert rl.orienting_reflections[1].name == "r2"


def test_only_or1_set(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.setor1("r1")
    assert len(rl.orienting_reflections) == 1


def test_setor1_replaces_previous(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.add("r2", hkl=(0, 1, 0), angles={})
    rl.setor1("r1")
    rl.setor1("r2")
    ors = rl.orienting_reflections
    assert len(ors) == 1
    assert ors[0].name == "r2"
    assert "r1" in rl


def test_setor2_replaces_previous(rl):
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.add("r2", hkl=(0, 1, 0), angles={})
    rl.add("r3", hkl=(0, 0, 1), angles={})
    rl.setor1("r1")
    rl.setor2("r2")
    rl.setor2("r3")
    ors = rl.orienting_reflections
    assert len(ors) == 2
    assert ors[1].name == "r3"
    assert "r2" in rl


def test_setor2_was_or1_clears_or1(rl):
    """Moving a reflection from or1 to or2 clears the primary slot."""
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.setor1("r1")
    rl.setor2("r1")  # r1 moves to or2; or1 becomes None
    ors = rl.orienting_reflections
    assert len(ors) == 1
    assert ors[0].name == "r1"


def test_setor1_was_or2_clears_or2(rl):
    """Moving a reflection from or2 to or1 clears the secondary slot."""
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.add("r2", hkl=(0, 1, 0), angles={})
    rl.setor1("r1")
    rl.setor2("r2")
    rl.setor1("r2")  # r2 moves from or2 to or1
    ors = rl.orienting_reflections
    assert len(ors) == 1
    assert ors[0].name == "r2"


def test_setor1_unknown_raises(rl):
    with pytest.raises(KeyError):
        rl.setor1("missing")


def test_setor2_unknown_raises(rl):
    with pytest.raises(KeyError):
        rl.setor2("missing")


def test_delete_or2_clears_secondary_slot(rl):
    """Deleting the reflection designated as or2 clears the secondary slot."""
    rl.add("r1", hkl=(1, 0, 0), angles={})
    rl.add("r2", hkl=(0, 1, 0), angles={})
    rl.setor1("r1")
    rl.setor2("r2")
    del rl["r2"]
    assert "r2" not in rl
    ors = rl.orienting_reflections
    # or1 remains; or2 slot is now None
    assert len(ors) == 1
    assert ors[0].name == "r1"


# ---------------------------------------------------------------------------
# Cross-geometry safety
# ---------------------------------------------------------------------------


def test_cross_geometry_angle_keys_rejected():
    """psic stage names are not valid in kappa6c."""
    g = kappa6c()
    with pytest.raises(ValueError, match=re.escape("not stage names")):
        g.reflections.add(
            "r1", hkl=(1, 0, 0), angles={"eta": 20.0, "chi": 90.0, "phi": 0.0}
        )


def test_cross_geometry_geometry_name_differs():
    g1 = psic()
    g2 = kappa6c()
    r = g1.add_reflection("r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES)
    assert r.geometry_name == "psic"
    assert r.geometry_name != g2.name


# ---------------------------------------------------------------------------
# AdHocDiffractometer.add_reflection() convenience wrapper
# ---------------------------------------------------------------------------


def test_add_reflection_inherits_wavelength(psic_geom):
    psic_geom.wavelength = 1.5406
    r = psic_geom.add_reflection("r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES)
    assert r.wavelength == pytest.approx(1.5406)


def test_add_reflection_explicit_wavelength_overrides(psic_geom):
    psic_geom.wavelength = 1.5406
    r = psic_geom.add_reflection(
        "r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES, wavelength=0.7107
    )
    assert r.wavelength == pytest.approx(0.7107)


def test_add_reflection_stored_in_active_sample(psic_geom):
    psic_geom.add_reflection("r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES)
    assert "r1" in psic_geom.sample.reflections


def test_add_reflection_geometry_name_set(psic_geom):
    r = psic_geom.add_reflection("r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES)
    assert r.geometry_name == "psic"


# ---------------------------------------------------------------------------
# Reflection.__eq__ branch coverage
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "r1_kwargs, r2_kwargs, expected_equal, context",
    [
        pytest.param(
            {"name": "a", "hkl": (1, 0, 0), "angles": {"omega": 10.0}},
            {"name": "b", "hkl": (1, 0, 0), "angles": {"omega": 10.0}},
            False,
            does_not_raise(),
            id="different-name",
        ),
        pytest.param(
            {
                "name": "r",
                "hkl": (1, 0, 0),
                "angles": {"omega": 10.0},
                "geometry_name": "psic",
            },
            {
                "name": "r",
                "hkl": (1, 0, 0),
                "angles": {"omega": 10.0},
                "geometry_name": "fourcv",
            },
            False,
            does_not_raise(),
            id="different-geometry-name",
        ),
        pytest.param(
            {"name": "r", "hkl": (1, 0, 0), "angles": {"omega": 10.0}},
            {"name": "r", "hkl": (1, 0, 0), "angles": {"omega": 10.0, "chi": 0.0}},
            False,
            does_not_raise(),
            id="different-angle-keys",
        ),
        pytest.param(
            {
                "name": "r",
                "hkl": (1, 0, 0),
                "angles": {"omega": 10.0},
                "wavelength": 1.5,
            },
            {
                "name": "r",
                "hkl": (1, 0, 0),
                "angles": {"omega": 10.0},
                "wavelength": None,
            },
            False,
            does_not_raise(),
            id="one-wavelength-none",
        ),
    ],
)
def test_reflection_eq_branches(r1_kwargs, r2_kwargs, expected_equal, context):
    with context:
        r1 = Reflection(**r1_kwargs)
        r2 = Reflection(**r2_kwargs)
        assert (r1 == r2) == expected_equal


def test_reflection_eq_non_reflection():
    """__eq__ returns NotImplemented for non-Reflection objects."""
    r = Reflection("r", (1, 0, 0), {"omega": 10.0})
    assert r.__eq__("not a reflection") is NotImplemented


# ---------------------------------------------------------------------------
# ReflectionList dict-interface methods
# ---------------------------------------------------------------------------


def test_reflectionlist_keys_values_items():
    """keys(), values(), and items() expose the underlying reflection data."""
    rl = ReflectionList(geometry_name="psic", valid_stages={"omega"})
    rl.add("r1", hkl=(1, 0, 0), angles={"omega": 10.0})
    assert "r1" in list(rl.keys())
    assert all(isinstance(v, Reflection) for v in rl.values())
    assert list(rl.items())[0][0] == "r1"


# ---------------------------------------------------------------------------
# Reflection.to_dict() / from_dict()
# ---------------------------------------------------------------------------

_REFL = Reflection(
    name="r1",
    hkl=(0.0, 0.0, 6.0),
    angles={"omega": 20.97, "chi": 90.0, "phi": 0.0, "ttheta": 41.94},
    wavelength=1.5498,
    geometry_name="fourcv",
)


def test_reflection_to_dict_structure():
    """to_dict() returns a JSON-serialisable dict with required keys."""
    d = _REFL.to_dict()
    assert isinstance(d, dict)
    assert {"name", "hkl", "angles", "wavelength"} <= set(d.keys())
    assert isinstance(d["hkl"], list)
    assert json.dumps(d)  # must not raise


@pytest.mark.parametrize(
    "key, expected, context",
    [
        pytest.param("name", "r1", does_not_raise(), id="name"),
        pytest.param("wavelength", 1.5498, does_not_raise(), id="wavelength"),
    ],
)
def test_reflection_to_dict_values(key, expected, context):
    """to_dict() stores the correct scalar value for each field."""
    with context:
        assert _REFL.to_dict()[key] == pytest.approx(expected)


@pytest.mark.parametrize(
    "attr, accessor, context",
    [
        pytest.param("name", lambda r: r.name, does_not_raise(), id="name"),
        pytest.param("hkl", lambda r: r.hkl, does_not_raise(), id="hkl"),
        pytest.param(
            "wavelength", lambda r: r.wavelength, does_not_raise(), id="wavelength"
        ),
        pytest.param("angles", lambda r: r.angles, does_not_raise(), id="angles"),
    ],
)
def test_reflection_from_dict_roundtrip(attr, accessor, context):
    """from_dict(to_dict()) recovers each attribute."""
    with context:
        restored = Reflection.from_dict(_REFL.to_dict())
        assert accessor(restored) == pytest.approx(accessor(_REFL))


# ---------------------------------------------------------------------------
# ReflectionList.to_dict() / from_dict()
# ---------------------------------------------------------------------------


def _make_rl():
    rl = ReflectionList(
        geometry_name="fourcv",
        valid_stages={"omega", "chi", "phi", "ttheta"},
    )
    rl.add(
        "r1",
        hkl=(0, 0, 6),
        angles={"omega": 20.97, "chi": 90.0, "phi": 0.0, "ttheta": 41.94},
    )
    rl.add(
        "r2",
        hkl=(1, 0, 0),
        angles={"omega": 30.0, "chi": 0.0, "phi": 0.0, "ttheta": 60.0},
    )
    rl.setor1("r1")
    rl.setor2("r2")
    return rl


def test_reflectionlist_to_dict_structure():
    """to_dict() returns a JSON-serialisable dict with a reflections list."""
    d = _make_rl().to_dict()
    assert isinstance(d, dict)
    assert isinstance(d["reflections"], list)
    assert json.dumps(d)  # must not raise


@pytest.mark.parametrize(
    "key, expected, context",
    [
        pytest.param("or1", "r1", does_not_raise(), id="or1"),
        pytest.param("or2", "r2", does_not_raise(), id="or2"),
    ],
)
def test_reflectionlist_to_dict_or_refs(key, expected, context):
    """to_dict() stores or1 and or2 names."""
    with context:
        assert _make_rl().to_dict()[key] == expected


def test_reflectionlist_to_dict_count():
    """to_dict() stores the correct number of reflections."""
    assert len(_make_rl().to_dict()["reflections"]) == 2


@pytest.mark.parametrize(
    "attr, expected, context",
    [
        pytest.param("names", {"r1", "r2"}, does_not_raise(), id="names"),
        pytest.param("or1", "r1", does_not_raise(), id="or1"),
        pytest.param("or2", "r2", does_not_raise(), id="or2"),
        pytest.param("r1-hkl", (0.0, 0.0, 6.0), does_not_raise(), id="r1-hkl"),
    ],
)
def test_reflectionlist_from_dict_roundtrip(attr, expected, context):
    """from_dict(to_dict()) recovers names, or1/or2, and hkl values."""
    with context:
        rl2 = ReflectionList.from_dict(_make_rl().to_dict())
        if attr == "names":
            assert set(rl2._data.keys()) == expected
        elif attr == "or1":
            assert rl2._or1_name == expected
        elif attr == "or2":
            assert rl2._or2_name == expected
        elif attr == "r1-hkl":
            assert rl2["r1"].hkl == pytest.approx(expected)
