"""
Unit tests for ad_hoc_diffractometer.reflection.

Covers:
  - Reflection dataclass: construction, normalisation, validation, __eq__, __repr__
  - ReflectionList: add, remove, clear, dict-like interface,
    setor1/setor2, orienting_reflections, cross-geometry safety
  - AdHocDiffractometer.add_reflection() convenience wrapper
"""

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


def test_add_reflection_stored_in_reflections(psic_geom):
    psic_geom.add_reflection("r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES)
    assert "r1" in psic_geom.reflections


def test_add_reflection_geometry_name_set(psic_geom):
    r = psic_geom.add_reflection("r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES)
    assert r.geometry_name == "psic"
