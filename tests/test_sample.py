# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.sample.

Covers:
  - Sample construction, __repr__, __eq__
  - SampleDict: type guard, active-sample guard, pop/clear/del protection,
    direct _data assignment bypass blocked by read-only _samples property
  - AdHocDiffractometer.samples, .sample, add_sample(), remove_sample()
  - add_reflection() targets the active sample
  - Switching active sample
  - Sample.to_dict() / from_dict(): lattice, reflections, U, UB, name
"""

import json
import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

from ad_hoc_diffractometer import Lattice
from ad_hoc_diffractometer import Sample
from ad_hoc_diffractometer.reflection import ReflectionList

_PSIC_ANGLES = {
    "mu": 0.0,
    "eta": 20.0,
    "chi": 90.0,
    "phi": 0.0,
    "nu": 0.0,
    "delta": 40.0,
}


# ---------------------------------------------------------------------------
# Sample dataclass
# ---------------------------------------------------------------------------


def test_sample_default_U_UB_none(psic_geom):
    assert psic_geom.sample.U is None
    assert psic_geom.sample.UB is None


def test_sample_name(psic_geom):
    assert psic_geom.sample.name == "test"


def test_sample_default_lattice_cubic_1angstrom(psic_geom):
    lat = psic_geom.sample.lattice
    assert lat.a == pytest.approx(1.0)
    assert lat.b == pytest.approx(1.0)
    assert lat.c == pytest.approx(1.0)
    assert lat.system == "cubic"


def test_sample_repr(psic_geom):
    r = repr(psic_geom.sample)
    assert "test" in r
    assert "U=None" in r
    assert "UB=None" in r


def test_sample_eq_identical():
    from ad_hoc_diffractometer.reflection import ReflectionList

    rl = ReflectionList(geometry_name="psic", valid_stages=set())
    s1 = Sample(name="a", lattice=Lattice(a=1.0), reflections=rl)
    s2 = Sample(name="a", lattice=Lattice(a=1.0), reflections=rl)
    assert s1 == s2


def test_sample_eq_different_name():
    from ad_hoc_diffractometer.reflection import ReflectionList

    rl = ReflectionList(geometry_name="psic", valid_stages=set())
    s1 = Sample(name="a", lattice=Lattice(a=1.0), reflections=rl)
    s2 = Sample(name="b", lattice=Lattice(a=1.0), reflections=rl)
    assert s1 != s2


def test_sample_eq_not_implemented_for_non_sample():
    from ad_hoc_diffractometer.reflection import ReflectionList

    rl = ReflectionList(geometry_name="psic", valid_stages=set())
    s = Sample(name="a", lattice=Lattice(a=1.0), reflections=rl)
    assert s.__eq__("not a sample") is NotImplemented


# ---------------------------------------------------------------------------
# SampleDict guards
# ---------------------------------------------------------------------------


def test_sample_dict_rejects_none(psic_geom):
    with pytest.raises(TypeError, match=re.escape("only accepts Sample")):
        psic_geom.samples["test"] = None  # type: ignore


def test_sample_dict_rejects_arbitrary_object(psic_geom):
    with pytest.raises(TypeError, match=re.escape("only accepts Sample")):
        psic_geom.samples["test"] = 42  # type: ignore


def test_sample_dict_rejects_replace_active(psic_geom):
    """Replacing the active sample's value is blocked."""
    from ad_hoc_diffractometer.reflection import ReflectionList

    rl = ReflectionList(geometry_name="psic", valid_stages=set())
    impostor = Sample(name="test", lattice=Lattice(a=2.0), reflections=rl)
    with pytest.raises(ValueError, match=re.escape("Cannot replace")):
        psic_geom.samples["test"] = impostor


def test_sample_dict_del_active_raises(psic_geom):
    with pytest.raises(ValueError, match=re.escape("Cannot remove")):
        del psic_geom.samples["test"]


def test_sample_dict_pop_active_raises(psic_geom):
    with pytest.raises(ValueError, match=re.escape("Cannot pop")):
        psic_geom.samples.pop("test")


def test_sample_dict_clear_raises(psic_geom):
    with pytest.raises(ValueError, match=re.escape("not permitted")):
        psic_geom.samples.clear()


def test_sample_dict_replace_samples_attr_raises(psic_geom):
    """Reassigning _samples or samples is blocked (read-only property)."""
    with pytest.raises(AttributeError):
        psic_geom._samples = {}  # type: ignore


def test_sample_dict_pop_nonexistent_returns_default(psic_geom):
    result = psic_geom.samples.pop("nonexistent", "fallback")
    assert result == "fallback"


def test_sample_dict_pop_nonexistent_raises_without_default(psic_geom):
    with pytest.raises(KeyError):
        psic_geom.samples.pop("nonexistent")


def test_sample_dict_del_non_active_succeeds(psic_geom):
    psic_geom.add_sample("silicon", Lattice(a=5.431))
    del psic_geom.samples["silicon"]
    assert "silicon" not in psic_geom.samples


# ---------------------------------------------------------------------------
# AdHocDiffractometer.samples / .sample
# ---------------------------------------------------------------------------


def test_samples_contains_test_by_default(psic_geom):
    assert "test" in psic_geom.samples


def test_samples_length_default(psic_geom):
    assert len(psic_geom.samples) == 1


def test_sample_property_returns_active(psic_geom):
    assert psic_geom.sample.name == "test"


def test_sample_setter_by_name(psic_geom):
    psic_geom.add_sample("silicon", Lattice(a=5.431))
    psic_geom.sample = "silicon"
    assert psic_geom.sample.name == "silicon"


def test_sample_setter_by_object(psic_geom):
    s = psic_geom.add_sample("silicon", Lattice(a=5.431))
    psic_geom.sample = s
    assert psic_geom.sample.name == "silicon"


def test_sample_setter_unknown_raises(psic_geom):
    with pytest.raises(KeyError, match="missing"):
        psic_geom.sample = "missing"


# ---------------------------------------------------------------------------
# add_sample() / remove_sample()
# ---------------------------------------------------------------------------


def test_add_sample_returns_sample(psic_geom):
    s = psic_geom.add_sample("silicon", Lattice(a=5.431))
    assert isinstance(s, Sample)
    assert s.name == "silicon"


def test_add_sample_default_lattice(psic_geom):
    s = psic_geom.add_sample("unknown")
    assert s.lattice.a == pytest.approx(1.0)
    assert s.lattice.system == "cubic"


def test_add_sample_duplicate_raises(psic_geom):
    psic_geom.add_sample("silicon", Lattice(a=5.431))
    with pytest.raises(ValueError, match=re.escape("already exists")):
        psic_geom.add_sample("silicon", Lattice(a=5.431))


def test_add_sample_reflections_use_geometry_stages(psic_geom):
    s = psic_geom.add_sample("silicon", Lattice(a=5.431))
    assert s.reflections.valid_stages == set(psic_geom._stages)


def test_remove_sample_non_active(psic_geom):
    psic_geom.add_sample("silicon", Lattice(a=5.431))
    psic_geom.remove_sample("silicon")
    assert "silicon" not in psic_geom.samples


def test_remove_sample_active_raises(psic_geom):
    with pytest.raises(ValueError, match=re.escape("Cannot remove")):
        psic_geom.remove_sample("test")


def test_remove_sample_then_reselect(psic_geom):
    psic_geom.add_sample("silicon", Lattice(a=5.431))
    psic_geom.sample = "silicon"
    psic_geom.sample = "test"
    psic_geom.remove_sample("silicon")
    assert "silicon" not in psic_geom.samples
    assert psic_geom.sample.name == "test"


# ---------------------------------------------------------------------------
# add_reflection() targets the active sample
# ---------------------------------------------------------------------------


def test_add_reflection_goes_to_active_sample(psic_geom):
    psic_geom.add_reflection("r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES)
    assert "r1" in psic_geom.sample.reflections


def test_add_reflection_not_in_other_sample(psic_geom):
    psic_geom.add_sample("silicon", Lattice(a=5.431))
    psic_geom.add_reflection("r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES)
    # r1 is in "test", not in "silicon"
    assert "r1" in psic_geom.samples["test"].reflections
    assert "r1" not in psic_geom.samples["silicon"].reflections


def test_add_reflection_switches_with_active_sample(psic_geom):
    psic_geom.add_sample("silicon", Lattice(a=5.431))
    psic_geom.add_reflection("r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES)
    psic_geom.sample = "silicon"
    psic_geom.add_reflection("s1", hkl=(1, 1, 1), angles=_PSIC_ANGLES)
    assert "r1" in psic_geom.samples["test"].reflections
    assert "s1" in psic_geom.samples["silicon"].reflections
    assert "s1" not in psic_geom.samples["test"].reflections


def test_add_reflection_inherits_wavelength_from_geometry(psic_geom):
    psic_geom.wavelength = 1.5406
    r = psic_geom.add_reflection("r1", hkl=(1, 0, 0), angles=_PSIC_ANGLES)
    assert r.wavelength == pytest.approx(1.5406)


# ---------------------------------------------------------------------------
# U / UB on Sample
# ---------------------------------------------------------------------------


def test_sample_U_can_be_set(psic_geom):
    psic_geom.sample.U = np.eye(3)
    np.testing.assert_array_equal(psic_geom.sample.U, np.eye(3))


def test_sample_UB_can_be_set(psic_geom):
    psic_geom.sample.UB = np.eye(3)
    np.testing.assert_array_equal(psic_geom.sample.UB, np.eye(3))


def test_sample_U_UB_independent_per_sample(psic_geom):
    psic_geom.add_sample("silicon", Lattice(a=5.431))
    psic_geom.sample.U = np.eye(3)  # set on "test"
    psic_geom.sample = "silicon"
    assert psic_geom.sample.U is None  # "silicon" still has None


# ---------------------------------------------------------------------------
# Sample.parent
# ---------------------------------------------------------------------------


def test_default_sample_parent_is_geometry(psic_geom):
    """Default 'test' sample's parent is the owning geometry."""
    assert psic_geom.sample.parent is psic_geom


def test_add_sample_sets_parent(psic_geom):
    s = psic_geom.add_sample("silicon", Lattice(a=5.431))
    assert s.parent is psic_geom


def test_remove_sample_clears_parent(psic_geom):
    s = psic_geom.add_sample("silicon", Lattice(a=5.431))
    psic_geom.remove_sample("silicon")
    assert s.parent is None


def test_standalone_sample_parent_is_none():
    """Sample constructed without a geometry has parent=None."""
    from ad_hoc_diffractometer.reflection import ReflectionList

    rl = ReflectionList(geometry_name="test", valid_stages=set())
    s = Sample(name="standalone", lattice=Lattice(a=1.0), reflections=rl)
    assert s.parent is None


def test_sample_parent_excluded_from_eq(psic_geom):
    """Two samples with same content but different parents compare equal."""
    from ad_hoc_diffractometer.reflection import ReflectionList

    rl = ReflectionList(geometry_name="test", valid_stages=set())
    s_standalone = Sample(name="test", lattice=Lattice(a=1.0), reflections=rl)
    # psic_geom.sample has parent=psic_geom; standalone has parent=None
    # They should still be equal (content-based)
    assert psic_geom.sample == s_standalone


def test_sample_repr_shows_geometry_name(psic_geom):
    assert "psic" in repr(psic_geom.sample)


def test_sample_repr_shows_no_parent_for_standalone():
    from ad_hoc_diffractometer.reflection import ReflectionList

    rl = ReflectionList(geometry_name="test", valid_stages=set())
    s = Sample(name="standalone", lattice=Lattice(a=1.0), reflections=rl)
    assert "(no parent)" in repr(s)


# ---------------------------------------------------------------------------
# SampleDict guard branches
# ---------------------------------------------------------------------------


def _two_sample_geom():
    """fourcv with 'test' (active) and 's2' samples."""
    from helpers import fourcv

    g = fourcv()
    g.add_sample("s2", g.sample.lattice)
    return g


@pytest.mark.parametrize(
    "op, context",
    [
        pytest.param(
            "replace_active",
            pytest.raises(ValueError, match="active"),
            id="replace-active-sample",
        ),
        pytest.param(
            "del_missing",
            pytest.raises(KeyError),
            id="del-missing-sample",
        ),
        pytest.param(
            "clear",
            pytest.raises(ValueError, match="clear"),
            id="clear-not-permitted",
        ),
        pytest.param(
            "pop_missing_no_default",
            pytest.raises(KeyError),
            id="pop-missing-no-default",
        ),
        pytest.param(
            "pop_active",
            pytest.raises(ValueError, match="active"),
            id="pop-active",
        ),
    ],
)
def test_sampledict_guards(op, context):
    """SampleDict enforces type and active-sample invariants."""
    from ad_hoc_diffractometer.reflection import ReflectionList
    from ad_hoc_diffractometer.sample import Sample as _Sample

    g = _two_sample_geom()
    with context:
        if op == "replace_active":
            new = _Sample(
                name="test",
                lattice=Lattice(a=5.0),
                reflections=ReflectionList(geometry_name="fourcv", valid_stages=set()),
            )
            g.samples["test"] = new
        elif op == "del_missing":
            del g.samples["no_such"]
        elif op == "clear":
            g.samples.clear()
        elif op == "pop_missing_no_default":
            g.samples.pop("no_such")
        elif op == "pop_active":
            g.samples.pop("test")


def test_sampledict_replace_non_active():
    """SampleDict.__setitem__ replaces a non-active sample without error."""
    from ad_hoc_diffractometer.reflection import ReflectionList
    from ad_hoc_diffractometer.sample import Sample as _Sample

    g = _two_sample_geom()
    g.sample = "s2"
    new = _Sample(
        name="test",
        lattice=Lattice(a=5.0),
        reflections=ReflectionList(geometry_name="fourcv", valid_stages=set()),
    )
    g.samples["test"] = new
    assert g.samples["test"].lattice.a == pytest.approx(5.0)


def test_sampledict_pop_missing_with_default():
    """SampleDict.pop() returns the supplied default when name is missing."""
    g = _two_sample_geom()
    assert g.samples.pop("no_such", "sentinel") == "sentinel"


def test_sampledict_pop_non_active_succeeds():
    """SampleDict.pop() removes and returns a non-active sample."""
    g = _two_sample_geom()
    g.sample = "s2"  # make s2 active, so 'test' is non-active
    removed = g.samples.pop("test")
    assert removed.name == "test"
    assert "test" not in g.samples


def test_sampledict_iter():
    """SampleDict supports iteration over sample names."""
    g = _two_sample_geom()
    names = list(g.samples)
    assert "test" in names
    assert "s2" in names


def test_sampledict_repr():
    """SampleDict.__repr__ lists the sample names."""
    g = _two_sample_geom()
    r = repr(g.samples)
    assert "test" in r


def test_sampledict_keys_values_items():
    """SampleDict.keys(), values(), and items() expose the underlying data."""
    g = _two_sample_geom()
    assert "test" in g.samples.keys()
    assert "s2" in g.samples.keys()
    assert all(hasattr(v, "lattice") for v in g.samples.values())
    assert any(k == "test" for k, _ in g.samples.items())


# ---------------------------------------------------------------------------
# Sample.to_dict() / from_dict()
# ---------------------------------------------------------------------------


def _make_sample():
    rl = ReflectionList(
        geometry_name="fourcv", valid_stages={"omega", "chi", "phi", "ttheta"}
    )
    rl.add(
        "r1",
        hkl=(0, 0, 6),
        angles={"omega": 20.97, "chi": 90.0, "phi": 0.0, "ttheta": 41.94},
    )
    s = Sample(
        name="sapphire", lattice=Lattice(a=4.785, c=12.991, gamma=120), reflections=rl
    )
    s.U = np.eye(3)
    s.UB = np.eye(3) * 0.5
    return s


def test_sample_to_dict_structure():
    """to_dict() returns a JSON-serialisable dict with all required keys."""
    d = _make_sample().to_dict()
    assert isinstance(d, dict)
    assert {"name", "lattice", "reflections", "U", "UB"} <= set(d.keys())
    assert isinstance(d["U"], list)
    assert json.dumps(d)  # must not raise


@pytest.mark.parametrize(
    "attr, accessor, context",
    [
        pytest.param("name", lambda s: s.name, does_not_raise(), id="name"),
        pytest.param("lattice", lambda s: s.lattice, does_not_raise(), id="lattice"),
    ],
)
def test_sample_from_dict_roundtrip_scalars(attr, accessor, context):
    """from_dict(to_dict()) recovers name and lattice."""
    with context:
        s = _make_sample()
        s2 = Sample.from_dict(s.to_dict())
        assert accessor(s2) == pytest.approx(accessor(s))


@pytest.mark.parametrize(
    "matrix_attr, context",
    [
        pytest.param("U", does_not_raise(), id="U"),
        pytest.param("UB", does_not_raise(), id="UB"),
    ],
)
def test_sample_from_dict_roundtrip_matrices(matrix_attr, context):
    """from_dict(to_dict()) recovers U and UB matrices."""
    with context:
        s = _make_sample()
        s2 = Sample.from_dict(s.to_dict())
        np.testing.assert_allclose(getattr(s2, matrix_attr), getattr(s, matrix_attr))


def test_sample_from_dict_roundtrip_none_matrices():
    """from_dict(to_dict()) correctly restores None for unset U/UB."""
    rl = ReflectionList(geometry_name="g", valid_stages=set())
    s = Sample(name="t", lattice=Lattice(a=1.0), reflections=rl)
    s2 = Sample.from_dict(s.to_dict())
    with does_not_raise():
        assert s2.U is None
        assert s2.UB is None
