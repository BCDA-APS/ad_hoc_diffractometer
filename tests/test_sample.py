"""
Unit tests for ad_hoc_diffractometer.sample.

Covers:
  - Sample construction, __repr__, __eq__
  - SampleDict: type guard, active-sample guard, pop/clear/del protection,
    direct _data assignment bypass blocked by read-only _samples property
  - AdHocDiffractometer.samples, .sample, add_sample(), remove_sample()
  - add_reflection() targets the active sample
  - Switching active sample
"""

import re

import numpy as np
import pytest

from ad_hoc_diffractometer import Lattice
from ad_hoc_diffractometer import Sample

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
