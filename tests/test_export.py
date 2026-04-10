"""
Unit tests for to_dict() / from_dict() serialisation (#52).

Covers:
  - Stage.to_dict() / from_dict(): name, axis, role, parent, angle, limits
  - Lattice.to_dict() / from_dict(): all six params, round-trip
  - Reflection.to_dict() / from_dict(): hkl, angles, wavelength, name
  - ReflectionList.to_dict() / from_dict(): ordering, or1/or2 preserved
  - Sample.to_dict() / from_dict(): lattice, reflections, U, UB, name
  - AdHocDiffractometer.to_dict(): delegates stage serialisation to Stage.to_dict();
    _meta keys, JSON-serialisable
  - AdHocDiffractometer.from_dict(): delegates stage construction to Stage.from_dict();
    full round-trip for fourcv and psic
  - Round-trip invariants: name, wavelength, lattice params, reflections,
    stages (names, roles, angles, limits), active sample, U/UB, azimuthal ref
  - JSON round-trip (json.dumps → json.loads → from_dict reproduces geometry)
  - _meta: software name, version string present
"""

import json
import math

import numpy as np
import pytest

from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer import Lattice
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer import psic
from ad_hoc_diffractometer import ub_from_two_reflections_bl1967
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer.reflection import Reflection
from ad_hoc_diffractometer.reflection import ReflectionList
from ad_hoc_diffractometer.sample import Sample
from ad_hoc_diffractometer.stage import Stage

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_TWO_PI = 2.0 * math.pi


def _sapphire_fourcv():
    """fourcv with sapphire, two reflections, UB set."""
    g = fourcv()
    g.wavelength = 1.549802558
    g.azimuthal_reference = (0, 0, 1)
    g.add_sample("sapphire", Lattice(a=4.785, c=12.991, gamma=120))
    g.sample = "sapphire"
    g.add_reflection(
        "or1",
        hkl=(0, 0, 6),
        angles={"omega": 20.97, "chi": 90.0, "phi": 0.0, "two_theta": 41.94188},
    )
    g.add_reflection(
        "or2",
        hkl=(1, 0, 0),
        angles={"omega": 30.0, "chi": 0.0, "phi": 0.0, "two_theta": 60.0},
    )
    g.sample.reflections.setor1("or1")
    g.sample.reflections.setor2("or2")
    ub_from_two_reflections_bl1967(g.sample)
    g.set_angle("omega", 20.97)
    g.set_angle("chi", 90.0)
    return g


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


class TestStageSerialization:
    @pytest.fixture
    def stage(self):
        s = Stage(
            name="omega",
            axis=np.array([0.0, 0.0, -1.0]),
            role="sample",
            parent=None,
            limits=(-180.0, 180.0),
        )
        s.angle = 20.97
        return s

    @pytest.fixture
    def child_stage(self):
        s = Stage(
            name="chi",
            axis=np.array([0.0, 1.0, 0.0]),
            role="sample",
            parent="omega",
            limits=(-90.0, 90.0),
        )
        s.angle = 45.0
        return s

    def test_to_dict_returns_dict(self, stage):
        assert isinstance(stage.to_dict(), dict)

    def test_to_dict_has_required_keys(self, stage):
        d = stage.to_dict()
        assert {"name", "axis", "role", "parent", "angle", "limits"} <= set(d.keys())

    def test_to_dict_name(self, stage):
        assert stage.to_dict()["name"] == "omega"

    def test_to_dict_axis_is_list(self, stage):
        assert isinstance(stage.to_dict()["axis"], list)
        assert len(stage.to_dict()["axis"]) == 3

    def test_to_dict_axis_values(self, stage):
        assert stage.to_dict()["axis"] == pytest.approx([0.0, 0.0, -1.0])

    def test_to_dict_role(self, stage):
        assert stage.to_dict()["role"] == "sample"

    def test_to_dict_parent_none(self, stage):
        assert stage.to_dict()["parent"] is None

    def test_to_dict_parent_name(self, child_stage):
        assert child_stage.to_dict()["parent"] == "omega"

    def test_to_dict_angle(self, stage):
        assert stage.to_dict()["angle"] == pytest.approx(20.97)

    def test_to_dict_limits(self, stage):
        assert stage.to_dict()["limits"] == pytest.approx([-180.0, 180.0])

    def test_to_dict_json_serialisable(self, stage):
        assert json.dumps(stage.to_dict())

    def test_from_dict_roundtrip_name(self, stage):
        assert Stage.from_dict(stage.to_dict()).name == stage.name

    def test_from_dict_roundtrip_axis(self, stage):
        np.testing.assert_allclose(Stage.from_dict(stage.to_dict()).axis, stage.axis)

    def test_from_dict_roundtrip_role(self, stage):
        assert Stage.from_dict(stage.to_dict()).role == stage.role

    def test_from_dict_roundtrip_parent(self, child_stage):
        assert Stage.from_dict(child_stage.to_dict()).parent == child_stage.parent

    def test_from_dict_roundtrip_angle(self, stage):
        assert Stage.from_dict(stage.to_dict()).angle == pytest.approx(stage.angle)

    def test_from_dict_roundtrip_limits(self, stage):
        assert Stage.from_dict(stage.to_dict()).limits == pytest.approx(stage.limits)

    def test_geometry_to_dict_uses_stage_to_dict(self):
        """Each stage dict in geometry.to_dict() matches Stage.to_dict()."""
        g = fourcv()
        g.set_angle("omega", 20.97)
        geo_stages = {sd["name"]: sd for sd in g.to_dict()["stages"]}
        for name, stage_obj in g._stages.items():
            assert geo_stages[name] == stage_obj.to_dict()

    def test_geometry_from_dict_uses_stage_from_dict(self):
        """Stages restored by AdHocDiffractometer.from_dict() match Stage.from_dict()."""
        g = fourcv()
        g.set_angle("chi", 45.0)
        g2 = AdHocDiffractometer.from_dict(g.to_dict())
        for name in g._stages:
            assert g2._stages[name].angle == pytest.approx(g._stages[name].angle)
            np.testing.assert_allclose(g2._stages[name].axis, g._stages[name].axis)
            assert g2._stages[name].role == g._stages[name].role
            assert g2._stages[name].limits == pytest.approx(g._stages[name].limits)


# ---------------------------------------------------------------------------
# Lattice
# ---------------------------------------------------------------------------


class TestLatticeSerialization:
    def test_to_dict_returns_dict(self):
        assert isinstance(Lattice(a=5.0).to_dict(), dict)

    def test_to_dict_has_all_six_keys(self):
        d = Lattice(a=5.0).to_dict()
        assert set(d.keys()) == {"a", "b", "c", "alpha", "beta", "gamma"}

    def test_to_dict_cubic(self):
        d = Lattice(a=5.431).to_dict()
        assert d["a"] == pytest.approx(5.431)
        assert d["b"] == pytest.approx(5.431)
        assert d["c"] == pytest.approx(5.431)
        assert d["alpha"] == pytest.approx(90.0)

    def test_to_dict_hexagonal(self):
        d = Lattice(a=4.785, c=12.991, gamma=120).to_dict()
        assert d["gamma"] == pytest.approx(120.0)
        assert d["a"] == pytest.approx(4.785)
        assert d["c"] == pytest.approx(12.991)

    def test_from_dict_roundtrip_cubic(self):
        lat = Lattice(a=5.431)
        lat2 = Lattice.from_dict(lat.to_dict())
        assert lat == lat2

    def test_from_dict_roundtrip_hexagonal(self):
        lat = Lattice(a=4.785, c=12.991, gamma=120)
        lat2 = Lattice.from_dict(lat.to_dict())
        assert lat == lat2

    def test_from_dict_roundtrip_triclinic(self):
        lat = Lattice(a=5.0, b=6.0, c=7.0, alpha=80.0, beta=90.0, gamma=100.0)
        lat2 = Lattice.from_dict(lat.to_dict())
        assert lat == lat2

    def test_to_dict_json_serialisable(self):
        d = Lattice(a=4.785, c=12.991, gamma=120).to_dict()
        assert json.dumps(d)


# ---------------------------------------------------------------------------
# Reflection
# ---------------------------------------------------------------------------


class TestReflectionSerialization:
    @pytest.fixture
    def refl(self):
        return Reflection(
            name="r1",
            hkl=(0.0, 0.0, 6.0),
            angles={"omega": 20.97, "chi": 90.0, "phi": 0.0, "two_theta": 41.94},
            wavelength=1.5498,
            geometry_name="fourcv",
        )

    def test_to_dict_returns_dict(self, refl):
        assert isinstance(refl.to_dict(), dict)

    def test_to_dict_has_required_keys(self, refl):
        d = refl.to_dict()
        assert "name" in d and "hkl" in d and "angles" in d

    def test_to_dict_hkl_is_list(self, refl):
        assert isinstance(refl.to_dict()["hkl"], list)

    def test_to_dict_wavelength(self, refl):
        assert refl.to_dict()["wavelength"] == pytest.approx(1.5498)

    def test_roundtrip_name(self, refl):
        assert Reflection.from_dict(refl.to_dict()).name == refl.name

    def test_roundtrip_hkl(self, refl):
        assert Reflection.from_dict(refl.to_dict()).hkl == pytest.approx(refl.hkl)

    def test_roundtrip_angles(self, refl):
        d = Reflection.from_dict(refl.to_dict())
        assert d.angles == pytest.approx(refl.angles)

    def test_roundtrip_wavelength(self, refl):
        d = Reflection.from_dict(refl.to_dict())
        assert d.wavelength == pytest.approx(refl.wavelength)

    def test_to_dict_json_serialisable(self, refl):
        assert json.dumps(refl.to_dict())


# ---------------------------------------------------------------------------
# ReflectionList
# ---------------------------------------------------------------------------


class TestReflectionListSerialization:
    @pytest.fixture
    def rl(self):
        r = ReflectionList(
            geometry_name="fourcv",
            valid_stages={"omega", "chi", "phi", "two_theta"},
        )
        r.add(
            "r1",
            hkl=(0, 0, 6),
            angles={"omega": 20.97, "chi": 90.0, "phi": 0.0, "two_theta": 41.94},
        )
        r.add(
            "r2",
            hkl=(1, 0, 0),
            angles={"omega": 30.0, "chi": 0.0, "phi": 0.0, "two_theta": 60.0},
        )
        r.setor1("r1")
        r.setor2("r2")
        return r

    def test_to_dict_returns_dict(self, rl):
        assert isinstance(rl.to_dict(), dict)

    def test_to_dict_has_reflections_list(self, rl):
        assert isinstance(rl.to_dict()["reflections"], list)

    def test_to_dict_reflection_count(self, rl):
        assert len(rl.to_dict()["reflections"]) == 2

    def test_to_dict_or1_or2_stored(self, rl):
        d = rl.to_dict()
        assert d["or1"] == "r1"
        assert d["or2"] == "r2"

    def test_roundtrip_names(self, rl):
        rl2 = ReflectionList.from_dict(rl.to_dict())
        assert set(rl2._data.keys()) == {"r1", "r2"}

    def test_roundtrip_or1_or2(self, rl):
        rl2 = ReflectionList.from_dict(rl.to_dict())
        assert rl2._or1_name == "r1"
        assert rl2._or2_name == "r2"

    def test_roundtrip_hkl(self, rl):
        rl2 = ReflectionList.from_dict(rl.to_dict())
        assert rl2["r1"].hkl == pytest.approx((0.0, 0.0, 6.0))

    def test_to_dict_json_serialisable(self, rl):
        assert json.dumps(rl.to_dict())


# ---------------------------------------------------------------------------
# Sample
# ---------------------------------------------------------------------------


class TestSampleSerialization:
    @pytest.fixture
    def sample(self):
        rl = ReflectionList(
            geometry_name="fourcv",
            valid_stages={"omega", "chi", "phi", "two_theta"},
        )
        rl.add(
            "r1",
            hkl=(0, 0, 6),
            angles={"omega": 20.97, "chi": 90.0, "phi": 0.0, "two_theta": 41.94},
        )
        s = Sample(
            name="sapphire",
            lattice=Lattice(a=4.785, c=12.991, gamma=120),
            reflections=rl,
        )
        s.U = np.eye(3)
        s.UB = np.eye(3) * 0.5
        return s

    def test_to_dict_returns_dict(self, sample):
        assert isinstance(sample.to_dict(), dict)

    def test_to_dict_has_required_keys(self, sample):
        d = sample.to_dict()
        assert {"name", "lattice", "reflections", "U", "UB"} <= set(d.keys())

    def test_to_dict_U_is_list(self, sample):
        assert isinstance(sample.to_dict()["U"], list)

    def test_roundtrip_name(self, sample):
        s2 = Sample.from_dict(sample.to_dict())
        assert s2.name == sample.name

    def test_roundtrip_lattice(self, sample):
        s2 = Sample.from_dict(sample.to_dict())
        assert s2.lattice == sample.lattice

    def test_roundtrip_U(self, sample):
        s2 = Sample.from_dict(sample.to_dict())
        np.testing.assert_allclose(s2.U, sample.U)

    def test_roundtrip_UB(self, sample):
        s2 = Sample.from_dict(sample.to_dict())
        np.testing.assert_allclose(s2.UB, sample.UB)

    def test_roundtrip_U_none(self):
        rl = ReflectionList(geometry_name="g", valid_stages=set())
        s = Sample(name="t", lattice=Lattice(a=1.0), reflections=rl)
        s2 = Sample.from_dict(s.to_dict())
        assert s2.U is None
        assert s2.UB is None

    def test_to_dict_json_serialisable(self, sample):
        assert json.dumps(sample.to_dict())


# ---------------------------------------------------------------------------
# AdHocDiffractometer.to_dict()
# ---------------------------------------------------------------------------


class TestGeometryToDict:
    @pytest.fixture
    def geom(self):
        return _sapphire_fourcv()

    def test_returns_dict(self, geom):
        assert isinstance(geom.to_dict(), dict)

    def test_json_serialisable(self, geom):
        d = geom.to_dict()
        assert json.dumps(d)  # must not raise

    def test_meta_software_key(self, geom):
        assert geom.to_dict()["_meta"]["software"] == "ad_hoc_diffractometer"

    def test_meta_version_key(self, geom):
        v = geom.to_dict()["_meta"]["version"]
        assert isinstance(v, str) and len(v) > 0

    def test_meta_created_key(self, geom):
        created = geom.to_dict()["_meta"]["created"]
        assert isinstance(created, str) and "T" in created  # ISO-8601

    def test_name(self, geom):
        assert geom.to_dict()["name"] == "fourcv"

    def test_wavelength(self, geom):
        assert geom.to_dict()["wavelength"] == pytest.approx(1.549802558)

    def test_azimuthal_reference(self, geom):
        assert geom.to_dict()["azimuthal_reference"] == pytest.approx([0.0, 0.0, 1.0])

    def test_stages_list(self, geom):
        stages = geom.to_dict()["stages"]
        assert isinstance(stages, list)
        assert {s["name"] for s in stages} == {"omega", "chi", "phi", "two_theta"}

    def test_stage_has_required_keys(self, geom):
        for sd in geom.to_dict()["stages"]:
            assert {"name", "axis", "role", "parent", "angle", "limits"} <= set(
                sd.keys()
            )

    def test_stage_angle_preserved(self, geom):
        """omega was set to 20.97 before to_dict()."""
        stages = {s["name"]: s for s in geom.to_dict()["stages"]}
        assert stages["omega"]["angle"] == pytest.approx(20.97)

    def test_samples_dict(self, geom):
        d = geom.to_dict()
        assert "test" in d["samples"]
        assert "sapphire" in d["samples"]

    def test_active_sample(self, geom):
        assert geom.to_dict()["active_sample"] == "sapphire"

    def test_basis_stored(self, geom):
        basis = geom.to_dict()["basis"]
        assert isinstance(basis, dict)
        assert len(basis) == 3


# ---------------------------------------------------------------------------
# AdHocDiffractometer.from_dict() — full round-trip
# ---------------------------------------------------------------------------


class TestGeometryRoundTrip:
    @pytest.fixture
    def original(self):
        return _sapphire_fourcv()

    @pytest.fixture
    def restored(self, original):
        return AdHocDiffractometer.from_dict(original.to_dict())

    def test_name(self, original, restored):
        assert restored.name == original.name

    def test_wavelength(self, original, restored):
        assert restored.wavelength == pytest.approx(original.wavelength)

    def test_azimuthal_reference(self, original, restored):
        assert restored.azimuthal_reference == pytest.approx(
            original.azimuthal_reference
        )

    def test_stage_names(self, original, restored):
        assert set(restored._stages.keys()) == set(original._stages.keys())

    def test_stage_roles(self, original, restored):
        for name in original._stages:
            assert restored._stages[name].role == original._stages[name].role

    def test_stage_angles(self, original, restored):
        for name in original._stages:
            assert restored._stages[name].angle == pytest.approx(
                original._stages[name].angle
            )

    def test_stage_limits(self, original, restored):
        for name in original._stages:
            assert restored._stages[name].limits == pytest.approx(
                original._stages[name].limits
            )

    def test_active_sample_name(self, original, restored):
        assert restored._active_ref[0] == original._active_ref[0]

    def test_sample_names(self, original, restored):
        assert set(restored.samples._data.keys()) == set(original.samples._data.keys())

    def test_sample_lattice_a(self, original, restored):
        assert restored.sample.lattice.a == pytest.approx(original.sample.lattice.a)

    def test_sample_lattice_c(self, original, restored):
        assert restored.sample.lattice.c == pytest.approx(original.sample.lattice.c)

    def test_sample_lattice_gamma(self, original, restored):
        assert restored.sample.lattice.gamma == pytest.approx(
            original.sample.lattice.gamma
        )

    def test_sample_reflections_names(self, original, restored):
        orig_names = set(original.sample.reflections._data.keys())
        rest_names = set(restored.sample.reflections._data.keys())
        assert rest_names == orig_names

    def test_sample_or1(self, original, restored):
        assert restored.sample.reflections._or1_name == "or1"

    def test_sample_or2(self, original, restored):
        assert restored.sample.reflections._or2_name == "or2"

    def test_sample_UB_preserved(self, original, restored):
        np.testing.assert_allclose(restored.sample.UB, original.sample.UB, atol=1e-12)

    def test_sample_U_preserved(self, original, restored):
        np.testing.assert_allclose(restored.sample.U, original.sample.U, atol=1e-12)

    def test_json_roundtrip(self, original):
        """Full JSON serialize → deserialize → from_dict round-trip."""
        d = json.loads(json.dumps(original.to_dict()))
        g2 = AdHocDiffractometer.from_dict(d)
        assert g2.name == original.name
        assert g2.wavelength == pytest.approx(original.wavelength)
        assert g2.sample.lattice.a == pytest.approx(original.sample.lattice.a)

    def test_psic_roundtrip(self):
        """Round-trip for psic geometry (6 stages, different basis)."""
        g = psic()
        g.wavelength = _TWO_PI
        ub_identity(g.sample)
        d = g.to_dict()
        g2 = AdHocDiffractometer.from_dict(d)
        assert g2.name == "psic"
        assert set(g2._stages.keys()) == set(g._stages.keys())

    def test_no_azimuthal_reference(self):
        """Geometry with no azimuthal reference round-trips cleanly."""
        g = fourcv()
        g.wavelength = 1.5406
        g2 = AdHocDiffractometer.from_dict(g.to_dict())
        assert g2.azimuthal_reference is None

    def test_kappa_alpha_deg_preserved(self):
        """kappa_alpha_deg is preserved in the round-trip."""
        from ad_hoc_diffractometer import kappa4cv

        g = kappa4cv(alpha_deg=50.0)
        g.wavelength = 1.5406
        g2 = AdHocDiffractometer.from_dict(g.to_dict())
        assert g2.kappa_alpha_deg == pytest.approx(50.0)

    def test_default_test_sample_removed_when_not_in_export(self):
        """
        The default 'test' sample created by __init__ must not appear in the
        restored geometry if it was not present in the exported dict.
        """
        g = fourcv()
        g.wavelength = 1.5406
        # Add a named sample and make it active, then remove the default
        g.add_sample("sapphire", Lattice(a=4.785, c=12.991, gamma=120))
        g.sample = "sapphire"
        g.remove_sample("test")  # 'test' is gone before export
        assert "test" not in g.samples._data

        d = g.to_dict()
        g2 = AdHocDiffractometer.from_dict(d)
        # The default 'test' sample must NOT have crept back in
        assert "test" not in g2.samples._data
        assert "sapphire" in g2.samples._data

    def test_only_saved_samples_present_after_restore(self):
        """
        After from_dict(), the sample dict contains exactly the samples
        that were saved — no extras, no stale defaults.
        """
        g = _sapphire_fourcv()
        d = g.to_dict()
        g2 = AdHocDiffractometer.from_dict(d)
        assert set(g2.samples._data.keys()) == set(d["samples"].keys())

    def test_active_sample_pointer_correct_after_restore(self):
        """
        The active-sample pointer must reference the restored active name,
        not the stale default, even when 'test' is removed during restore.
        """
        g = fourcv()
        g.wavelength = 1.5406
        g.add_sample("mycrystal", Lattice(a=5.0))
        g.sample = "mycrystal"
        g.remove_sample("test")

        g2 = AdHocDiffractometer.from_dict(g.to_dict())
        assert g2._active_ref[0] == "mycrystal"
        assert g2.sample.name == "mycrystal"
