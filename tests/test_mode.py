# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.mode and mode-related geometry features.

Covers:
  - DiffractionMode ABC (constrained_stages, apply_cut_point, __repr__, __eq__)
  - FixedAngleMode: construction, frozen_angles, constrained_stages, __repr__, __eq__
  - BisectingMode: construction, constrained_stages, __repr__, __eq__
  - ModeDict: construction, set/get/delete, type guard, __len__, __iter__,
    keys/values/items, __repr__, __eq__, __contains__
  - AdHocDiffractometer: modes, default_mode, mode_name, mode property,
    cut_points, mode_name setter validation, mode=None when no mode set
  - Factory modes: psic, fourcv, fourch, kappa4cv, kappa4ch, kappa6c
    each pre-populate modes and have a bisecting default
  - Serialisation round-trip: to_dict / from_dict for FixedAngleMode,
    BisectingMode, ModeDict on geometry, and cut_points
  - Export/import: mode state survives to_dict / from_dict
"""

import json
import re
from contextlib import nullcontext as does_not_raise

import pytest

from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer import BisectingMode
from ad_hoc_diffractometer import FixedAngleMode
from ad_hoc_diffractometer import ModeDict
from ad_hoc_diffractometer import Stage
from ad_hoc_diffractometer import fivec
from ad_hoc_diffractometer import fourch
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer import kappa4ch
from ad_hoc_diffractometer import kappa4cv
from ad_hoc_diffractometer import kappa6c
from ad_hoc_diffractometer import psic
from ad_hoc_diffractometer import s2d2
from ad_hoc_diffractometer import sixc
from ad_hoc_diffractometer import zaxis
from ad_hoc_diffractometer.constants import XHAT
from ad_hoc_diffractometer.constants import YHAT
from ad_hoc_diffractometer.constants import ZHAT

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _simple_geometry(**kwargs):
    """Minimal 2-stage geometry for testing mode features."""
    stages = [
        Stage("omega", -ZHAT, parent=None, role="sample"),
        Stage("ttheta", -ZHAT, parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name="test_geom",
        stages=stages,
        basis={"vertical": XHAT, "longitudinal": YHAT, "lateral": ZHAT},
        **kwargs,
    )


# ---------------------------------------------------------------------------
# FixedAngleMode
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "stage, value, cut_points, expected_frozen, expected_constrained, context",
    [
        pytest.param(
            "chi",
            90.0,
            None,
            {"chi": 90.0},
            ["chi"],
            does_not_raise(),
            id="basic-fixed-chi-90",
        ),
        pytest.param(
            "phi",
            0.0,
            None,
            {"phi": 0.0},
            ["phi"],
            does_not_raise(),
            id="basic-fixed-phi-0",
        ),
        pytest.param(
            "mu",
            -45.0,
            {"mu": -180.0},
            {"mu": -45.0},
            ["mu"],
            does_not_raise(),
            id="with-cut-point",
        ),
    ],
)
def test_fixed_angle_mode_construction(
    stage, value, cut_points, expected_frozen, expected_constrained, context
):
    with context:
        mode = FixedAngleMode(stage=stage, value=value, cut_points=cut_points)
        assert mode.frozen_angles == expected_frozen
        assert mode.constrained_stages == expected_constrained
        assert mode.cut_points == (cut_points or {})


def test_fixed_angle_mode_repr():
    mode = FixedAngleMode(stage="chi", value=90.0)
    r = repr(mode)
    assert "FixedAngleMode" in r
    assert "chi" in r
    assert "90.0" in r


def test_fixed_angle_mode_eq_same():
    m1 = FixedAngleMode(stage="chi", value=90.0)
    m2 = FixedAngleMode(stage="chi", value=90.0)
    assert m1 == m2


def test_fixed_angle_mode_eq_different_value():
    m1 = FixedAngleMode(stage="chi", value=90.0)
    m2 = FixedAngleMode(stage="chi", value=0.0)
    assert m1 != m2


def test_fixed_angle_mode_eq_different_stage():
    m1 = FixedAngleMode(stage="chi", value=90.0)
    m2 = FixedAngleMode(stage="phi", value=90.0)
    assert m1 != m2


def test_fixed_angle_mode_eq_different_type():
    m1 = FixedAngleMode(stage="chi", value=90.0)
    m2 = BisectingMode(sample_stage="omega", detector_stage="ttheta")
    assert m1 != m2


def test_fixed_angle_mode_eq_non_mode():
    m = FixedAngleMode(stage="chi", value=90.0)
    assert m != "not a mode"


# ---------------------------------------------------------------------------
# BisectingMode
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "sample_stage, detector_stage, frozen_angles, cut_points, "
    "expected_constrained, context",
    [
        pytest.param(
            "omega",
            "ttheta",
            None,
            None,
            ["omega"],
            does_not_raise(),
            id="basic-bisecting-omega-ttheta",
        ),
        pytest.param(
            "eta",
            "delta",
            {"mu": 0.0, "nu": 0.0},
            None,
            ["eta", "mu", "nu"],
            does_not_raise(),
            id="psic-bisecting-with-frozen",
        ),
        pytest.param(
            "komega",
            "delta",
            {"mu": 0.0, "nu": 0.0},
            {"komega": -180.0},
            ["komega", "mu", "nu"],
            does_not_raise(),
            id="bisecting-with-cut-point",
        ),
    ],
)
def test_bisecting_mode_construction(
    sample_stage,
    detector_stage,
    frozen_angles,
    cut_points,
    expected_constrained,
    context,
):
    with context:
        mode = BisectingMode(
            sample_stage=sample_stage,
            detector_stage=detector_stage,
            frozen_angles=frozen_angles,
            cut_points=cut_points,
        )
        assert mode.sample_stage == sample_stage
        assert mode.detector_stage == detector_stage
        assert mode.frozen_angles == (frozen_angles or {})
        assert mode.cut_points == (cut_points or {})
        # constrained_stages preserves insertion order
        assert mode.constrained_stages == expected_constrained


def test_bisecting_mode_repr():
    mode = BisectingMode(sample_stage="omega", detector_stage="ttheta")
    r = repr(mode)
    assert "BisectingMode" in r
    assert "omega" in r
    assert "ttheta" in r


def test_bisecting_mode_eq_same():
    m1 = BisectingMode(
        sample_stage="eta",
        detector_stage="delta",
        frozen_angles={"mu": 0.0, "nu": 0.0},
    )
    m2 = BisectingMode(
        sample_stage="eta",
        detector_stage="delta",
        frozen_angles={"mu": 0.0, "nu": 0.0},
    )
    assert m1 == m2


def test_bisecting_mode_eq_different_sample_stage():
    m1 = BisectingMode(sample_stage="eta", detector_stage="delta")
    m2 = BisectingMode(sample_stage="omega", detector_stage="delta")
    assert m1 != m2


def test_bisecting_mode_eq_different_detector_stage():
    m1 = BisectingMode(sample_stage="eta", detector_stage="delta")
    m2 = BisectingMode(sample_stage="eta", detector_stage="ttheta")
    assert m1 != m2


def test_bisecting_mode_eq_different_frozen():
    m1 = BisectingMode(
        sample_stage="eta", detector_stage="delta", frozen_angles={"mu": 0.0}
    )
    m2 = BisectingMode(sample_stage="eta", detector_stage="delta")
    assert m1 != m2


def test_bisecting_mode_eq_different_type():
    m1 = BisectingMode(sample_stage="omega", detector_stage="ttheta")
    m2 = FixedAngleMode(stage="omega", value=0.0)
    assert m1 != m2


def test_bisecting_mode_eq_non_mode():
    m = BisectingMode(sample_stage="omega", detector_stage="ttheta")
    assert m != 42


def test_bisecting_mode_constrained_stages_no_duplicate():
    """
    When a frozen_angles key equals the sample_stage name, it must NOT appear
    twice in constrained_stages — the 'if name not in constrained' guard works.
    """
    mode = BisectingMode(
        sample_stage="omega",
        detector_stage="ttheta",
        frozen_angles={"omega": 0.0, "chi": 90.0},  # omega also in frozen_angles
    )
    constrained = mode.constrained_stages
    # omega must appear only once (deduplicated by the guard)
    assert constrained.count("omega") == 1
    assert "chi" in constrained
    assert constrained[0] == "omega"  # sample_stage still first


# ---------------------------------------------------------------------------
# DiffractionMode.apply_cut_point
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cut_points, stage, angle, expected, context",
    [
        pytest.param(
            {"omega": -180.0},
            "omega",
            0.0,
            0.0 + 0.0,
            does_not_raise(),
            id="cut-at-neg180-angle-0",
        ),
        pytest.param(
            {"omega": -180.0},
            "omega",
            -180.0,
            -180.0,
            does_not_raise(),
            id="cut-at-neg180-angle-is-cut",
        ),
        pytest.param(
            {"omega": -180.0},
            "omega",
            185.0,
            -175.0,
            does_not_raise(),
            id="cut-at-neg180-wraps-185",
        ),
        pytest.param(
            {"omega": 0.0},
            "omega",
            350.0,
            350.0,
            does_not_raise(),
            id="cut-at-0-350-stays",
        ),
        pytest.param(
            {"omega": 0.0},
            "omega",
            -10.0,
            350.0,
            does_not_raise(),
            id="cut-at-0-neg10-wraps",
        ),
        pytest.param(
            {},
            "omega",
            350.0,
            350.0,
            does_not_raise(),
            id="no-cut-point-passes-through",
        ),
        pytest.param(
            {"chi": -180.0},
            "omega",
            350.0,
            350.0,
            does_not_raise(),
            id="cut-for-different-stage-passes-through",
        ),
    ],
)
def test_apply_cut_point(cut_points, stage, angle, expected, context):
    mode = FixedAngleMode(stage="phi", value=0.0, cut_points=cut_points)
    with context:
        result = mode.apply_cut_point(stage, angle)
        assert abs(result - expected) < 1e-10


# ---------------------------------------------------------------------------
# ModeDict
# ---------------------------------------------------------------------------


def test_mode_dict_construction_empty():
    md = ModeDict()
    assert len(md) == 0


def test_mode_dict_construction_from_dict():
    modes = {
        "bisecting": BisectingMode("omega", "ttheta"),
        "fixed_chi": FixedAngleMode("chi", 90.0),
    }
    md = ModeDict(modes)
    assert len(md) == 2
    assert "bisecting" in md
    assert "fixed_chi" in md


def test_mode_dict_setitem_valid():
    md = ModeDict()
    md["m1"] = FixedAngleMode("chi", 90.0)
    assert "m1" in md


def test_mode_dict_setitem_invalid_type():
    md = ModeDict()
    with pytest.raises(TypeError, match=re.escape("DiffractionMode instances")):
        md["bad"] = "not a mode"


def test_mode_dict_getitem():
    md = ModeDict({"m": FixedAngleMode("chi", 0.0)})
    assert isinstance(md["m"], FixedAngleMode)


def test_mode_dict_getitem_missing():
    md = ModeDict()
    with pytest.raises(KeyError):
        _ = md["missing"]


def test_mode_dict_delitem():
    md = ModeDict({"m": FixedAngleMode("chi", 0.0)})
    del md["m"]
    assert "m" not in md


def test_mode_dict_iter():
    keys = ["a", "b", "c"]
    md = ModeDict({k: FixedAngleMode("chi", 0.0) for k in keys})
    assert list(md) == keys


def test_mode_dict_keys():
    md = ModeDict({"a": FixedAngleMode("chi", 0.0), "b": FixedAngleMode("phi", 0.0)})
    assert set(md.keys()) == {"a", "b"}


def test_mode_dict_values():
    mode = FixedAngleMode("chi", 0.0)
    md = ModeDict({"m": mode})
    assert list(md.values()) == [mode]


def test_mode_dict_items():
    mode = FixedAngleMode("chi", 0.0)
    md = ModeDict({"m": mode})
    assert list(md.items()) == [("m", mode)]


def test_mode_dict_repr():
    md = ModeDict({"m": FixedAngleMode("chi", 0.0)})
    r = repr(md)
    assert "ModeDict" in r


def test_mode_dict_eq_same():
    m1 = ModeDict({"m": FixedAngleMode("chi", 0.0)})
    m2 = ModeDict({"m": FixedAngleMode("chi", 0.0)})
    assert m1 == m2


def test_mode_dict_eq_different():
    m1 = ModeDict({"m": FixedAngleMode("chi", 0.0)})
    m2 = ModeDict({"m": FixedAngleMode("chi", 90.0)})
    assert m1 != m2


def test_mode_dict_eq_non_mode_dict():
    md = ModeDict()
    assert md != {}


# ---------------------------------------------------------------------------
# AdHocDiffractometer — mode integration
# ---------------------------------------------------------------------------


def test_geometry_no_modes_by_default():
    """A geometry with no modes argument has an empty ModeDict."""
    g = _simple_geometry()
    assert isinstance(g.modes, ModeDict)
    assert len(g.modes) == 0
    assert g.mode_name is None
    assert g.mode is None


def test_geometry_modes_from_dict():
    modes = {
        "bisecting": BisectingMode("omega", "ttheta"),
        "fixed_chi": FixedAngleMode("chi", 90.0),
    }
    g = _simple_geometry(modes=modes, default_mode="bisecting")
    assert "bisecting" in g.modes
    assert "fixed_chi" in g.modes
    assert g.mode_name == "bisecting"
    assert isinstance(g.mode, BisectingMode)


def test_geometry_modes_from_mode_dict():
    md = ModeDict({"bisecting": BisectingMode("omega", "ttheta")})
    g = _simple_geometry(modes=md, default_mode="bisecting")
    assert g.mode_name == "bisecting"


def test_geometry_default_mode_none():
    modes = {"bisecting": BisectingMode("omega", "ttheta")}
    g = _simple_geometry(modes=modes)  # no default_mode
    assert g.mode_name is None
    assert g.mode is None


def test_geometry_default_mode_invalid():
    modes = {"bisecting": BisectingMode("omega", "ttheta")}
    with pytest.raises(ValueError, match=re.escape("default_mode")):
        _simple_geometry(modes=modes, default_mode="nonexistent")


def test_geometry_mode_name_setter_valid():
    modes = {
        "bisecting": BisectingMode("omega", "ttheta"),
        "fixed_chi": FixedAngleMode("chi", 90.0),
    }
    g = _simple_geometry(modes=modes, default_mode="bisecting")
    g.mode_name = "fixed_chi"
    assert g.mode_name == "fixed_chi"
    assert isinstance(g.mode, FixedAngleMode)


def test_geometry_mode_name_setter_none():
    modes = {"bisecting": BisectingMode("omega", "ttheta")}
    g = _simple_geometry(modes=modes, default_mode="bisecting")
    g.mode_name = None
    assert g.mode_name is None
    assert g.mode is None


def test_geometry_mode_name_setter_invalid():
    modes = {"bisecting": BisectingMode("omega", "ttheta")}
    g = _simple_geometry(modes=modes, default_mode="bisecting")
    with pytest.raises(ValueError, match=re.escape("not available")):
        g.mode_name = "nonexistent"


def test_geometry_cut_points_default_empty():
    g = _simple_geometry()
    assert g.cut_points == {}


def test_geometry_cut_points_set_at_construction():
    g = _simple_geometry(cut_points={"omega": -180.0, "ttheta": -180.0})
    assert g.cut_points == {"omega": -180.0, "ttheta": -180.0}


def test_geometry_cut_points_mutable():
    g = _simple_geometry()
    g.cut_points["omega"] = -170.0
    assert g.cut_points["omega"] == -170.0


# ---------------------------------------------------------------------------
# Factory canonical modes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, expected_modes, expected_default",
    [
        pytest.param(
            psic,
            {"bisecting", "fixed_chi", "fixed_phi", "fixed_mu"},
            "bisecting",
            id="psic-modes",
        ),
        pytest.param(
            fourcv,
            {"bisecting", "fixed_chi", "fixed_phi"},
            "bisecting",
            id="fourcv-modes",
        ),
        pytest.param(
            fourch,
            {"bisecting", "fixed_chi", "fixed_phi"},
            "bisecting",
            id="fourch-modes",
        ),
        pytest.param(
            kappa4cv,
            {"bisecting", "fixed_kphi"},
            "bisecting",
            id="kappa4cv-modes",
        ),
        pytest.param(
            kappa4ch,
            {"bisecting", "fixed_kphi"},
            "bisecting",
            id="kappa4ch-modes",
        ),
        pytest.param(
            kappa6c,
            {"bisecting", "fixed_kphi", "fixed_mu"},
            "bisecting",
            id="kappa6c-modes",
        ),
    ],
)
def test_factory_canonical_modes(factory, expected_modes, expected_default):
    g = factory()
    assert set(g.modes.keys()) == expected_modes
    assert g.mode_name == expected_default


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(sixc, id="sixc-no-modes"),
        pytest.param(zaxis, id="zaxis-no-modes"),
        pytest.param(s2d2, id="s2d2-no-modes"),
        pytest.param(fivec, id="fivec-no-modes"),
    ],
)
def test_factory_no_modes(factory):
    """Geometries without canonical modes have empty ModeDict and no active mode."""
    g = factory()
    assert len(g.modes) == 0
    assert g.mode_name is None


def test_psic_bisecting_mode_correct_stages():
    g = psic()
    mode = g.modes["bisecting"]
    assert isinstance(mode, BisectingMode)
    assert mode.sample_stage == "eta"
    assert mode.detector_stage == "delta"
    assert mode.frozen_angles == {"mu": 0.0, "nu": 0.0}


def test_fourcv_bisecting_mode_correct_stages():
    g = fourcv()
    mode = g.modes["bisecting"]
    assert isinstance(mode, BisectingMode)
    assert mode.sample_stage == "omega"
    assert mode.detector_stage == "ttheta"
    assert mode.frozen_angles == {}


def test_kappa6c_bisecting_mode_correct_stages():
    g = kappa6c()
    mode = g.modes["bisecting"]
    assert isinstance(mode, BisectingMode)
    assert mode.sample_stage == "komega"
    assert mode.detector_stage == "delta"
    assert mode.frozen_angles == {"mu": 0.0, "nu": 0.0}


# ---------------------------------------------------------------------------
# Serialisation — mode round-trip via to_dict / from_dict
# ---------------------------------------------------------------------------


def test_fixed_angle_mode_to_dict_from_dict():
    """FixedAngleMode round-trips through the geometry serialisation."""
    modes = {"fixed_chi": FixedAngleMode("chi", 90.0, cut_points={"chi": -180.0})}
    g = _simple_geometry(modes=modes, default_mode="fixed_chi")

    d = g.to_dict()
    assert "modes" in d
    assert "fixed_chi" in d["modes"]
    assert d["modes"]["fixed_chi"]["type"] == "FixedAngleMode"
    assert d["modes"]["fixed_chi"]["stage"] == "chi"
    assert d["modes"]["fixed_chi"]["value"] == 90.0
    assert d["modes"]["fixed_chi"]["cut_points"] == {"chi": -180.0}
    assert d["mode_name"] == "fixed_chi"

    # Verify JSON-serialisable
    assert json.dumps(d)

    g2 = AdHocDiffractometer.from_dict(d)
    assert "fixed_chi" in g2.modes
    assert g2.mode_name == "fixed_chi"
    mode2 = g2.modes["fixed_chi"]
    assert isinstance(mode2, FixedAngleMode)
    assert mode2.frozen_angles == {"chi": 90.0}
    assert mode2.cut_points == {"chi": -180.0}


def test_bisecting_mode_to_dict_from_dict():
    """BisectingMode round-trips through the geometry serialisation."""
    modes = {
        "bisecting": BisectingMode(
            "omega",
            "ttheta",
            frozen_angles={"chi": 0.0},
            cut_points={"omega": -180.0},
        )
    }
    g = _simple_geometry(modes=modes, default_mode="bisecting")

    d = g.to_dict()
    md = d["modes"]["bisecting"]
    assert md["type"] == "BisectingMode"
    assert md["sample_stage"] == "omega"
    assert md["detector_stage"] == "ttheta"
    assert md["frozen_angles"] == {"chi": 0.0}
    assert md["cut_points"] == {"omega": -180.0}

    # Verify JSON-serialisable
    assert json.dumps(d)

    g2 = AdHocDiffractometer.from_dict(d)
    mode2 = g2.modes["bisecting"]
    assert isinstance(mode2, BisectingMode)
    assert mode2.sample_stage == "omega"
    assert mode2.detector_stage == "ttheta"
    assert mode2.frozen_angles == {"chi": 0.0}
    assert mode2.cut_points == {"omega": -180.0}


def test_cut_points_round_trip():
    """Geometry-level cut_points survive to_dict / from_dict."""
    g = _simple_geometry(cut_points={"omega": -180.0, "ttheta": -170.0})
    d = g.to_dict()
    assert d["cut_points"] == {"omega": -180.0, "ttheta": -170.0}
    g2 = AdHocDiffractometer.from_dict(d)
    assert g2.cut_points == {"omega": -180.0, "ttheta": -170.0}


def test_mode_name_none_round_trip():
    """mode_name=None is stored and restored."""
    modes = {"bisecting": BisectingMode("omega", "ttheta")}
    g = _simple_geometry(modes=modes)  # no default_mode
    d = g.to_dict()
    assert d["mode_name"] is None
    g2 = AdHocDiffractometer.from_dict(d)
    assert g2.mode_name is None


def test_psic_full_round_trip():
    """psic() with all its modes survives to_dict / from_dict."""
    g = psic()
    d = g.to_dict()
    assert json.dumps(d)
    g2 = AdHocDiffractometer.from_dict(d)
    assert set(g2.modes.keys()) == {"bisecting", "fixed_chi", "fixed_phi", "fixed_mu"}
    assert g2.mode_name == "bisecting"
    assert isinstance(g2.modes["bisecting"], BisectingMode)
    assert isinstance(g2.modes["fixed_chi"], FixedAngleMode)


def test_geometry_no_modes_round_trip():
    """Geometry with no modes serialises with empty modes dict."""
    g = sixc()
    d = g.to_dict()
    assert d["modes"] == {}
    assert d["mode_name"] is None
    g2 = AdHocDiffractometer.from_dict(d)
    assert len(g2.modes) == 0
    assert g2.mode_name is None


# ---------------------------------------------------------------------------
# DiffractionMode base-class repr
# ---------------------------------------------------------------------------


def test_diffraction_mode_base_repr():
    """DiffractionMode.__repr__ includes class name and frozen/cut fields."""
    mode = FixedAngleMode("chi", 90.0, cut_points={"chi": -180.0})
    r = repr(mode)
    # FixedAngleMode overrides repr — check it is informative
    assert "chi" in r
    assert "90.0" in r
