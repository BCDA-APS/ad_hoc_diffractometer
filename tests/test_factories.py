"""
Unit tests for ad_hoc_diffractometer.factories.

Covers:
  - list_geometries()
  - get_geometry()
  - make_geometry()
  - All geometry factory functions: psic, fourc_v, fourc_h, sixc,
    kappa4c, kappa4c_h, kappa6c, zaxis, s2d2, fivec
"""

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

from ad_hoc_diffractometer import XHAT
from ad_hoc_diffractometer import YHAT
from ad_hoc_diffractometer import ZHAT
from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer import fivec
from ad_hoc_diffractometer import fourc_h
from ad_hoc_diffractometer import fourc_v
from ad_hoc_diffractometer import get_geometry
from ad_hoc_diffractometer import kappa4c
from ad_hoc_diffractometer import kappa4c_h
from ad_hoc_diffractometer import kappa6c
from ad_hoc_diffractometer import list_geometries
from ad_hoc_diffractometer import make_geometry
from ad_hoc_diffractometer import psic
from ad_hoc_diffractometer import s2d2
from ad_hoc_diffractometer import sixc
from ad_hoc_diffractometer import zaxis

# ---------------------------------------------------------------------------
# list_geometries()
# ---------------------------------------------------------------------------


def test_list_geometries_returns_all_factories():
    geoms = list_geometries()
    expected = {
        "psic",
        "fourc_v",
        "fourc_h",
        "sixc",
        "kappa4c",
        "kappa4c_h",
        "kappa6c",
        "zaxis",
        "s2d2",
        "fivec",
    }
    assert set(geoms.keys()) == expected


def test_list_geometries_values_are_callable():
    for name, func in list_geometries().items():
        assert callable(func), f"{name} is not callable"


def test_list_geometries_instantiate_all():
    """Every registered factory must instantiate without error."""
    for name, func in list_geometries().items():
        g = func()
        assert g.name == name, f"Factory {name!r} returned name {g.name!r}"


def test_list_geometries_returns_copy():
    """Mutating the returned dict must not affect the registry."""
    geoms = list_geometries()
    geoms.clear()
    assert len(list_geometries()) > 0


# ---------------------------------------------------------------------------
# get_geometry()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, context",
    [
        pytest.param("psic", does_not_raise(), id="get-psic"),
        pytest.param("fourc_v", does_not_raise(), id="get-fourc_v"),
        pytest.param("fourc_h", does_not_raise(), id="get-fourc_h"),
        pytest.param("sixc", does_not_raise(), id="get-sixc"),
        pytest.param("kappa4c", does_not_raise(), id="get-kappa4c"),
        pytest.param("kappa4c_h", does_not_raise(), id="get-kappa4c_h"),
        pytest.param("kappa6c", does_not_raise(), id="get-kappa6c"),
        pytest.param("zaxis", does_not_raise(), id="get-zaxis"),
        pytest.param("s2d2", does_not_raise(), id="get-s2d2"),
        pytest.param("fivec", does_not_raise(), id="get-fivec"),
        pytest.param(
            "bogus",
            pytest.raises(ValueError, match=re.escape("No geometry named")),
            id="get-invalid-name-raises",
        ),
        pytest.param(
            "",
            pytest.raises(ValueError, match=re.escape("No geometry named")),
            id="get-empty-string-raises",
        ),
    ],
)
def test_get_geometry(name, context):
    with context:
        factory = get_geometry(name)
        assert callable(factory)
        assert factory.__name__ == name


def test_get_geometry_error_lists_available():
    """Error message for unknown name must list available geometry names."""
    with pytest.raises(ValueError, match=re.escape("Available geometries")):
        get_geometry("bogus")


# ---------------------------------------------------------------------------
# make_geometry()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, kwargs, context",
    [
        pytest.param("psic", {}, does_not_raise(), id="make-psic-no-kwargs"),
        pytest.param("fourc_v", {}, does_not_raise(), id="make-fourc_v-no-kwargs"),
        pytest.param("kappa4c", {}, does_not_raise(), id="make-kappa4c-default-alpha"),
        pytest.param(
            "kappa4c",
            {"alpha_deg": 50.0},
            does_not_raise(),
            id="make-kappa4c-explicit-alpha-50",
        ),
        pytest.param(
            "kappa6c",
            {"alpha_deg": 55.0},
            does_not_raise(),
            id="make-kappa6c-explicit-alpha-55",
        ),
        pytest.param(
            "bogus",
            {},
            pytest.raises(ValueError, match=re.escape("No geometry named")),
            id="make-invalid-name-raises",
        ),
    ],
)
def test_make_geometry(name, kwargs, context):
    with context:
        g = make_geometry(name, **kwargs)
        assert isinstance(g, AdHocDiffractometer)
        assert g.name == name


def test_make_geometry_returns_correct_instance():
    """make_geometry and the factory function must return equivalent instances."""
    g1 = make_geometry("psic")
    g2 = psic()
    assert g1.name == g2.name
    assert [s.name for s in g1.sample_stages] == [s.name for s in g2.sample_stages]
    assert [s.name for s in g1.detector_stages] == [s.name for s in g2.detector_stages]


def test_make_geometry_kappa_alpha_forwarded():
    """Keyword args must be forwarded to the factory (kappa alpha test)."""
    g = make_geometry("kappa4c", alpha_deg=45.0)
    expected = np.cos(np.deg2rad(45)) * np.array([0, 0, 1]) + np.sin(
        np.deg2rad(45)
    ) * np.array([1, 0, 0])
    np.testing.assert_allclose(g.stage("kappa").axis, expected, atol=1e-12)


# ---------------------------------------------------------------------------
# Factory stage lists
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, expected_name, sample_names, detector_names, context",
    [
        pytest.param(
            psic,
            "psic",
            ["mu", "eta", "chi", "phi"],
            ["nu", "delta"],
            does_not_raise(),
            id="psic-stage-lists",
        ),
        pytest.param(
            fourc_v,
            "fourc_v",
            ["omega", "chi", "phi"],
            ["two_theta"],
            does_not_raise(),
            id="fourc_v-stage-lists",
        ),
        pytest.param(
            fourc_h,
            "fourc_h",
            ["omega", "chi", "phi"],
            ["two_theta"],
            does_not_raise(),
            id="fourc_h-stage-lists",
        ),
        pytest.param(
            sixc,
            "sixc",
            ["alpha", "omega", "chi", "phi"],
            ["delta", "gamma"],
            does_not_raise(),
            id="sixc-stage-lists",
        ),
        pytest.param(
            kappa4c,
            "kappa4c",
            ["komega", "kappa", "kphi"],
            ["two_theta"],
            does_not_raise(),
            id="kappa4c-stage-lists",
        ),
        pytest.param(
            kappa4c_h,
            "kappa4c_h",
            ["komega", "kappa", "kphi"],
            ["two_theta"],
            does_not_raise(),
            id="kappa4c_h-stage-lists",
        ),
        pytest.param(
            kappa6c,
            "kappa6c",
            ["mu", "komega", "kappa", "kphi"],
            ["nu", "delta"],
            does_not_raise(),
            id="kappa6c-stage-lists",
        ),
        pytest.param(
            zaxis,
            "zaxis",
            ["alpha", "Z"],
            ["delta", "gamma"],
            does_not_raise(),
            id="zaxis-stage-lists",
        ),
        pytest.param(
            s2d2,
            "s2d2",
            ["mu", "Z"],
            ["nu", "delta"],
            does_not_raise(),
            id="s2d2-stage-lists",
        ),
        pytest.param(
            fivec,
            "fivec",
            ["mu", "omega", "chi", "phi"],
            ["two_theta"],
            does_not_raise(),
            id="fivec-stage-lists",
        ),
    ],
)
def test_geometry_factories(
    factory, expected_name, sample_names, detector_names, context
):
    with context:
        g = factory()
        assert g.name == expected_name
        assert [s.name for s in g.sample_stages] == sample_names
        assert [s.name for s in g.detector_stages] == detector_names


# ---------------------------------------------------------------------------
# Factory parent chains
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, stage_name, expected_parent, context",
    [
        # psic
        pytest.param(psic, "mu", None, does_not_raise(), id="psic-mu-on-floor"),
        pytest.param(psic, "nu", None, does_not_raise(), id="psic-nu-on-floor"),
        pytest.param(psic, "eta", "mu", does_not_raise(), id="psic-eta-on-mu"),
        pytest.param(psic, "chi", "eta", does_not_raise(), id="psic-chi-on-eta"),
        pytest.param(psic, "phi", "chi", does_not_raise(), id="psic-phi-on-chi"),
        pytest.param(psic, "delta", "nu", does_not_raise(), id="psic-delta-on-nu"),
        # fourc_v
        pytest.param(
            fourc_v, "omega", None, does_not_raise(), id="fourc_v-omega-on-floor"
        ),
        pytest.param(
            fourc_v,
            "two_theta",
            None,
            does_not_raise(),
            id="fourc_v-two_theta-decoupled",
        ),
        pytest.param(
            fourc_v, "chi", "omega", does_not_raise(), id="fourc_v-chi-on-omega"
        ),
        pytest.param(fourc_v, "phi", "chi", does_not_raise(), id="fourc_v-phi-on-chi"),
        # fourc_h
        pytest.param(
            fourc_h,
            "two_theta",
            None,
            does_not_raise(),
            id="fourc_h-two_theta-decoupled",
        ),
        # sixc
        pytest.param(sixc, "alpha", None, does_not_raise(), id="sixc-alpha-on-floor"),
        pytest.param(
            sixc, "omega", "alpha", does_not_raise(), id="sixc-omega-on-alpha"
        ),
        pytest.param(
            sixc, "delta", "alpha", does_not_raise(), id="sixc-delta-on-alpha"
        ),
        pytest.param(
            sixc, "gamma", "delta", does_not_raise(), id="sixc-gamma-on-delta"
        ),
        # kappa4c
        pytest.param(
            kappa4c, "komega", None, does_not_raise(), id="kappa4c-komega-on-floor"
        ),
        pytest.param(
            kappa4c,
            "two_theta",
            None,
            does_not_raise(),
            id="kappa4c-two_theta-decoupled",
        ),
        pytest.param(
            kappa4c, "kappa", "komega", does_not_raise(), id="kappa4c-kappa-on-komega"
        ),
        pytest.param(
            kappa4c, "kphi", "kappa", does_not_raise(), id="kappa4c-kphi-on-kappa"
        ),
        # kappa6c
        pytest.param(kappa6c, "mu", None, does_not_raise(), id="kappa6c-mu-on-floor"),
        pytest.param(kappa6c, "nu", None, does_not_raise(), id="kappa6c-nu-on-floor"),
        pytest.param(
            kappa6c, "komega", "mu", does_not_raise(), id="kappa6c-komega-on-mu"
        ),
        pytest.param(
            kappa6c, "kappa", "komega", does_not_raise(), id="kappa6c-kappa-on-komega"
        ),
        pytest.param(
            kappa6c, "kphi", "kappa", does_not_raise(), id="kappa6c-kphi-on-kappa"
        ),
        # zaxis
        pytest.param(zaxis, "alpha", None, does_not_raise(), id="zaxis-alpha-on-floor"),
        pytest.param(zaxis, "Z", "alpha", does_not_raise(), id="zaxis-Z-on-alpha"),
        pytest.param(
            zaxis, "delta", "alpha", does_not_raise(), id="zaxis-delta-on-alpha"
        ),
        pytest.param(
            zaxis, "gamma", "delta", does_not_raise(), id="zaxis-gamma-on-delta"
        ),
        # s2d2
        pytest.param(s2d2, "mu", None, does_not_raise(), id="s2d2-mu-on-floor"),
        pytest.param(s2d2, "nu", None, does_not_raise(), id="s2d2-nu-on-floor"),
        pytest.param(s2d2, "Z", "mu", does_not_raise(), id="s2d2-Z-on-mu"),
        pytest.param(s2d2, "delta", "nu", does_not_raise(), id="s2d2-delta-on-nu"),
        # fivec
        pytest.param(fivec, "mu", None, does_not_raise(), id="fivec-mu-on-floor"),
        pytest.param(fivec, "omega", "mu", does_not_raise(), id="fivec-omega-on-mu"),
        pytest.param(
            fivec, "two_theta", "mu", does_not_raise(), id="fivec-two_theta-on-mu"
        ),
        pytest.param(fivec, "chi", "omega", does_not_raise(), id="fivec-chi-on-omega"),
        pytest.param(fivec, "phi", "chi", does_not_raise(), id="fivec-phi-on-chi"),
    ],
)
def test_geometry_parent_chain(factory, stage_name, expected_parent, context):
    with context:
        g = factory()
        assert g.stage(stage_name).parent == expected_parent


# ---------------------------------------------------------------------------
# Factory axis vectors
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, stage_name, expected_axis, context",
    [
        pytest.param(
            psic, "mu", +XHAT, does_not_raise(), id="psic-mu-vertical-right-handed"
        ),
        pytest.param(
            psic, "eta", -ZHAT, does_not_raise(), id="psic-eta-lateral-left-handed"
        ),
        pytest.param(
            psic,
            "chi",
            +YHAT,
            does_not_raise(),
            id="psic-chi-longitudinal-right-handed",
        ),
        pytest.param(
            psic, "phi", -ZHAT, does_not_raise(), id="psic-phi-lateral-left-handed"
        ),
        pytest.param(
            psic, "nu", +XHAT, does_not_raise(), id="psic-nu-vertical-right-handed"
        ),
        pytest.param(
            psic, "delta", -ZHAT, does_not_raise(), id="psic-delta-lateral-left-handed"
        ),
        pytest.param(
            kappa4c,
            "kappa",
            np.cos(np.deg2rad(50)) * np.array([0, 0, 1])
            + np.sin(np.deg2rad(50)) * np.array([1, 0, 0]),
            does_not_raise(),
            id="kappa4c-kappa-axis-tilted",
        ),
    ],
)
def test_geometry_axes(factory, stage_name, expected_axis, context):
    with context:
        g = factory()
        np.testing.assert_allclose(g.stage(stage_name).axis, expected_axis, atol=1e-12)


# ---------------------------------------------------------------------------
# kappa_alpha_deg on kappa factory instances
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, alpha_deg, context",
    [
        pytest.param(kappa4c, 50.0, does_not_raise(), id="kappa4c-default-alpha"),
        pytest.param(kappa4c_h, 50.0, does_not_raise(), id="kappa4c_h-default-alpha"),
        pytest.param(kappa6c, 50.0, does_not_raise(), id="kappa6c-default-alpha"),
        pytest.param(kappa4c, 45.0, does_not_raise(), id="kappa4c-custom-alpha-45"),
        pytest.param(kappa4c_h, 35.0, does_not_raise(), id="kappa4c_h-custom-alpha-35"),
        pytest.param(kappa6c, 55.0, does_not_raise(), id="kappa6c-custom-alpha-55"),
    ],
)
def test_kappa_alpha_deg_stored(factory, alpha_deg, context):
    """kappa_alpha_deg on the instance matches the alpha_deg passed to the factory."""
    with context:
        g = factory(alpha_deg=alpha_deg)
        assert g.kappa_alpha_deg == pytest.approx(alpha_deg)


@pytest.mark.parametrize(
    "factory, alpha_deg, context",
    [
        pytest.param(kappa4c, 50.0, does_not_raise(), id="kappa4c-axis-matches-50"),
        pytest.param(kappa4c, 45.0, does_not_raise(), id="kappa4c-axis-matches-45"),
        pytest.param(kappa6c, 55.0, does_not_raise(), id="kappa6c-axis-matches-55"),
    ],
)
def test_kappa_alpha_deg_matches_axis_vector(factory, alpha_deg, context):
    """kappa_alpha_deg is consistent with the kappa stage axis vector."""
    from ad_hoc_diffractometer import kappa_axis
    from ad_hoc_diffractometer.factories import _BASIS_BL
    from ad_hoc_diffractometer.factories import _BASIS_YOU

    with context:
        g = factory(alpha_deg=alpha_deg)
        kax = g.stage("kappa").axis
        basis = _BASIS_YOU if factory is kappa6c else _BASIS_BL
        expected_axis = kappa_axis(g.kappa_alpha_deg, basis=basis)
        np.testing.assert_allclose(kax, expected_axis, atol=1e-12)
