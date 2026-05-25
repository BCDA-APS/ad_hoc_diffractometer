# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
Unit tests for ad_hoc_diffractometer.factories.

Covers:
  - list_geometries()
  - get_geometry()
  - make_geometry()
  - All geometry factory functions: psic, fourcv, fourch, sixc,
    kappa4cv, kappa4ch, kappa6c, zaxis, s2d2, fivec
  - Entry-point extensibility (#37):
    - GEOMETRY_ENTRY_POINT_GROUP constant
    - All 10 built-in factories declared as entry points in pyproject.toml
      and discoverable via importlib.metadata
    - list_geometries() / get_geometry() load plugins via entry points
    - Mock-based test: a simulated third-party plugin is picked up
    - Broken plugin silently skipped
"""

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from helpers import fivec
from helpers import fourch
from helpers import fourcv
from helpers import kappa4ch
from helpers import kappa4cv
from helpers import kappa6c
from helpers import psic
from helpers import s2d2
from helpers import sixc
from helpers import zaxis

from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer import get_geometry
from ad_hoc_diffractometer import list_geometries
from ad_hoc_diffractometer import make_geometry
from ad_hoc_diffractometer.constants import XHAT
from ad_hoc_diffractometer.constants import YHAT
from ad_hoc_diffractometer.constants import ZHAT

# ---------------------------------------------------------------------------
# list_geometries()
# ---------------------------------------------------------------------------


def test_list_geometries_returns_all_factories():
    geoms = list_geometries()
    expected = {
        "psic",
        "fourcv",
        "fourch",
        "sixc",
        "kappa4cv",
        "kappa4ch",
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
        pytest.param("fourcv", does_not_raise(), id="get-fourcv"),
        pytest.param("fourch", does_not_raise(), id="get-fourch"),
        pytest.param("sixc", does_not_raise(), id="get-sixc"),
        pytest.param("kappa4cv", does_not_raise(), id="get-kappa4cv"),
        pytest.param("kappa4ch", does_not_raise(), id="get-kappa4ch"),
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
        pytest.param("fourcv", {}, does_not_raise(), id="make-fourcv-no-kwargs"),
        pytest.param(
            "kappa4cv", {}, does_not_raise(), id="make-kappa4cv-default-alpha"
        ),
        pytest.param(
            "kappa4cv",
            {"alpha_deg": 50.0},
            does_not_raise(),
            id="make-kappa4cv-explicit-alpha-50",
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
    """Keyword args must be forwarded to the factory (kappa alpha test).

    The kappa axis lies in the transverse-vertical plane, tilted
    ``alpha`` degrees from +T toward +V (Walko 2016 Fig. 3 and
    Thorkildsen et al. 2006 Table 1; see issue #252 for the
    correction of the kappa-arm tilt direction).  For ``kappa4cv``
    (BL basis: T=+x, L=+y, V=+z) the kappa axis is therefore
    ``+x̂·cos(α) + ẑ·sin(α)``.

    Issue #284 separated the kappa-arm tilt direction (which the YAML
    ``kappa_chi_eq`` field still controls; tested here) from the
    equivalent-Eulerian chi pseudo-angle axis (controlled by the
    new YAML ``kappa_eulerian_chi`` field; auto-derived to
    ``+longitudinal`` for every shipped kappa preset).
    """
    g = make_geometry("kappa4cv", alpha_deg=45.0)
    expected = np.cos(np.deg2rad(45)) * np.array([1, 0, 0]) + np.sin(
        np.deg2rad(45)
    ) * np.array([0, 0, 1])
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
            fourcv,
            "fourcv",
            ["omega", "chi", "phi"],
            ["ttheta"],
            does_not_raise(),
            id="fourcv-stage-lists",
        ),
        pytest.param(
            fourch,
            "fourch",
            ["omega", "chi", "phi"],
            ["ttheta"],
            does_not_raise(),
            id="fourch-stage-lists",
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
            kappa4cv,
            "kappa4cv",
            ["komega", "kappa", "kphi"],
            ["ttheta"],
            does_not_raise(),
            id="kappa4cv-stage-lists",
        ),
        pytest.param(
            kappa4ch,
            "kappa4ch",
            ["komega", "kappa", "kphi"],
            ["ttheta"],
            does_not_raise(),
            id="kappa4ch-stage-lists",
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
            ["ttheta"],
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
        # fourcv
        pytest.param(
            fourcv, "omega", None, does_not_raise(), id="fourcv-omega-on-floor"
        ),
        pytest.param(
            fourcv,
            "ttheta",
            None,
            does_not_raise(),
            id="fourcv-ttheta-decoupled",
        ),
        pytest.param(
            fourcv, "chi", "omega", does_not_raise(), id="fourcv-chi-on-omega"
        ),
        pytest.param(fourcv, "phi", "chi", does_not_raise(), id="fourcv-phi-on-chi"),
        # fourch
        pytest.param(
            fourch,
            "ttheta",
            None,
            does_not_raise(),
            id="fourch-ttheta-decoupled",
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
        # kappa4cv
        pytest.param(
            kappa4cv, "komega", None, does_not_raise(), id="kappa4cv-komega-on-floor"
        ),
        pytest.param(
            kappa4cv,
            "ttheta",
            None,
            does_not_raise(),
            id="kappa4cv-ttheta-decoupled",
        ),
        pytest.param(
            kappa4cv, "kappa", "komega", does_not_raise(), id="kappa4cv-kappa-on-komega"
        ),
        pytest.param(
            kappa4cv, "kphi", "kappa", does_not_raise(), id="kappa4cv-kphi-on-kappa"
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
        pytest.param(fivec, "ttheta", "mu", does_not_raise(), id="fivec-ttheta-on-mu"),
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
            psic, "eta", -ZHAT, does_not_raise(), id="psic-eta-transverse-left-handed"
        ),
        pytest.param(
            psic,
            "chi",
            +YHAT,
            does_not_raise(),
            id="psic-chi-longitudinal-right-handed",
        ),
        pytest.param(
            psic, "phi", -ZHAT, does_not_raise(), id="psic-phi-transverse-left-handed"
        ),
        pytest.param(
            psic, "nu", +XHAT, does_not_raise(), id="psic-nu-vertical-right-handed"
        ),
        pytest.param(
            psic,
            "delta",
            -ZHAT,
            does_not_raise(),
            id="psic-delta-transverse-left-handed",
        ),
        # fourcv — BL1967, vertical scattering plane; basis: transverse=+x, long=+y, vert=+z
        pytest.param(
            fourcv,
            "omega",
            -XHAT,
            does_not_raise(),
            id="fourcv-omega-transverse-left-handed",
        ),
        pytest.param(
            fourcv,
            "chi",
            +YHAT,
            does_not_raise(),
            id="fourcv-chi-longitudinal-right-handed",
        ),
        pytest.param(
            fourcv,
            "phi",
            -XHAT,
            does_not_raise(),
            id="fourcv-phi-transverse-left-handed",
        ),
        pytest.param(
            fourcv,
            "ttheta",
            -XHAT,
            does_not_raise(),
            id="fourcv-ttheta-transverse-left-handed",
        ),
        # fourch — BL1967, horizontal scattering plane; basis: transverse=+x, long=+y, vert=+z
        pytest.param(
            fourch,
            "omega",
            -ZHAT,
            does_not_raise(),
            id="fourch-omega-vertical-left-handed",
        ),
        pytest.param(
            fourch,
            "chi",
            +YHAT,
            does_not_raise(),
            id="fourch-chi-longitudinal-right-handed",
        ),
        pytest.param(
            fourch, "phi", -ZHAT, does_not_raise(), id="fourch-phi-vertical-left-handed"
        ),
        pytest.param(
            fourch,
            "ttheta",
            -ZHAT,
            does_not_raise(),
            id="fourch-ttheta-vertical-left-handed",
        ),
        # sixc — LV1993 Fig. 1 and §2.1:
        #   alpha, gamma: vertical (+x), right-handed
        #   omega, phi, delta: transverse (-z), left-handed
        #   chi: longitudinal (+y), right-handed
        pytest.param(
            sixc,
            "alpha",
            +XHAT,
            does_not_raise(),
            id="sixc-alpha-vertical-right-handed",
        ),
        pytest.param(
            sixc,
            "omega",
            -ZHAT,
            does_not_raise(),
            id="sixc-omega-transverse-left-handed",
        ),
        pytest.param(
            sixc,
            "chi",
            +YHAT,
            does_not_raise(),
            id="sixc-chi-longitudinal-right-handed",
        ),
        pytest.param(
            sixc, "phi", -ZHAT, does_not_raise(), id="sixc-phi-transverse-left-handed"
        ),
        pytest.param(
            sixc,
            "delta",
            -ZHAT,
            does_not_raise(),
            id="sixc-delta-transverse-left-handed",
        ),
        pytest.param(
            sixc,
            "gamma",
            +XHAT,
            does_not_raise(),
            id="sixc-gamma-vertical-right-handed",
        ),
        # kappa4cv — BL1967 basis: transverse=+x, long=+y, vert=+z; vertical scattering plane
        pytest.param(
            kappa4cv,
            "komega",
            -XHAT,
            does_not_raise(),
            id="kappa4cv-komega-transverse-left-handed",
        ),
        pytest.param(
            kappa4cv,
            "kphi",
            -XHAT,
            does_not_raise(),
            id="kappa4cv-kphi-transverse-left-handed",
        ),
        pytest.param(
            kappa4cv,
            "ttheta",
            -XHAT,
            does_not_raise(),
            id="kappa4cv-ttheta-transverse-left-handed",
        ),
        pytest.param(
            kappa4cv,
            "kappa",
            # The kappa axis lies in the transverse-vertical plane,
            # tilted α from +T toward +V (Walko 2016 Fig. 3 and
            # Thorkildsen 2006 Table 1; see issue #252).  For kappa4cv
            # (BL: T=+x, L=+y, V=+z): n_kappa = +x̂·cos α + ẑ·sin α.
            np.cos(np.deg2rad(50)) * np.array([1, 0, 0])
            + np.sin(np.deg2rad(50)) * np.array([0, 0, 1]),
            does_not_raise(),
            id="kappa4cv-kappa-axis-tilted",
        ),
        # kappa4ch — BL1967 basis; horizontal scattering plane
        pytest.param(
            kappa4ch,
            "komega",
            -ZHAT,
            does_not_raise(),
            id="kappa4ch-komega-vertical-left-handed",
        ),
        pytest.param(
            kappa4ch,
            "kphi",
            -ZHAT,
            does_not_raise(),
            id="kappa4ch-kphi-vertical-left-handed",
        ),
        pytest.param(
            kappa4ch,
            "ttheta",
            -ZHAT,
            does_not_raise(),
            id="kappa4ch-ttheta-vertical-left-handed",
        ),
        # kappa6c — You basis: vert=+x, long=+y, lat=+z; psic-style outer axes
        pytest.param(
            kappa6c,
            "mu",
            +XHAT,
            does_not_raise(),
            id="kappa6c-mu-vertical-right-handed",
        ),
        pytest.param(
            kappa6c,
            "komega",
            -ZHAT,
            does_not_raise(),
            id="kappa6c-komega-transverse-left-handed",
        ),
        pytest.param(
            kappa6c,
            "kphi",
            -ZHAT,
            does_not_raise(),
            id="kappa6c-kphi-transverse-left-handed",
        ),
        pytest.param(
            kappa6c,
            "nu",
            +XHAT,
            does_not_raise(),
            id="kappa6c-nu-vertical-right-handed",
        ),
        pytest.param(
            kappa6c,
            "delta",
            -ZHAT,
            does_not_raise(),
            id="kappa6c-delta-transverse-left-handed",
        ),
        # zaxis — You basis: vert=+x, long=+y, lat=+z
        pytest.param(
            zaxis,
            "alpha",
            +XHAT,
            does_not_raise(),
            id="zaxis-alpha-vertical-right-handed",
        ),
        pytest.param(
            zaxis,
            "Z",
            +YHAT,
            does_not_raise(),
            id="zaxis-Z-longitudinal-right-handed",
        ),
        pytest.param(
            zaxis,
            "delta",
            -ZHAT,
            does_not_raise(),
            id="zaxis-delta-transverse-left-handed",
        ),
        pytest.param(
            zaxis,
            "gamma",
            +XHAT,
            does_not_raise(),
            id="zaxis-gamma-vertical-right-handed",
        ),
        # s2d2 — You basis: vert=+x, long=+y, lat=+z
        pytest.param(
            s2d2, "mu", +XHAT, does_not_raise(), id="s2d2-mu-vertical-right-handed"
        ),
        pytest.param(
            s2d2, "Z", +YHAT, does_not_raise(), id="s2d2-Z-longitudinal-right-handed"
        ),
        pytest.param(
            s2d2, "nu", +XHAT, does_not_raise(), id="s2d2-nu-vertical-right-handed"
        ),
        pytest.param(
            s2d2,
            "delta",
            -ZHAT,
            does_not_raise(),
            id="s2d2-delta-transverse-left-handed",
        ),
        # fivec — You basis: vert=+x, long=+y, lat=+z; fourcv on vertical base
        pytest.param(
            fivec, "mu", +XHAT, does_not_raise(), id="fivec-mu-vertical-right-handed"
        ),
        pytest.param(
            fivec,
            "omega",
            -ZHAT,
            does_not_raise(),
            id="fivec-omega-transverse-left-handed",
        ),
        pytest.param(
            fivec,
            "chi",
            +YHAT,
            does_not_raise(),
            id="fivec-chi-longitudinal-right-handed",
        ),
        pytest.param(
            fivec, "phi", -ZHAT, does_not_raise(), id="fivec-phi-transverse-left-handed"
        ),
        pytest.param(
            fivec,
            "ttheta",
            -ZHAT,
            does_not_raise(),
            id="fivec-ttheta-transverse-left-handed",
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
        pytest.param(kappa4cv, 50.0, does_not_raise(), id="kappa4cv-default-alpha"),
        pytest.param(kappa4ch, 50.0, does_not_raise(), id="kappa4ch-default-alpha"),
        pytest.param(kappa6c, 50.0, does_not_raise(), id="kappa6c-default-alpha"),
        pytest.param(kappa4cv, 45.0, does_not_raise(), id="kappa4cv-custom-alpha-45"),
        pytest.param(kappa4ch, 35.0, does_not_raise(), id="kappa4ch-custom-alpha-35"),
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
        pytest.param(kappa4cv, 50.0, does_not_raise(), id="kappa4cv-axis-matches-50"),
        pytest.param(kappa4cv, 45.0, does_not_raise(), id="kappa4cv-axis-matches-45"),
        pytest.param(kappa6c, 55.0, does_not_raise(), id="kappa6c-axis-matches-55"),
    ],
)
def test_kappa_alpha_deg_matches_axis_vector(factory, alpha_deg, context):
    """kappa_alpha_deg is consistent with the kappa stage axis vector.

    Per Walko (2016) Fig. 3, Wyckoff (1985) Fig. 2(b), and Thorkildsen
    et al. (2006) Table 1, the kappa axis lies in the plane spanned
    by two basis directions and is tilted ``alpha`` degrees from one
    toward the other (see issue #252):

        kappa4cv: in T–V plane, between +T and +V (α from +T toward +V)
        kappa4ch: in V–L plane, between +V and +L (α from +V toward +L)
        kappa6c:  same as kappa4cv

    The angle between the kappa stage axis and each of the two
    spanning basis directions must be ``alpha`` and ``90°−alpha``
    respectively.
    """
    with context:
        g = factory(alpha_deg=alpha_deg)
        kax = np.asarray(g.stage("kappa").axis, dtype=float)
        kax = kax / np.linalg.norm(kax)

        plane_directions = {
            "kappa4cv": ("transverse", "vertical"),
            "kappa4ch": ("vertical", "longitudinal"),
            "kappa6c": ("transverse", "vertical"),
        }
        from_name, to_name = plane_directions[g.name]
        from_vec = np.asarray(g.basis[from_name], dtype=float)
        to_vec = np.asarray(g.basis[to_name], dtype=float)

        cos_a = np.cos(np.deg2rad(alpha_deg))
        sin_a = np.sin(np.deg2rad(alpha_deg))
        expected_axis = cos_a * from_vec + sin_a * to_vec
        np.testing.assert_allclose(kax, expected_axis, atol=1e-12)


# ---------------------------------------------------------------------------
# Entry-point extensibility (#37)
# ---------------------------------------------------------------------------

_BUILTIN_NAMES = frozenset(
    [
        "psic",
        "fourcv",
        "fourch",
        "sixc",
        "kappa4cv",
        "kappa4ch",
        "kappa6c",
        "zaxis",
        "s2d2",
        "fivec",
    ]
)


class TestEntryPointExtensibility:
    """Tests for the entry-point plugin mechanism (issue #37)."""

    # --- GEOMETRY_ENTRY_POINT_GROUP constant --------------------------------

    def test_entry_point_group_constant_value(self):
        """GEOMETRY_ENTRY_POINT_GROUP must be the expected string."""
        from ad_hoc_diffractometer.factories import GEOMETRY_ENTRY_POINT_GROUP

        assert GEOMETRY_ENTRY_POINT_GROUP == "ad_hoc_diffractometer.geometries"

    def test_entry_point_group_exported(self):
        """GEOMETRY_ENTRY_POINT_GROUP must be accessible from ad_hoc_diffractometer.factories."""
        import ad_hoc_diffractometer.factories as fac

        assert hasattr(fac, "GEOMETRY_ENTRY_POINT_GROUP")

    # --- Built-in factories declared via the declarative-YAML loader -----
    #
    # Issue #267 removed the entry-point declarations for the 10 demo
    # geometries from ``pyproject.toml`` and replaced the legacy
    # ``ad_hoc_diffractometer.presets`` Python factories with declarative
    # YAML files.  Built-ins are now registered exclusively by the
    # loader scanning :mod:`ad_hoc_diffractometer.geometries`.  The
    # entry-point group remains supported for *third-party* plugins.

    def test_no_builtins_declared_as_entry_points(self):
        """Built-in geometries are no longer declared as entry points."""
        from importlib.metadata import entry_points

        from ad_hoc_diffractometer.factories import GEOMETRY_ENTRY_POINT_GROUP

        eps = entry_points(group=GEOMETRY_ENTRY_POINT_GROUP)
        names = {ep.name for ep in eps}
        # Any of the 10 built-in names appearing here would mean a stale
        # ``pyproject.toml`` declaration leaked through; a clean install
        # produces an empty set.
        assert _BUILTIN_NAMES.isdisjoint(names), (
            f"Built-in names should not appear as entry points: "
            f"{sorted(_BUILTIN_NAMES & names)}"
        )

    # --- list_geometries() discovers entry points --------------------------

    def test_list_geometries_contains_all_builtins(self):
        """list_geometries() must return all 10 built-in geometry names."""
        geoms = list_geometries()
        assert _BUILTIN_NAMES <= set(geoms.keys()), (
            f"Missing: {_BUILTIN_NAMES - set(geoms.keys())}"
        )

    def test_list_geometries_picks_up_plugin_via_mock(self):
        """
        A simulated third-party plugin installed as an entry point must appear
        in list_geometries() after the registry is reset.
        """
        from unittest.mock import MagicMock
        from unittest.mock import patch

        from helpers import fourcv

        import ad_hoc_diffractometer.factories as fac

        # Build a fake entry point that loads a trivial factory
        fake_factory = fourcv  # reuse an existing factory as the "plugin"
        fake_ep = MagicMock()
        fake_ep.name = "my_plugin_geom"
        fake_ep.load.return_value = fake_factory

        # Reset the EP-loaded flag so discovery runs again
        original_ep_loaded = fac._EP_LOADED
        original_registry = dict(fac._GEOMETRY_REGISTRY)
        fac._EP_LOADED = False

        try:
            with patch(
                "ad_hoc_diffractometer.factories.entry_points",
                return_value=[fake_ep],
            ):
                geoms = list_geometries()
            assert "my_plugin_geom" in geoms
            assert geoms["my_plugin_geom"] is fake_factory
        finally:
            # Restore original state
            fac._EP_LOADED = original_ep_loaded
            fac._GEOMETRY_REGISTRY.clear()
            fac._GEOMETRY_REGISTRY.update(original_registry)

    # --- get_geometry() discovers entry points -----------------------------

    def test_get_geometry_finds_builtin_via_entry_point(self):
        """get_geometry() must resolve all built-in names."""
        for name in _BUILTIN_NAMES:
            factory = get_geometry(name)
            assert callable(factory)
            assert factory.__name__ == name

    # --- Broken plugin silently skipped ------------------------------------

    def test_broken_plugin_silently_skipped(self):
        """A plugin entry point whose load() raises must not crash list_geometries()."""
        from unittest.mock import MagicMock
        from unittest.mock import patch

        import ad_hoc_diffractometer.factories as fac

        broken_ep = MagicMock()
        broken_ep.name = "broken_plugin"
        broken_ep.load.side_effect = ImportError("simulated broken plugin")

        original_ep_loaded = fac._EP_LOADED
        original_registry = dict(fac._GEOMETRY_REGISTRY)
        fac._EP_LOADED = False

        try:
            with patch(
                "ad_hoc_diffractometer.factories.entry_points",
                return_value=[broken_ep],
            ):
                geoms = list_geometries()
            # broken plugin must not appear
            assert "broken_plugin" not in geoms
            # built-ins still present (registered via @register_geometry)
            assert "psic" in geoms
        finally:
            fac._EP_LOADED = original_ep_loaded
            fac._GEOMETRY_REGISTRY.clear()
            fac._GEOMETRY_REGISTRY.update(original_registry)

    # --- Duplicate name raises ValueError ----------------------------------

    def test_plugin_cannot_override_builtin_raises(self):
        """
        A plugin whose name collides with a built-in must raise ValueError,
        not silently win or lose.
        """
        import re
        from unittest.mock import MagicMock
        from unittest.mock import patch

        from helpers import fourcv

        import ad_hoc_diffractometer.factories as fac

        impostor_ep = MagicMock()
        impostor_ep.name = "psic"  # same name as built-in
        impostor_ep.value = "some_package.module:my_psic"
        impostor_ep.load.return_value = fourcv  # different callable

        original_ep_loaded = fac._EP_LOADED
        original_registry = dict(fac._GEOMETRY_REGISTRY)
        fac._EP_LOADED = False

        try:
            with patch(
                "ad_hoc_diffractometer.factories.entry_points",
                return_value=[impostor_ep],
            ):
                with pytest.raises(
                    ValueError,
                    match=re.escape("'psic' is already registered"),
                ):
                    list_geometries()
        finally:
            fac._EP_LOADED = original_ep_loaded
            fac._GEOMETRY_REGISTRY.clear()
            fac._GEOMETRY_REGISTRY.update(original_registry)

    def test_two_plugins_with_same_name_raises(self):
        """
        Two third-party plugins registering the same name must raise ValueError
        even if neither name collides with a built-in.
        """
        import re
        from unittest.mock import MagicMock
        from unittest.mock import patch

        from helpers import fourch
        from helpers import fourcv

        import ad_hoc_diffractometer.factories as fac

        ep1 = MagicMock()
        ep1.name = "my_custom_geom"
        ep1.value = "pkg_a.module:my_custom_geom"
        ep1.load.return_value = fourcv

        ep2 = MagicMock()
        ep2.name = "my_custom_geom"  # same name, different package
        ep2.value = "pkg_b.module:my_custom_geom"
        ep2.load.return_value = fourch

        original_ep_loaded = fac._EP_LOADED
        original_registry = dict(fac._GEOMETRY_REGISTRY)
        fac._EP_LOADED = False

        try:
            with patch(
                "ad_hoc_diffractometer.factories.entry_points",
                return_value=[ep1, ep2],
            ):
                with pytest.raises(
                    ValueError,
                    match=re.escape("'my_custom_geom' is already registered"),
                ):
                    list_geometries()
        finally:
            fac._EP_LOADED = original_ep_loaded
            fac._GEOMETRY_REGISTRY.clear()
            fac._GEOMETRY_REGISTRY.update(original_registry)

    def test_broken_plugin_logs_debug(self, caplog):
        """A broken entry point emits a DEBUG message to the package logger."""
        import logging
        from unittest.mock import MagicMock
        from unittest.mock import patch

        import ad_hoc_diffractometer.factories as _fac

        broken_ep = MagicMock()
        broken_ep.name = "broken_log_test"
        broken_ep.load.side_effect = ImportError("oops")

        original_ep_loaded = _fac._EP_LOADED
        original_registry = dict(_fac._GEOMETRY_REGISTRY)
        _fac._EP_LOADED = False
        try:
            with patch(
                "ad_hoc_diffractometer.factories.entry_points",
                return_value=[broken_ep],
            ):
                with caplog.at_level(
                    logging.DEBUG,
                    logger="ad_hoc_diffractometer.factories",
                ):
                    _fac._load_entry_point_geometries()
            assert any("broken_log_test" in r.message for r in caplog.records)
            assert any(r.levelno == logging.DEBUG for r in caplog.records)
        finally:
            _fac._EP_LOADED = original_ep_loaded
            _fac._GEOMETRY_REGISTRY.clear()
            _fac._GEOMETRY_REGISTRY.update(original_registry)

    def test_outer_importlib_exception_silently_ignored(self):
        """_load_entry_point_geometries swallows errors from entry_points() itself."""
        from unittest.mock import patch

        import ad_hoc_diffractometer.factories as _fac

        original_ep_loaded = _fac._EP_LOADED
        original_registry = dict(_fac._GEOMETRY_REGISTRY)
        _fac._EP_LOADED = False
        try:
            with patch(
                "ad_hoc_diffractometer.factories.entry_points",
                side_effect=Exception("metadata unavailable"),
            ):
                _fac._load_entry_point_geometries()
            assert "psic" in _fac._GEOMETRY_REGISTRY
        finally:
            _fac._EP_LOADED = original_ep_loaded
            _fac._GEOMETRY_REGISTRY.clear()
            _fac._GEOMETRY_REGISTRY.update(original_registry)


def test_packaged_geometry_loader_failure_is_non_fatal(monkeypatch):
    """If `_register_packaged_geometries` raises, registry init should
    not crash; the broken-loader exception is logged at debug level and
    swallowed (covers the broad except in factories.py)."""
    import ad_hoc_diffractometer.factories as fac

    fac._EP_LOADED = False  # noqa: SLF001 — force re-discovery
    original_registry = dict(fac._GEOMETRY_REGISTRY)  # noqa: SLF001
    try:
        # Monkey-patch the loader function the factories module imports
        # to raise when called.
        from ad_hoc_diffractometer import geometry_loader

        def _broken(*args, **kwargs):
            raise RuntimeError("simulated loader failure")

        monkeypatch.setattr(geometry_loader, "_register_packaged_geometries", _broken)
        # Also patch entry_points to an empty list so the rest of the
        # function runs without raising.
        monkeypatch.setattr(fac, "entry_points", lambda group: [])
        # Trigger the discovery path:
        fac._load_entry_point_geometries()  # noqa: SLF001
        # Built-ins from before this test must still be present
        # (the swallowed-exception branch did not corrupt the registry).
        for name in ("fourcv", "psic"):
            assert name in fac._GEOMETRY_REGISTRY  # noqa: SLF001
    finally:
        fac._EP_LOADED = True  # noqa: SLF001
        fac._GEOMETRY_REGISTRY.clear()  # noqa: SLF001
        fac._GEOMETRY_REGISTRY.update(original_registry)  # noqa: SLF001


def test_entry_point_redeclaring_existing_factory_is_skipped_silently():
    """An entry point whose .load() returns the same callable already in
    the registry is treated as a non-conflict and skipped (covers the
    dedup `continue` branch in _load_entry_point_geometries)."""
    from unittest.mock import MagicMock
    from unittest.mock import patch

    import ad_hoc_diffractometer.factories as fac
    from ad_hoc_diffractometer import geometry_loader

    # Trigger registry population (no-op if already populated).
    fac.list_geometries()
    original_ep_loaded = fac._EP_LOADED  # noqa: SLF001
    original_registry = dict(fac._GEOMETRY_REGISTRY)  # noqa: SLF001
    fac._EP_LOADED = False  # noqa: SLF001

    # Capture the existing registry callable for ``fourcv``.  The mock
    # entry point will return that exact callable, so the dedup
    # ``existing is factory`` check fires and the loop ``continue``s.
    duplicate_factory = fac._GEOMETRY_REGISTRY["fourcv"]
    fake_ep = MagicMock()
    fake_ep.name = "fourcv"
    fake_ep.value = "ad_hoc_diffractometer.geometries:fourcv.yml"
    fake_ep.load.return_value = duplicate_factory

    try:
        # Patch _register_packaged_geometries to a no-op so the
        # second invocation does not rebuild a new closure for fourcv;
        # we want the entry-point loop to see the original callable.
        with (
            patch.object(
                geometry_loader, "_register_packaged_geometries", lambda: None
            ),
            patch(
                "ad_hoc_diffractometer.factories.entry_points",
                return_value=[fake_ep],
            ),
        ):
            fac.list_geometries()  # triggers _load_entry_point_geometries
        # Registry still has fourcv pointing at the same factory
        assert fac._GEOMETRY_REGISTRY["fourcv"] is duplicate_factory
    finally:
        fac._EP_LOADED = original_ep_loaded
        fac._GEOMETRY_REGISTRY.clear()
        fac._GEOMETRY_REGISTRY.update(original_registry)
