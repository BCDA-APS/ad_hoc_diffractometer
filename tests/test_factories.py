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

from ad_hoc_diffractometer import XHAT
from ad_hoc_diffractometer import YHAT
from ad_hoc_diffractometer import ZHAT
from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer import fivec
from ad_hoc_diffractometer import fourch
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer import get_geometry
from ad_hoc_diffractometer import kappa4ch
from ad_hoc_diffractometer import kappa4cv
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
    """Keyword args must be forwarded to the factory (kappa alpha test)."""
    g = make_geometry("kappa4cv", alpha_deg=45.0)
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
            fourcv,
            "fourcv",
            ["omega", "chi", "phi"],
            ["two_theta"],
            does_not_raise(),
            id="fourcv-stage-lists",
        ),
        pytest.param(
            fourch,
            "fourch",
            ["omega", "chi", "phi"],
            ["two_theta"],
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
            ["two_theta"],
            does_not_raise(),
            id="kappa4cv-stage-lists",
        ),
        pytest.param(
            kappa4ch,
            "kappa4ch",
            ["komega", "kappa", "kphi"],
            ["two_theta"],
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
        # fourcv
        pytest.param(
            fourcv, "omega", None, does_not_raise(), id="fourcv-omega-on-floor"
        ),
        pytest.param(
            fourcv,
            "two_theta",
            None,
            does_not_raise(),
            id="fourcv-two_theta-decoupled",
        ),
        pytest.param(
            fourcv, "chi", "omega", does_not_raise(), id="fourcv-chi-on-omega"
        ),
        pytest.param(fourcv, "phi", "chi", does_not_raise(), id="fourcv-phi-on-chi"),
        # fourch
        pytest.param(
            fourch,
            "two_theta",
            None,
            does_not_raise(),
            id="fourch-two_theta-decoupled",
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
            "two_theta",
            None,
            does_not_raise(),
            id="kappa4cv-two_theta-decoupled",
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
            kappa4cv,
            "kappa",
            np.cos(np.deg2rad(50)) * np.array([0, 0, 1])
            + np.sin(np.deg2rad(50)) * np.array([1, 0, 0]),
            does_not_raise(),
            id="kappa4cv-kappa-axis-tilted",
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
    """kappa_alpha_deg is consistent with the kappa stage axis vector."""
    from ad_hoc_diffractometer import kappa_axis
    from ad_hoc_diffractometer.factories import _BASIS_BL
    from ad_hoc_diffractometer.factories import _BASIS_YOU

    with context:
        g = factory(alpha_deg=alpha_deg)
        kax = g.stage("kappa").axis
        basis = (
            _BASIS_YOU if factory is kappa6c else _BASIS_BL
        )  # kappa4cv/kappa4ch use BL
        expected_axis = kappa_axis(g.kappa_alpha_deg, basis=basis)
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
        from ad_hoc_diffractometer import GEOMETRY_ENTRY_POINT_GROUP

        assert GEOMETRY_ENTRY_POINT_GROUP == "ad_hoc_diffractometer.geometries"

    def test_entry_point_group_exported(self):
        """GEOMETRY_ENTRY_POINT_GROUP must be in __all__."""
        import ad_hoc_diffractometer as ahd

        assert "GEOMETRY_ENTRY_POINT_GROUP" in ahd.__all__

    # --- Built-in factories declared as entry points -----------------------

    def test_all_builtins_declared_as_entry_points(self):
        """All 10 built-in factories must appear in the installed entry points."""
        from importlib.metadata import entry_points

        from ad_hoc_diffractometer import GEOMETRY_ENTRY_POINT_GROUP

        eps = entry_points(group=GEOMETRY_ENTRY_POINT_GROUP)
        names = {ep.name for ep in eps}
        assert _BUILTIN_NAMES <= names, (
            f"Missing from entry points: {_BUILTIN_NAMES - names}"
        )

    def test_entry_point_loads_correct_factory(self):
        """Each built-in entry point must load to the same callable as the module."""
        from importlib.metadata import entry_points

        import ad_hoc_diffractometer.factories as fac
        from ad_hoc_diffractometer import GEOMETRY_ENTRY_POINT_GROUP

        eps = {ep.name: ep for ep in entry_points(group=GEOMETRY_ENTRY_POINT_GROUP)}
        for name in _BUILTIN_NAMES:
            assert name in eps
            loaded = eps[name].load()
            assert loaded is getattr(fac, name), (
                f"Entry point '{name}' loaded {loaded!r}, "
                f"expected {getattr(fac, name)!r}"
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

        import ad_hoc_diffractometer.factories as fac
        from ad_hoc_diffractometer import fourcv

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

    # --- Plugin does not override built-in ---------------------------------

    def test_plugin_cannot_override_builtin(self):
        """A plugin named 'psic' must not overwrite the built-in psic factory."""
        from unittest.mock import MagicMock
        from unittest.mock import patch

        import ad_hoc_diffractometer.factories as fac
        from ad_hoc_diffractometer import fourcv

        impostor_ep = MagicMock()
        impostor_ep.name = "psic"  # same name as built-in
        impostor_ep.load.return_value = fourcv  # but different callable

        original_ep_loaded = fac._EP_LOADED
        original_registry = dict(fac._GEOMETRY_REGISTRY)
        fac._EP_LOADED = False

        try:
            with patch(
                "ad_hoc_diffractometer.factories.entry_points",
                return_value=[impostor_ep],
            ):
                geoms = list_geometries()
            # The built-in psic must win — plugin cannot override it
            import ad_hoc_diffractometer.factories as fac2

            assert geoms["psic"] is fac2.psic
        finally:
            fac._EP_LOADED = original_ep_loaded
            fac._GEOMETRY_REGISTRY.clear()
            fac._GEOMETRY_REGISTRY.update(original_registry)
