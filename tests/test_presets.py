# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.presets.

Covers:
  - All 10 preset geometry functions are importable from ad_hoc_diffractometer.presets
  - presets.__all__ contains exactly the 10 preset names
  - Attribute-style access via ahd.presets.<name>() works
  - Preset functions are registered in the geometry registry
  - Presets are NOT importable from ad_hoc_diffractometer (top-level)

The preset *functions* themselves (factory output, stages, modes, etc.)
are tested exhaustively in test_factories.py.
"""

import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer import list_geometries

PRESET_NAMES = [
    "fivec",
    "fourch",
    "fourcv",
    "kappa4ch",
    "kappa4cv",
    "kappa6c",
    "psic",
    "s2d2",
    "sixc",
    "zaxis",
]


def test_presets_all_contains_exactly_10():
    """presets.__all__ has exactly the 10 preset geometry names."""
    import ad_hoc_diffractometer.presets as presets_mod

    assert sorted(presets_mod.__all__) == PRESET_NAMES


@pytest.mark.parametrize(
    "name",
    PRESET_NAMES,
    ids=PRESET_NAMES,
)
def test_preset_importable_from_presets_module(name):
    """Each preset is importable from ad_hoc_diffractometer.presets."""
    import ad_hoc_diffractometer.presets as presets_mod

    assert hasattr(presets_mod, name)
    factory = getattr(presets_mod, name)
    assert callable(factory)


@pytest.mark.parametrize(
    "name",
    PRESET_NAMES,
    ids=PRESET_NAMES,
)
def test_preset_not_in_top_level_namespace(name):
    """Presets are NOT importable from ad_hoc_diffractometer top-level."""
    assert name not in dir(ahd)
    assert not hasattr(ahd, name)


@pytest.mark.parametrize(
    "name",
    PRESET_NAMES,
    ids=PRESET_NAMES,
)
def test_preset_attribute_access_returns_diffractometer(name):
    """ahd.presets.<name>() returns an AdHocDiffractometer."""
    factory = getattr(ahd.presets, name)
    g = factory()
    assert isinstance(g, AdHocDiffractometer)
    assert g.name == name


@pytest.mark.parametrize(
    "name",
    PRESET_NAMES,
    ids=PRESET_NAMES,
)
def test_preset_registered_in_geometry_registry(name):
    """Each preset is registered and discoverable via list_geometries()."""
    geoms = list_geometries()
    assert name in geoms
