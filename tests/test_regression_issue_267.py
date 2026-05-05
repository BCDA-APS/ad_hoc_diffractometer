# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Regression tests for issue #267 — declarative geometries.

For every demo geometry that has been migrated from a Python factory to
a declarative YAML file, this module verifies that the loader-built
geometry is **structurally equivalent** to a hand-built reference
geometry (the legacy Python factory).  This is the primary safety net
for the migration.

The current commit migrates only ``fourcv``.  Subsequent commits will
extend the parametrisation as the remaining nine demo geometries are
moved to YAML.
"""

from __future__ import annotations

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer.factories import BASIS_BL
from ad_hoc_diffractometer.factories import BASIS_YOU
from ad_hoc_diffractometer.mode import REQUIRED
from ad_hoc_diffractometer.mode import BisectConstraint
from ad_hoc_diffractometer.mode import ConstraintSet
from ad_hoc_diffractometer.mode import ReferenceConstraint
from ad_hoc_diffractometer.mode import SampleConstraint
from ad_hoc_diffractometer.stage import Stage

# ---------------------------------------------------------------------------
# Hand-built reference geometries
# ---------------------------------------------------------------------------


def _reference_fourcv() -> AdHocDiffractometer:
    """Hand-built fourcv (the pre-#267 ``presets.fourcv`` body verbatim)."""
    basis = BASIS_BL
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("omega", -TRANSVERSE, parent=None, role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -TRANSVERSE, parent="chi", role="sample"),
        Stage("ttheta", -TRANSVERSE, parent=None, role="detector"),
    ]
    modes = {
        "bisecting": ConstraintSet(
            [BisectConstraint("omega", "ttheta")],
            computed=["omega", "chi", "phi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["omega", "phi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["omega", "chi", "ttheta"],
        ),
        "fixed_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["chi", "phi", "ttheta"],
        ),
        "fixed_psi": ConstraintSet(
            [ReferenceConstraint("psi", 0.0)],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction": ConstraintSet(
            [],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
    }
    return AdHocDiffractometer(
        name="fourcv",
        stages=stages,
        basis=BASIS_BL,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(vertical scattering plane, transverse ttheta, synchrotron)"
        ),
        modes=modes,
        default_mode="bisecting",
    )


def _reference_fourch() -> AdHocDiffractometer:
    """Hand-built fourch (the pre-#267 ``presets.fourch`` body verbatim)."""
    basis = BASIS_BL
    VERTICAL = basis["vertical"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("omega", -VERTICAL, parent=None, role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -VERTICAL, parent="chi", role="sample"),
        Stage("ttheta", -VERTICAL, parent=None, role="detector"),
    ]
    modes = {
        "bisecting": ConstraintSet(
            [BisectConstraint("omega", "ttheta")],
            computed=["omega", "chi", "phi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["omega", "phi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["omega", "chi", "ttheta"],
        ),
        "fixed_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["chi", "phi", "ttheta"],
        ),
        "fixed_psi": ConstraintSet(
            [ReferenceConstraint("psi", 0.0)],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction": ConstraintSet(
            [],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
    }
    return AdHocDiffractometer(
        name="fourch",
        stages=stages,
        basis=basis,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(horizontal scattering plane, vertical ttheta, laboratory)"
        ),
        modes=modes,
        default_mode="bisecting",
    )


# ---------------------------------------------------------------------------
# Equivalence checker
# ---------------------------------------------------------------------------


def _assert_geometries_equivalent(
    declarative: AdHocDiffractometer,
    reference: AdHocDiffractometer,
) -> None:
    """Compare two geometries field-by-field for structural equivalence."""
    assert declarative.name == reference.name
    # basis
    for key in ("vertical", "longitudinal", "transverse"):
        np.testing.assert_array_equal(
            declarative.basis[key],
            reference.basis[key],
            err_msg=f"basis['{key}'] differs",
        )
    # stages
    decl_stages = list(declarative._stages.values())  # noqa: SLF001
    ref_stages = list(reference._stages.values())  # noqa: SLF001
    assert [s.name for s in decl_stages] == [s.name for s in ref_stages]
    for ds, rs in zip(decl_stages, ref_stages, strict=False):
        assert ds.name == rs.name
        np.testing.assert_array_equal(
            ds.axis,
            rs.axis,
            err_msg=f"stage {ds.name!r} axis differs",
        )
        assert ds.parent == rs.parent, f"stage {ds.name!r} parent differs"
        assert ds.role == rs.role, f"stage {ds.name!r} role differs"
    # modes
    assert list(declarative.modes) == list(reference.modes)
    for mname, decl_cs in declarative.modes.items():
        ref_cs = reference.modes[mname]
        assert decl_cs.constraints == ref_cs.constraints, (
            f"mode {mname!r} constraints differ: "
            f"{decl_cs.constraints!r} vs {ref_cs.constraints!r}"
        )
        assert decl_cs.computed == ref_cs.computed, f"mode {mname!r} computed differs"
        # extras: same keys, identical sentinel mapping for REQUIRED
        assert set(decl_cs.extras) == set(ref_cs.extras), (
            f"mode {mname!r} extras keys differ"
        )
        for ek in decl_cs.extras:
            assert decl_cs.extras[ek] is ref_cs.extras[ek] or (
                decl_cs.extras[ek] == ref_cs.extras[ek]
            ), f"mode {mname!r} extras[{ek!r}] differs"
    # default mode
    assert declarative.mode_name == reference.mode_name


# ---------------------------------------------------------------------------
# Per-geometry equivalence tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "geom_name, reference_factory",
    [
        pytest.param("fourcv", _reference_fourcv, id="fourcv"),
        pytest.param("fourch", _reference_fourch, id="fourch"),
    ],
)
def test_declarative_matches_reference(geom_name, reference_factory):
    """The declarative YAML build must equal the hand-built reference."""
    declarative = ahd.make_geometry(geom_name)
    reference = reference_factory()
    _assert_geometries_equivalent(declarative, reference)


@pytest.mark.parametrize(
    "geom_name",
    [
        pytest.param("fourcv", id="fourcv"),
        pytest.param("fourch", id="fourch"),
    ],
)
def test_declarative_basis_override_round_trips(geom_name):
    """Caller-supplied ``basis=`` reaches the constructed geometry."""
    g_default = ahd.make_geometry(geom_name)
    # Pick a basis that is *different* from the file's declared default
    # (which is BASIS_BL for fourcv); BASIS_YOU is a clean test choice.
    other = BASIS_YOU
    if all(
        np.array_equal(g_default.basis[k], other[k])
        for k in ("vertical", "longitudinal", "transverse")
    ):
        pytest.skip("default basis already equals override")
    g_override = ahd.make_geometry(geom_name, basis=other)
    for k in ("vertical", "longitudinal", "transverse"):
        np.testing.assert_array_equal(g_override.basis[k], other[k])
