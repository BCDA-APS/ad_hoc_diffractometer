# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.mode and mode-related geometry features.

Covers:
  - BisectConstraint: construction, evaluate, is_implemented (both stages present /
    one missing / wrong category), to_dict/from_dict, __repr__, __eq__, __hash__
  - SampleConstraint: construction, evaluate, is_implemented,
    to_dict/from_dict, __repr__, __eq__
  - DetectorConstraint: construction, evaluate, is_qaz, is_implemented,
    to_dict/from_dict, __repr__, __eq__
  - ReferenceConstraint: construction, invalid names, is_implemented,
    to_dict/from_dict, __repr__, __eq__
  - ConstraintSet: construction, taxonomy validation (at most 1 bisect/det/ref),
    has_bisect, bisect_constraint, sample_constraints, detector_constraint,
    reference_constraint, fixed_sample_constraints, is_fully_constrained,
    is_implemented, constrained_stages, bisect_stages(), apply_cut_point,
    to_dict/from_dict, __repr__, __eq__, __len__
  - ModeDict: construction, set/get/delete, type guard, __len__, __iter__,
    keys/values/items, __repr__, __eq__, __contains__
  - AdHocDiffractometer: modes, default_mode, mode_name, mode property,
    free_dof_after_bragg, cut_points, mode_name setter validation
  - Factory modes: psic, fourcv, fourch, kappa4cv, kappa4ch, kappa6c
  - Geometries without modes: sixc, zaxis, s2d2, fivec
  - Serialisation round-trip: to_dict / from_dict for ConstraintSet,
    ModeDict on geometry, cut_points; error on old-format types
  - REQUIRED/OPTIONAL sentinels in extras survive round-trip
"""

import re
from contextlib import nullcontext as does_not_raise

import pytest

from ad_hoc_diffractometer import OPTIONAL
from ad_hoc_diffractometer import REFERENCE_NAMES
from ad_hoc_diffractometer import REQUIRED
from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer import BisectConstraint
from ad_hoc_diffractometer import ConstraintSet
from ad_hoc_diffractometer import ConstraintViolation
from ad_hoc_diffractometer import DetectorConstraint
from ad_hoc_diffractometer import ModeDict
from ad_hoc_diffractometer import ReferenceConstraint
from ad_hoc_diffractometer import SampleConstraint
from ad_hoc_diffractometer import Stage
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer.constants import XHAT
from ad_hoc_diffractometer.constants import YHAT
from ad_hoc_diffractometer.constants import ZHAT

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _simple_geometry(**kwargs):
    """Minimal 2-stage geometry for testing mode features (N=2, free_dof=-1)."""
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


def _fourcv_like(**kwargs):
    """3-sample + 1-detector geometry (N=4, free_dof=1)."""
    stages = [
        Stage("omega", -ZHAT, parent=None, role="sample"),
        Stage("chi", +YHAT, parent="omega", role="sample"),
        Stage("phi", -ZHAT, parent="chi", role="sample"),
        Stage("ttheta", -ZHAT, parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name="fourcv_like",
        stages=stages,
        basis={"vertical": XHAT, "longitudinal": YHAT, "lateral": ZHAT},
        **kwargs,
    )


def _psic_like(**kwargs):
    """4-sample + 2-detector geometry (N=6, free_dof=3)."""
    stages = [
        Stage("mu", +XHAT, parent=None, role="sample"),
        Stage("eta", -ZHAT, parent="mu", role="sample"),
        Stage("chi", +YHAT, parent="eta", role="sample"),
        Stage("phi", -ZHAT, parent="chi", role="sample"),
        Stage("nu", +XHAT, parent=None, role="detector"),
        Stage("delta", -ZHAT, parent="nu", role="detector"),
    ]
    return AdHocDiffractometer(
        name="psic_like",
        stages=stages,
        basis={"vertical": XHAT, "longitudinal": YHAT, "lateral": ZHAT},
        **kwargs,
    )


# ---------------------------------------------------------------------------
# SampleConstraint
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, value, context",
    [
        pytest.param("chi", 90.0, does_not_raise(), id="fixed-chi-90"),
        pytest.param("phi", 0.0, does_not_raise(), id="fixed-phi-0"),
        pytest.param("mu", -45.0, does_not_raise(), id="fixed-mu-neg45"),
    ],
)
def test_sample_constraint_construction(name, value, context):
    with context:
        sc = SampleConstraint(name, value)
        assert sc.name == name
        assert sc.is_bisect is False
        assert sc.category == "sample"


def test_sample_constraint_value_coerced_to_float():
    sc = SampleConstraint("chi", 90)
    assert isinstance(sc.value, float)
    assert sc.value == 90.0


def test_sample_constraint_repr_fixed():
    sc = SampleConstraint("chi", 90.0)
    r = repr(sc)
    assert "SampleConstraint" in r
    assert "chi" in r
    assert "90.0" in r


def test_sample_constraint_eq_same():
    assert SampleConstraint("chi", 90.0) == SampleConstraint("chi", 90.0)


def test_sample_constraint_eq_different_name():
    assert SampleConstraint("chi", 90.0) != SampleConstraint("phi", 90.0)


def test_sample_constraint_eq_different_value():
    assert SampleConstraint("chi", 90.0) != SampleConstraint("chi", 0.0)


def test_sample_constraint_eq_different_type():
    assert SampleConstraint("chi", 90.0) != DetectorConstraint("chi", 90.0)


def test_sample_constraint_to_dict_from_dict_fixed():
    sc = SampleConstraint("chi", 90.0)
    d = sc.to_dict()
    assert d == {"type": "SampleConstraint", "name": "chi", "value": 90.0}
    sc2 = SampleConstraint.from_dict(d)
    assert sc2 == sc


def test_bisect_constraint_to_dict_from_dict():
    bc = BisectConstraint("omega", "ttheta")
    d = bc.to_dict()
    assert d == {
        "type": "BisectConstraint",
        "sample_stage": "omega",
        "detector_stage": "ttheta",
    }
    bc2 = BisectConstraint.from_dict(d)
    assert bc2 == bc


def test_sample_constraint_is_implemented_real_stage():
    g = _fourcv_like()
    assert SampleConstraint("chi", 90.0).is_implemented(g) is True
    assert SampleConstraint("phi", 0.0).is_implemented(g) is True


def test_sample_constraint_is_implemented_missing_stage():
    g = _fourcv_like()
    assert SampleConstraint("mu", 0.0).is_implemented(g) is False


def test_bisect_constraint_is_implemented_both_stages_exist():
    g = _fourcv_like()
    assert BisectConstraint("omega", "ttheta").is_implemented(g) is True


def test_bisect_constraint_is_implemented_sample_stage_missing():
    g = _fourcv_like()
    assert BisectConstraint("eta", "ttheta").is_implemented(g) is False


def test_bisect_constraint_is_implemented_detector_stage_missing():
    g = _fourcv_like()
    assert BisectConstraint("omega", "delta").is_implemented(g) is False


def test_bisect_constraint_is_implemented_psic():
    g = _psic_like()
    assert BisectConstraint("eta", "delta").is_implemented(g) is True


def test_sample_constraint_evaluate_fixed(tmp_path):
    g = _fourcv_like()
    angles = {"omega": 10.0, "chi": 90.0, "phi": 0.0, "ttheta": 20.0}
    sc = SampleConstraint("chi", 90.0)
    assert sc.evaluate(angles, g) == pytest.approx(0.0)
    sc2 = SampleConstraint("chi", 45.0)
    assert sc2.evaluate(angles, g) == pytest.approx(45.0)


def test_bisect_constraint_evaluate_satisfied():
    g = _fourcv_like()
    angles = {"omega": 10.0, "chi": 90.0, "phi": 0.0, "ttheta": 20.0}
    bc = BisectConstraint("omega", "ttheta")
    # residual = omega - ttheta/2 = 10 - 10 = 0
    assert bc.evaluate(angles, g) == pytest.approx(0.0)


def test_bisect_constraint_evaluate_not_satisfied():
    g = _fourcv_like()
    angles = {"omega": 5.0, "chi": 90.0, "phi": 0.0, "ttheta": 20.0}
    bc = BisectConstraint("omega", "ttheta")
    assert bc.evaluate(angles, g) == pytest.approx(-5.0)


def test_bisect_constraint_is_satisfied():
    g = _fourcv_like()
    angles = {"omega": 10.0, "chi": 90.0, "phi": 0.0, "ttheta": 20.0}
    assert BisectConstraint("omega", "ttheta").is_satisfied(angles, g)
    angles2 = dict(angles, omega=5.0)
    assert not BisectConstraint("omega", "ttheta").is_satisfied(angles2, g)


def test_bisect_constraint_repr():
    bc = BisectConstraint("eta", "delta")
    r = repr(bc)
    assert "BisectConstraint" in r
    assert "eta" in r
    assert "delta" in r


def test_bisect_constraint_eq_same():
    assert BisectConstraint("omega", "ttheta") == BisectConstraint("omega", "ttheta")


def test_bisect_constraint_eq_different():
    assert BisectConstraint("omega", "ttheta") != BisectConstraint("eta", "delta")


def test_bisect_constraint_eq_different_type():
    assert BisectConstraint("omega", "ttheta") != SampleConstraint("omega", 0.0)


def test_bisect_constraint_hash():
    bc1 = BisectConstraint("omega", "ttheta")
    bc2 = BisectConstraint("omega", "ttheta")
    assert hash(bc1) == hash(bc2)
    assert len({bc1, bc2}) == 1


def test_bisect_constraint_category():
    assert BisectConstraint("omega", "ttheta").category == "sample"


def test_bisect_constraint_is_bisect_flag():
    assert BisectConstraint("omega", "ttheta").is_bisect is True


def test_bisect_constraint_name():
    assert BisectConstraint("omega", "ttheta").name == "bisect"


def test_sample_constraint_is_satisfied():
    g = _fourcv_like()
    angles = {"omega": 10.0, "chi": 90.0, "phi": 0.0, "ttheta": 20.0}
    assert SampleConstraint("chi", 90.0).is_satisfied(angles, g)
    assert not SampleConstraint("chi", 45.0).is_satisfied(angles, g)


# ---------------------------------------------------------------------------
# DetectorConstraint
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, value, expected_is_qaz, context",
    [
        pytest.param("delta", 0.0, False, does_not_raise(), id="fixed-delta"),
        pytest.param("nu", 0.0, False, does_not_raise(), id="fixed-nu"),
        pytest.param("ttheta", 20.0, False, does_not_raise(), id="fixed-ttheta"),
        pytest.param("qaz", 90.0, True, does_not_raise(), id="qaz-pseudo"),
    ],
)
def test_detector_constraint_construction(name, value, expected_is_qaz, context):
    with context:
        dc = DetectorConstraint(name, value)
        assert dc.name == name
        assert dc.value == value
        assert dc.is_qaz == expected_is_qaz
        assert dc.category == "detector"


def test_detector_constraint_repr():
    dc = DetectorConstraint("nu", 0.0)
    r = repr(dc)
    assert "DetectorConstraint" in r
    assert "nu" in r


def test_detector_constraint_eq_same():
    assert DetectorConstraint("nu", 0.0) == DetectorConstraint("nu", 0.0)


def test_detector_constraint_eq_different():
    assert DetectorConstraint("nu", 0.0) != DetectorConstraint("delta", 0.0)


def test_detector_constraint_eq_different_type():
    assert DetectorConstraint("nu", 0.0) != SampleConstraint("nu", 0.0)


def test_detector_constraint_to_dict_from_dict():
    dc = DetectorConstraint("nu", 0.0)
    d = dc.to_dict()
    assert d == {"type": "DetectorConstraint", "name": "nu", "value": 0.0}
    dc2 = DetectorConstraint.from_dict(d)
    assert dc2 == dc


def test_detector_constraint_is_implemented_real_stage():
    g = _psic_like()
    assert DetectorConstraint("nu", 0.0).is_implemented(g) is True
    assert DetectorConstraint("delta", 0.0).is_implemented(g) is True


def test_detector_constraint_is_implemented_missing():
    g = _psic_like()
    assert DetectorConstraint("gamma", 0.0).is_implemented(g) is False


def test_detector_constraint_is_implemented_qaz_not_implemented():
    g = _psic_like()
    assert DetectorConstraint("qaz", 90.0).is_implemented(g) is False


def test_detector_constraint_evaluate_fixed():
    g = _psic_like()
    angles = {"mu": 0.0, "eta": 10.0, "chi": 90.0, "phi": 0.0, "nu": 0.0, "delta": 20.0}
    dc = DetectorConstraint("nu", 0.0)
    assert dc.evaluate(angles, g) == pytest.approx(0.0)
    dc2 = DetectorConstraint("nu", 5.0)
    assert dc2.evaluate(angles, g) == pytest.approx(-5.0)


def test_detector_constraint_evaluate_qaz_raises():
    g = _psic_like()
    angles = {}
    dc = DetectorConstraint("qaz", 90.0)
    with pytest.raises(NotImplementedError, match="qaz"):
        dc.evaluate(angles, g)


# ---------------------------------------------------------------------------
# ReferenceConstraint
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, value, context",
    [
        pytest.param("psi", 90.0, does_not_raise(), id="psi"),
        pytest.param("alpha_i", 5.0, does_not_raise(), id="alpha_i"),
        pytest.param("beta_out", 5.0, does_not_raise(), id="beta_out"),
        pytest.param("a_eq_b", True, does_not_raise(), id="a_eq_b-true"),
        pytest.param("naz", 0.0, does_not_raise(), id="naz"),
        pytest.param(
            "a_eq_b",
            False,
            pytest.raises(ValueError, match=re.escape("must be True")),
            id="a_eq_b-false-raises",
        ),
        pytest.param(
            "unknown",
            0.0,
            pytest.raises(ValueError, match=re.escape("must be one of")),
            id="unknown-name-raises",
        ),
    ],
)
def test_reference_constraint_construction(name, value, context):
    with context:
        rc = ReferenceConstraint(name, value)
        assert rc.name == name
        assert rc.category == "reference"


def test_reference_constraint_value_coerced_to_float():
    rc = ReferenceConstraint("psi", 90)
    assert isinstance(rc.value, float)
    assert rc.value == 90.0


def test_reference_constraint_a_eq_b_value_is_true():
    rc = ReferenceConstraint("a_eq_b", True)
    assert rc.value is True


def test_reference_constraint_repr():
    rc = ReferenceConstraint("psi", 90.0)
    r = repr(rc)
    assert "ReferenceConstraint" in r
    assert "psi" in r


def test_reference_constraint_eq_same():
    assert ReferenceConstraint("psi", 90.0) == ReferenceConstraint("psi", 90.0)


def test_reference_constraint_eq_different():
    assert ReferenceConstraint("psi", 90.0) != ReferenceConstraint("psi", 0.0)


def test_reference_constraint_eq_different_type():
    assert ReferenceConstraint("psi", 90.0) != SampleConstraint("psi", 90.0)


def test_reference_constraint_to_dict_from_dict():
    rc = ReferenceConstraint("psi", 90.0)
    d = rc.to_dict()
    assert d == {"type": "ReferenceConstraint", "name": "psi", "value": 90.0}
    rc2 = ReferenceConstraint.from_dict(d)
    assert rc2 == rc


def test_reference_constraint_to_dict_from_dict_a_eq_b():
    rc = ReferenceConstraint("a_eq_b", True)
    d = rc.to_dict()
    assert d["value"] is True
    rc2 = ReferenceConstraint.from_dict(d)
    assert rc2 == rc


def test_reference_constraint_is_implemented_always_false():
    g = _psic_like()
    for name in REFERENCE_NAMES:
        value = True if name == "a_eq_b" else 0.0
        rc = ReferenceConstraint(name, value)
        assert rc.is_implemented(g) is False


def test_reference_constraint_evaluate_raises():
    g = _psic_like()
    rc = ReferenceConstraint("psi", 90.0)
    with pytest.raises(NotImplementedError, match="not yet implemented"):
        rc.evaluate({}, g)


def test_reference_constraint_is_satisfied_raises():
    g = _psic_like()
    rc = ReferenceConstraint("psi", 90.0)
    with pytest.raises(NotImplementedError):
        rc.is_satisfied({}, g)


# ---------------------------------------------------------------------------
# ConstraintSet — construction and taxonomy validation
# ---------------------------------------------------------------------------


def test_constraint_set_basic_construction():
    cs = ConstraintSet([BisectConstraint("omega", "ttheta")])
    assert len(cs) == 1
    assert cs.has_bisect


def test_constraint_set_three_constraints():
    cs = ConstraintSet(
        [
            BisectConstraint("omega", "ttheta"),
            SampleConstraint("mu", 0.0),
            DetectorConstraint("nu", 0.0),
        ]
    )
    assert len(cs) == 3


def test_constraint_set_two_detector_raises():
    with pytest.raises(ValueError, match=re.escape("at most one DetectorConstraint")):
        ConstraintSet(
            [
                DetectorConstraint("nu", 0.0),
                DetectorConstraint("delta", 0.0),
            ]
        )


def test_constraint_set_two_reference_raises():
    with pytest.raises(ValueError, match=re.escape("at most one ReferenceConstraint")):
        ConstraintSet(
            [
                ReferenceConstraint("psi", 0.0),
                ReferenceConstraint("alpha_i", 5.0),
            ]
        )


def test_constraint_set_two_bisect_raises():
    with pytest.raises(ValueError, match=re.escape("at most one BisectConstraint")):
        ConstraintSet(
            [
                BisectConstraint("omega", "ttheta"),
                BisectConstraint("eta", "delta"),
            ]
        )


def test_constraint_set_invalid_type_raises():
    """ConstraintSet raises ValueError for non-constraint items."""
    with pytest.raises(ValueError, match=re.escape("must be BisectConstraint")):
        ConstraintSet(["not_a_constraint"])  # type: ignore[list-item]


def test_constraint_set_constrained_stages_qaz_det_excluded():
    """qaz DetectorConstraint contributes no stage name to constrained_stages."""
    g = _psic_like()
    cs = ConstraintSet(
        [
            BisectConstraint("omega", "ttheta"),
            SampleConstraint("mu", 0.0),
            DetectorConstraint("qaz", 90.0),
        ]
    )
    stages = cs.constrained_stages(g)
    assert "qaz" not in stages
    # eta should be there from bisect; mu from fixed
    assert "mu" in stages


def test_bisect_constraint_explicit_stages_no_axis_heuristic():
    """BisectConstraint uses declared stage names directly, no axis-geometry heuristic.

    Even with anti-parallel axes, the declared names are used as-is.
    """
    # A geometry where sample axis is -ZHAT and detector is +ZHAT (anti-parallel)
    stages = [
        Stage("omega", -ZHAT, parent=None, role="sample"),
        Stage("ttheta", +ZHAT, parent=None, role="detector"),
    ]
    g = AdHocDiffractometer(
        name="antipar",
        stages=stages,
        basis={"vertical": XHAT, "longitudinal": YHAT, "lateral": ZHAT},
    )
    bc = BisectConstraint("omega", "ttheta")
    assert bc.is_implemented(g) is True
    assert bc.sample_stage == "omega"
    assert bc.detector_stage == "ttheta"


def test_reference_constraint_hash():
    """ReferenceConstraint is hashable."""
    rc1 = ReferenceConstraint("psi", 90.0)
    rc2 = ReferenceConstraint("psi", 90.0)
    assert hash(rc1) == hash(rc2)
    # a_eq_b with bool value also hashable
    rc3 = ReferenceConstraint("a_eq_b", True)
    assert isinstance(hash(rc3), int)


def test_constraint_set_constraints_property():
    """ConstraintSet.constraints property returns a copy of the constraint list."""
    cs = ConstraintSet([SampleConstraint("chi", 90.0)])
    c_list = cs.constraints
    assert len(c_list) == 1
    assert isinstance(c_list[0], SampleConstraint)
    # Modifying the returned list does not affect the ConstraintSet
    c_list.append(SampleConstraint("phi", 0.0))
    assert len(cs.constraints) == 1


def test_constraint_set_extras_plain_value_round_trip():
    """Extras with plain values (not sentinels) survive round-trip serialisation."""
    cs = ConstraintSet(
        [BisectConstraint("omega", "ttheta")],
        extras={"some_float": 42.0, "psi": None, "n_hat": REQUIRED},
    )
    d = cs.to_dict()
    # plain float value
    assert d["extras"]["some_float"] == 42.0
    cs2 = ConstraintSet.from_dict(d)
    assert cs2.extras["some_float"] == 42.0


def test_constraint_set_from_dict_with_reference_constraint():
    """ConstraintSet.from_dict handles ReferenceConstraint entries."""
    d = {
        "type": "ConstraintSet",
        "constraints": [
            {"type": "SampleConstraint", "name": "bisect", "value": True},
            {"type": "ReferenceConstraint", "name": "psi", "value": 90.0},
        ],
        "computed": None,
        "constant": {},
        "extras": {},
        "cut_points": {},
    }
    cs = ConstraintSet.from_dict(d)
    assert cs.reference_constraint is not None
    assert cs.reference_constraint.name == "psi"


def test_constraint_set_constrained_stages_duplicate_avoided():
    """constrained_stages does not duplicate bisect stage name."""
    g = _psic_like()
    # mu is both the bisect stage (mu+nu coaxial) and a fixed constraint
    # Since bisect is innermost-det-based, eta is the bisect sample for psic_like
    cs = ConstraintSet(
        [
            BisectConstraint("omega", "ttheta"),
            SampleConstraint("eta", 5.0),  # also fix eta explicitly
            DetectorConstraint("nu", 0.0),
        ]
    )
    stages = cs.constrained_stages(g)
    # eta should appear only once even if bisect returns eta too
    assert stages.count("eta") == 1


def test_constraint_set_constrained_stages_det_duplicate_avoided():
    """DetectorConstraint name already in names list is not duplicated."""
    g = _fourcv_like()
    cs = ConstraintSet(
        [
            SampleConstraint("omega", 0.0),
            SampleConstraint("chi", 90.0),
        ]
    )
    stages = cs.constrained_stages(g)
    assert stages.count("omega") == 1
    assert stages.count("chi") == 1


def test_constraint_set_constrained_stages_bisect_dedup():
    """bisect sample stage name not duplicated when it also appears in fixed constraints."""
    # In a psic-like geometry, bisect resolves to "eta".
    # If we also have SampleConstraint("eta", X), eta should only appear once.
    g = _psic_like()
    cs = ConstraintSet(
        [
            BisectConstraint("omega", "ttheta"),
            SampleConstraint("eta", 5.0),  # same name as bisect stage
            DetectorConstraint("nu", 0.0),
        ]
    )
    stages = cs.constrained_stages(g)
    # eta from bisect, but then also from fixed — dedup guard at line 726 should prevent double
    assert stages.count("eta") == 1


def test_constraint_set_constrained_stages_det_name_present():
    """DetectorConstraint name is added to constrained_stages list."""
    cs = ConstraintSet(
        [
            SampleConstraint("mu", 0.0),
            DetectorConstraint("nu", 0.0),
        ]
    )
    stages_list = cs.constrained_stages()
    assert "mu" in stages_list
    assert stages_list.count("nu") == 1


# ---------------------------------------------------------------------------
# constant_stages property
# ---------------------------------------------------------------------------


def test_constant_stages_bisect_only():
    """constant_stages includes the bisect sample stage."""
    cs = ConstraintSet([BisectConstraint("omega", "ttheta")])
    assert cs.constant_stages == ["omega"]


def test_constant_stages_psic_bisecting():
    """constant_stages for psic bisecting includes eta, mu, nu."""
    cs = ConstraintSet(
        [
            BisectConstraint("eta", "delta"),
            SampleConstraint("mu", 0.0),
            DetectorConstraint("nu", 0.0),
        ]
    )
    assert cs.constant_stages == ["eta", "mu", "nu"]


def test_constant_stages_no_bisect():
    """constant_stages with only fixed sample and detector."""
    cs = ConstraintSet(
        [
            SampleConstraint("chi", 90.0),
            DetectorConstraint("nu", 0.0),
        ]
    )
    assert "chi" in cs.constant_stages
    assert "nu" in cs.constant_stages


def test_constant_stages_reference_excluded():
    """ReferenceConstraint does not contribute to constant_stages."""
    cs = ConstraintSet(
        [
            BisectConstraint("omega", "ttheta"),
            ReferenceConstraint("psi", 90.0),
        ]
    )
    stages = cs.constant_stages
    assert "psi" not in stages
    assert "omega" in stages


def test_constant_stages_qaz_excluded():
    """DetectorConstraint qaz does not contribute to constant_stages."""
    cs = ConstraintSet(
        [
            SampleConstraint("chi", 90.0),
            DetectorConstraint("qaz", 90.0),
        ]
    )
    stages = cs.constant_stages
    assert "qaz" not in stages
    assert "chi" in stages


# ---------------------------------------------------------------------------
# _validate_solutions — post-solve constraint checking
# ---------------------------------------------------------------------------


def test_validate_solutions_constraint_violation_raises():
    """forward() raises ValueError when solver returns a solution violating a constraint."""
    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer import ub_identity

    g = fourcv()
    g.wavelength = 1.54
    g.sample.lattice = ahd.Lattice(a=4.0)
    ub_identity(g.sample)

    # Create a mode with SampleConstraint("chi", 90.0) but set chi to the
    # wrong value in the solver output by injecting a patched mode.
    # The easiest way: use a SampleConstraint with an impossible constraint
    # value that the solver will violate, by using a mode that is "implemented"
    # but whose constraint is incompatible with the solution.
    # We can't easily force this through the normal API, so test via
    # _validate_solutions directly.
    from ad_hoc_diffractometer.forward import _validate_solutions
    from ad_hoc_diffractometer.mode import BisectConstraint
    from ad_hoc_diffractometer.mode import ConstraintSet
    from ad_hoc_diffractometer.mode import SampleConstraint

    mode = ConstraintSet(
        [
            BisectConstraint("omega", "ttheta"),
            SampleConstraint("chi", 90.0),
        ]
    )

    # Solution that violates chi constraint (chi=45 instead of 90)
    bad_sol = {"omega": 10.0, "chi": 45.0, "phi": 0.0, "ttheta": 20.0}

    with pytest.raises(ConstraintViolation, match=re.escape("violates")):
        _validate_solutions([bad_sol], mode, g)


def test_validate_solutions_keyerror_skipped():
    """_validate_solutions silently skips constraints whose stage is absent from solution."""
    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer.forward import _validate_solutions
    from ad_hoc_diffractometer.mode import BisectConstraint
    from ad_hoc_diffractometer.mode import ConstraintSet

    g = fourcv()
    g.wavelength = 1.54
    g.sample.lattice = ahd.Lattice(a=4.0)

    mode = ConstraintSet([BisectConstraint("omega", "ttheta")])

    # Solution missing "ttheta" — KeyError should be caught and skipped
    incomplete_sol = {"omega": 10.0, "chi": 90.0, "phi": 0.0}  # no ttheta
    # Should not raise
    _validate_solutions([incomplete_sol], mode, g)


def test_validate_solutions_empty_list():
    """_validate_solutions with empty solution list does nothing."""
    from ad_hoc_diffractometer.forward import _validate_solutions
    from ad_hoc_diffractometer.mode import BisectConstraint
    from ad_hoc_diffractometer.mode import ConstraintSet

    g = fourcv()
    mode = ConstraintSet([BisectConstraint("omega", "ttheta")])
    # Empty solutions — should not raise
    _validate_solutions([], mode, g)


def test_validate_solutions_bisect_violation_message():
    """Violation message for BisectConstraint mentions sample and detector stage names."""
    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer.forward import _validate_solutions
    from ad_hoc_diffractometer.mode import BisectConstraint
    from ad_hoc_diffractometer.mode import ConstraintSet

    g = fourcv()
    g.wavelength = 1.54
    g.sample.lattice = ahd.Lattice(a=4.0)

    mode = ConstraintSet([BisectConstraint("omega", "ttheta")])
    # omega=5 but ttheta=20 → bisect residual = 5 - 10 = -5 (violation)
    bad_sol = {"omega": 5.0, "chi": 90.0, "phi": 0.0, "ttheta": 20.0}

    with pytest.raises(ConstraintViolation, match=re.escape("BisectConstraint")):
        _validate_solutions([bad_sol], mode, g)


def test_validate_solutions_sampleconstraint_violation_message():
    """Violation message for SampleConstraint mentions the stage name and value."""
    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer.forward import _validate_solutions
    from ad_hoc_diffractometer.mode import ConstraintSet
    from ad_hoc_diffractometer.mode import SampleConstraint

    g = fourcv()
    g.wavelength = 1.54
    g.sample.lattice = ahd.Lattice(a=4.0)

    mode = ConstraintSet([SampleConstraint("chi", 90.0)])
    bad_sol = {"omega": 10.0, "chi": 45.0, "phi": 0.0, "ttheta": 20.0}

    with pytest.raises(ConstraintViolation, match=re.escape("SampleConstraint")):
        _validate_solutions([bad_sol], mode, g)


# ---------------------------------------------------------------------------
# Coverage: remaining mode.py hash/is_satisfied/property paths
# ---------------------------------------------------------------------------


def test_sample_constraint_hash_used_in_set():
    sc1 = SampleConstraint("chi", 90.0)
    sc2 = SampleConstraint("chi", 90.0)
    assert len({sc1, sc2}) == 1


def test_detector_constraint_is_satisfied():
    g = _psic_like()
    angles = {"mu": 0.0, "eta": 10.0, "chi": 90.0, "phi": 0.0, "nu": 0.0, "delta": 20.0}
    dc = DetectorConstraint("nu", 0.0)
    assert dc.is_satisfied(angles, g) is True
    dc2 = DetectorConstraint("nu", 5.0)
    assert dc2.is_satisfied(angles, g) is False


def test_detector_constraint_hash_used_in_set():
    dc1 = DetectorConstraint("nu", 0.0)
    dc2 = DetectorConstraint("nu", 0.0)
    assert len({dc1, dc2}) == 1


def test_constraint_set_bisect_constraint_property_none():
    """bisect_constraint returns None when no BisectConstraint present."""
    cs = ConstraintSet([SampleConstraint("chi", 90.0)])
    assert cs.bisect_constraint is None


def test_constraint_set_bisect_constraint_property_present():
    """bisect_constraint returns the BisectConstraint when present."""
    bc = BisectConstraint("omega", "ttheta")
    cs = ConstraintSet([bc])
    assert cs.bisect_constraint is bc


def test_constraint_set_reference_constraint_absent_returns_none():
    cs = ConstraintSet([SampleConstraint("chi", 90.0)])
    assert cs.reference_constraint is None


def test_constant_stages_returns_list():
    cs = ConstraintSet([BisectConstraint("eta", "delta"), SampleConstraint("mu", 0.0)])
    result = cs.constant_stages
    assert isinstance(result, list)
    assert "eta" in result
    assert "mu" in result


def test_bisect_stages_returns_tuple():
    bc = BisectConstraint("eta", "delta")
    cs = ConstraintSet([bc])
    result = cs.bisect_stages()
    assert result == ("eta", "delta")


def test_bisect_stages_none_no_bisect():
    cs = ConstraintSet([SampleConstraint("chi", 90.0)])
    assert cs.bisect_stages() is None


def test_is_fully_constrained_returns_bool():
    g = _fourcv_like()
    cs = ConstraintSet([BisectConstraint("omega", "ttheta")])
    assert cs.is_fully_constrained(g) is True


def test_apply_cut_point_no_cut():
    cs = ConstraintSet([SampleConstraint("chi", 90.0)])
    assert cs.apply_cut_point("chi", 350.0) == 350.0


def test_apply_cut_point_with_cut():
    cs = ConstraintSet([SampleConstraint("chi", 90.0)], cut_points={"chi": 0.0})
    assert cs.apply_cut_point("chi", -10.0) == pytest.approx(350.0)


def test_constraint_set_from_dict_unknown_type_raises():
    """from_dict raises ValueError for unknown constraint type."""
    d = {
        "type": "ConstraintSet",
        "constraints": [{"type": "UnknownConstraint", "name": "x", "value": 0.0}],
        "computed": None,
        "extras": {},
        "cut_points": {},
    }
    with pytest.raises(ValueError, match=re.escape("unknown constraint type")):
        ConstraintSet.from_dict(d)


def test_constraint_set_from_dict_required_optional():
    """from_dict restores REQUIRED and OPTIONAL sentinels in extras."""
    cs = ConstraintSet(
        [BisectConstraint("omega", "ttheta")],
        extras={"n_hat": REQUIRED, "psi": None, "opt": OPTIONAL},
    )
    d = cs.to_dict()
    cs2 = ConstraintSet.from_dict(d)
    assert cs2.extras["n_hat"] is REQUIRED
    assert cs2.extras["psi"] is None
    assert cs2.extras["opt"] is OPTIONAL


def test_constraint_set_from_dict_with_reference():
    """from_dict handles ReferenceConstraint entries."""
    cs = ConstraintSet(
        [
            BisectConstraint("omega", "ttheta"),
            ReferenceConstraint("psi", 90.0),
        ]
    )
    d = cs.to_dict()
    cs2 = ConstraintSet.from_dict(d)
    assert cs2.reference_constraint is not None
    assert cs2.reference_constraint.name == "psi"


def test_constraint_set_repr_contains_bisect():
    cs = ConstraintSet([BisectConstraint("eta", "delta")])
    r = repr(cs)
    assert "ConstraintSet" in r
    assert "eta" in r
    assert "delta" in r


def test_constraint_set_eq_different_cutpoints():
    cs1 = ConstraintSet([SampleConstraint("chi", 90.0)], cut_points={"chi": 0.0})
    cs2 = ConstraintSet([SampleConstraint("chi", 90.0)])
    assert cs1 != cs2


def test_constraint_set_eq_different_type():
    cs = ConstraintSet([SampleConstraint("chi", 90.0)])
    assert cs != "not a constraint set"


def test_constraint_set_eq_different_constraints():
    cs1 = ConstraintSet([SampleConstraint("chi", 90.0)])
    cs2 = ConstraintSet([SampleConstraint("phi", 0.0)])
    assert cs1 != cs2


def test_mode_dict_setitem_invalid_raises():
    md = ModeDict()
    with pytest.raises(TypeError, match=re.escape("ConstraintSet instances")):
        md["bad"] = "not a mode"  # type: ignore[assignment]


def test_mode_dict_delitem():
    md = ModeDict({"m": ConstraintSet([SampleConstraint("chi", 0.0)])})
    del md["m"]
    assert "m" not in md


def test_mode_dict_len():
    md = ModeDict({"a": ConstraintSet([SampleConstraint("chi", 0.0)])})
    assert len(md) == 1


def test_mode_dict_iter():
    keys = ["a", "b"]
    md = ModeDict({k: ConstraintSet([SampleConstraint("chi", 0.0)]) for k in keys})
    assert list(md) == keys


def test_mode_dict_repr():
    md = ModeDict({"m": ConstraintSet([SampleConstraint("chi", 0.0)])})
    assert "ModeDict" in repr(md)


def test_mode_dict_eq_different():
    m1 = ModeDict({"m": ConstraintSet([SampleConstraint("chi", 0.0)])})
    m2 = ModeDict({"m": ConstraintSet([SampleConstraint("chi", 90.0)])})
    assert m1 != m2


def test_mode_dict_eq_non_mode_dict():
    md = ModeDict()
    assert md != {}


def test_mode_dict_values():
    cs = ConstraintSet([SampleConstraint("chi", 0.0)])
    md = ModeDict({"m": cs})
    assert list(md.values()) == [cs]
