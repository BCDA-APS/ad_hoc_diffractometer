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
from ad_hoc_diffractometer import fivec
from ad_hoc_diffractometer import fourch
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer import kappa4ch
from ad_hoc_diffractometer import kappa4cv
from ad_hoc_diffractometer import kappa6c
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
    "name, value",
    [
        pytest.param("chi", 90.0, id="fixed-chi-90"),
        pytest.param("phi", 0.0, id="fixed-phi-0"),
        pytest.param("mu", -45.0, id="fixed-mu-neg45"),
    ],
)
def test_sample_constraint_construction(name, value):
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


@pytest.mark.parametrize(
    "lhs, rhs, expected",
    [
        pytest.param(
            SampleConstraint("chi", 90.0),
            SampleConstraint("chi", 90.0),
            True,
            id="sample-eq-same",
        ),
        pytest.param(
            SampleConstraint("chi", 90.0),
            SampleConstraint("phi", 90.0),
            False,
            id="sample-eq-different-name",
        ),
        pytest.param(
            SampleConstraint("chi", 90.0),
            SampleConstraint("chi", 0.0),
            False,
            id="sample-eq-different-value",
        ),
        pytest.param(
            SampleConstraint("chi", 90.0),
            DetectorConstraint("chi", 90.0),
            False,
            id="sample-eq-different-type",
        ),
        pytest.param(
            BisectConstraint("omega", "ttheta"),
            BisectConstraint("omega", "ttheta"),
            True,
            id="bisect-eq-same",
        ),
        pytest.param(
            BisectConstraint("omega", "ttheta"),
            BisectConstraint("eta", "delta"),
            False,
            id="bisect-eq-different",
        ),
        pytest.param(
            BisectConstraint("omega", "ttheta"),
            SampleConstraint("omega", 0.0),
            False,
            id="bisect-eq-different-type",
        ),
        pytest.param(
            DetectorConstraint("nu", 0.0),
            DetectorConstraint("nu", 0.0),
            True,
            id="detector-eq-same",
        ),
        pytest.param(
            DetectorConstraint("nu", 0.0),
            DetectorConstraint("delta", 0.0),
            False,
            id="detector-eq-different",
        ),
        pytest.param(
            DetectorConstraint("nu", 0.0),
            SampleConstraint("nu", 0.0),
            False,
            id="detector-eq-different-type",
        ),
        pytest.param(
            ReferenceConstraint("psi", 90.0),
            ReferenceConstraint("psi", 90.0),
            True,
            id="reference-eq-same",
        ),
        pytest.param(
            ReferenceConstraint("psi", 90.0),
            ReferenceConstraint("psi", 0.0),
            False,
            id="reference-eq-different",
        ),
        pytest.param(
            ReferenceConstraint("psi", 90.0),
            SampleConstraint("psi", 90.0),
            False,
            id="reference-eq-different-type",
        ),
    ],
)
def test_constraint_eq(lhs, rhs, expected):
    assert (lhs == rhs) is expected


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


@pytest.mark.parametrize(
    "constraint_factory, args, geometry_factory, expected",
    [
        pytest.param(
            SampleConstraint,
            ("chi", 90.0),
            _fourcv_like,
            True,
            id="sample-real-stage-chi",
        ),
        pytest.param(
            SampleConstraint,
            ("phi", 0.0),
            _fourcv_like,
            True,
            id="sample-real-stage-phi",
        ),
        pytest.param(
            SampleConstraint,
            ("mu", 0.0),
            _fourcv_like,
            False,
            id="sample-missing-stage-mu",
        ),
        pytest.param(
            BisectConstraint,
            ("omega", "ttheta"),
            _fourcv_like,
            True,
            id="bisect-both-stages-exist-fourcv",
        ),
        pytest.param(
            BisectConstraint,
            ("eta", "delta"),
            _psic_like,
            True,
            id="bisect-both-stages-exist-psic",
        ),
        pytest.param(
            BisectConstraint,
            ("eta", "ttheta"),
            _fourcv_like,
            False,
            id="bisect-sample-stage-missing",
        ),
        pytest.param(
            BisectConstraint,
            ("omega", "delta"),
            _fourcv_like,
            False,
            id="bisect-detector-stage-missing",
        ),
        pytest.param(
            DetectorConstraint,
            ("nu", 0.0),
            _psic_like,
            True,
            id="detector-real-stage-nu",
        ),
        pytest.param(
            DetectorConstraint,
            ("delta", 0.0),
            _psic_like,
            True,
            id="detector-real-stage-delta",
        ),
        pytest.param(
            DetectorConstraint,
            ("gamma", 0.0),
            _psic_like,
            False,
            id="detector-missing-stage-gamma",
        ),
        pytest.param(
            DetectorConstraint,
            ("qaz", 90.0),
            _psic_like,
            False,
            id="detector-qaz-not-implemented",
        ),
    ],
)
def test_constraint_is_implemented(
    constraint_factory, args, geometry_factory, expected
):
    g = geometry_factory()
    assert constraint_factory(*args).is_implemented(g) is expected


_ANGLES_FOURCV = {"omega": 10.0, "chi": 90.0, "phi": 0.0, "ttheta": 20.0}
_ANGLES_PSIC = {
    "mu": 0.0,
    "eta": 10.0,
    "chi": 90.0,
    "phi": 0.0,
    "nu": 0.0,
    "delta": 20.0,
}


@pytest.mark.parametrize(
    "constraint_factory, args, angles, expected_residual",
    [
        pytest.param(
            SampleConstraint,
            ("chi", 90.0),
            _ANGLES_FOURCV,
            0.0,
            id="sample-chi-satisfied",
        ),
        pytest.param(
            SampleConstraint,
            ("chi", 45.0),
            _ANGLES_FOURCV,
            45.0,
            id="sample-chi-not-satisfied",
        ),
        pytest.param(
            BisectConstraint,
            ("omega", "ttheta"),
            _ANGLES_FOURCV,
            0.0,
            id="bisect-satisfied",
        ),
        pytest.param(
            BisectConstraint,
            ("omega", "ttheta"),
            {**_ANGLES_FOURCV, "omega": 5.0},
            -5.0,
            id="bisect-not-satisfied",
        ),
        pytest.param(
            DetectorConstraint,
            ("nu", 0.0),
            _ANGLES_PSIC,
            0.0,
            id="detector-nu-satisfied",
        ),
        pytest.param(
            DetectorConstraint,
            ("nu", 5.0),
            _ANGLES_PSIC,
            -5.0,
            id="detector-nu-not-satisfied",
        ),
    ],
)
def test_constraint_evaluate(constraint_factory, args, angles, expected_residual):
    g = _fourcv_like()
    assert constraint_factory(*args).evaluate(angles, g) == pytest.approx(
        expected_residual
    )


@pytest.mark.parametrize(
    "constraint_factory, args, angles, expected_satisfied",
    [
        pytest.param(
            SampleConstraint,
            ("chi", 90.0),
            _ANGLES_FOURCV,
            True,
            id="sample-chi-satisfied",
        ),
        pytest.param(
            SampleConstraint,
            ("chi", 45.0),
            _ANGLES_FOURCV,
            False,
            id="sample-chi-not-satisfied",
        ),
        pytest.param(
            BisectConstraint,
            ("omega", "ttheta"),
            _ANGLES_FOURCV,
            True,
            id="bisect-satisfied",
        ),
        pytest.param(
            BisectConstraint,
            ("omega", "ttheta"),
            {**_ANGLES_FOURCV, "omega": 5.0},
            False,
            id="bisect-not-satisfied",
        ),
        pytest.param(
            DetectorConstraint,
            ("nu", 0.0),
            _ANGLES_PSIC,
            True,
            id="detector-nu-satisfied",
        ),
        pytest.param(
            DetectorConstraint,
            ("nu", 5.0),
            _ANGLES_PSIC,
            False,
            id="detector-nu-not-satisfied",
        ),
    ],
)
def test_constraint_is_satisfied(constraint_factory, args, angles, expected_satisfied):
    g = _fourcv_like()
    assert constraint_factory(*args).is_satisfied(angles, g) is expected_satisfied


def test_bisect_constraint_repr():
    bc = BisectConstraint("eta", "delta")
    r = repr(bc)
    assert "BisectConstraint" in r
    assert "eta" in r
    assert "delta" in r


@pytest.mark.parametrize(
    "c1, c2",
    [
        pytest.param(
            BisectConstraint("omega", "ttheta"),
            BisectConstraint("omega", "ttheta"),
            id="bisect-hash",
        ),
        pytest.param(
            SampleConstraint("chi", 90.0),
            SampleConstraint("chi", 90.0),
            id="sample-hash",
        ),
        pytest.param(
            DetectorConstraint("nu", 0.0),
            DetectorConstraint("nu", 0.0),
            id="detector-hash",
        ),
        pytest.param(
            ReferenceConstraint("psi", 90.0),
            ReferenceConstraint("psi", 90.0),
            id="reference-hash",
        ),
    ],
)
def test_constraint_hash(c1, c2):
    assert hash(c1) == hash(c2)
    assert len({c1, c2}) == 1


def test_bisect_constraint_category():
    assert BisectConstraint("omega", "ttheta").category == "sample"


def test_bisect_constraint_is_bisect_flag():
    assert BisectConstraint("omega", "ttheta").is_bisect is True


def test_bisect_constraint_name():
    assert BisectConstraint("omega", "ttheta").name == "bisect"


# ---------------------------------------------------------------------------
# DetectorConstraint
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, value, expected_is_qaz",
    [
        pytest.param("delta", 0.0, False, id="fixed-delta"),
        pytest.param("nu", 0.0, False, id="fixed-nu"),
        pytest.param("ttheta", 20.0, False, id="fixed-ttheta"),
        pytest.param("qaz", 90.0, True, id="qaz-pseudo"),
    ],
)
def test_detector_constraint_construction(name, value, expected_is_qaz):
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


def test_detector_constraint_to_dict_from_dict():
    dc = DetectorConstraint("nu", 0.0)
    d = dc.to_dict()
    assert d == {"type": "DetectorConstraint", "name": "nu", "value": 0.0}
    dc2 = DetectorConstraint.from_dict(d)
    assert dc2 == dc


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


@pytest.mark.parametrize(
    "constraints, expected_in, expected_not_in",
    [
        pytest.param(
            [BisectConstraint("omega", "ttheta")],
            ["omega"],
            [],
            id="bisect-only",
        ),
        pytest.param(
            [
                BisectConstraint("eta", "delta"),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            ["eta", "mu", "nu"],
            [],
            id="psic-bisecting",
        ),
        pytest.param(
            [SampleConstraint("chi", 90.0), DetectorConstraint("nu", 0.0)],
            ["chi", "nu"],
            [],
            id="fixed-sample-and-detector",
        ),
        pytest.param(
            [BisectConstraint("omega", "ttheta"), ReferenceConstraint("psi", 90.0)],
            ["omega"],
            ["psi"],
            id="reference-excluded",
        ),
        pytest.param(
            [SampleConstraint("chi", 90.0), DetectorConstraint("qaz", 90.0)],
            ["chi"],
            ["qaz"],
            id="qaz-excluded",
        ),
    ],
)
def test_constant_stages(constraints, expected_in, expected_not_in):
    cs = ConstraintSet(constraints)
    stages = cs.constant_stages
    for name in expected_in:
        assert name in stages
    for name in expected_not_in:
        assert name not in stages


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


@pytest.mark.parametrize(
    "lhs, rhs, expected",
    [
        pytest.param(
            ConstraintSet([SampleConstraint("chi", 90.0)]),
            ConstraintSet([SampleConstraint("chi", 90.0)]),
            True,
            id="cs-eq-same",
        ),
        pytest.param(
            ConstraintSet([SampleConstraint("chi", 90.0)], cut_points={"chi": 0.0}),
            ConstraintSet([SampleConstraint("chi", 90.0)]),
            False,
            id="cs-eq-different-cutpoints",
        ),
        pytest.param(
            ConstraintSet([SampleConstraint("chi", 90.0)]),
            "not a constraint set",
            False,
            id="cs-eq-different-type",
        ),
        pytest.param(
            ConstraintSet([SampleConstraint("chi", 90.0)]),
            ConstraintSet([SampleConstraint("phi", 0.0)]),
            False,
            id="cs-eq-different-constraints",
        ),
    ],
)
def test_constraint_set_eq(lhs, rhs, expected):
    assert (lhs == rhs) is expected


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


# ---------------------------------------------------------------------------
# AdHocDiffractometer — geometry.py lines 126, 132-133, 627-628, 647, 1857
# ---------------------------------------------------------------------------


def test_geometry_modes_accepts_mode_dict_directly():
    """Passing a ModeDict directly to AdHocDiffractometer.modes covers line 126."""
    cs = ConstraintSet([BisectConstraint("omega", "ttheta")])
    md = ModeDict({"bisecting": cs})
    g = _fourcv_like(modes=md, default_mode="bisecting")
    assert g.mode_name == "bisecting"
    assert isinstance(g.modes, ModeDict)


def test_geometry_default_mode_invalid_raises():
    """Invalid default_mode raises ValueError (lines 132-133)."""
    cs = ConstraintSet([BisectConstraint("omega", "ttheta")])
    with pytest.raises(ValueError, match=re.escape("default_mode")):
        _fourcv_like(modes={"bisecting": cs}, default_mode="nonexistent")


def test_geometry_mode_name_setter_invalid_raises():
    """Setting mode_name to unknown name raises ValueError (lines 627-628)."""
    cs = ConstraintSet([BisectConstraint("omega", "ttheta")])
    g = _fourcv_like(modes={"bisecting": cs}, default_mode="bisecting")
    with pytest.raises(ValueError, match=re.escape("not available")):
        g.mode_name = "nonexistent"


def test_geometry_mode_property_returns_none_when_no_mode():
    """mode property returns None when mode_name is None (line 647)."""
    cs = ConstraintSet([BisectConstraint("omega", "ttheta")])
    g = _fourcv_like(modes={"bisecting": cs})
    assert g.mode_name is None
    assert g.mode is None


def test_geometry_from_dict_old_bisecting_mode_raises():
    """_mode_from_dict raises ValueError for old BisectingMode format (line 1857)."""
    g = _psic_like()
    d = g.to_dict()
    d["modes"]["bisecting"] = {
        "type": "BisectingMode",
        "sample_stage": "eta",
        "detector_stage": "delta",
        "frozen_angles": {"mu": 0.0, "nu": 0.0},
        "cut_points": {},
    }
    d["mode_name"] = "bisecting"
    with pytest.raises(ValueError, match=re.escape("BisectingMode")):
        AdHocDiffractometer.from_dict(d)


# ---------------------------------------------------------------------------
# Issue #149 — fourcv and fourch mode structure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, expected_modes",
    [
        pytest.param(
            fourcv,
            {
                "bisecting",
                "fixed_chi",
                "fixed_phi",
                "constant_omega",
                "psi_constant",
                "double_diffraction",
            },
            id="fourcv-six-modes",
        ),
        pytest.param(
            fourch,
            {
                "bisecting",
                "fixed_chi",
                "fixed_phi",
                "constant_omega",
                "psi_constant",
                "double_diffraction",
            },
            id="fourch-six-modes",
        ),
    ],
)
def test_four_circle_factory_mode_names(factory, expected_modes):
    """Both fourcv and fourch expose exactly the 6 declared mode names."""
    g = factory()
    assert set(g.modes.keys()) == expected_modes


@pytest.mark.parametrize(
    "factory, mode_name, expected_type, expected_has_bisect",
    [
        pytest.param(
            fourcv, "bisecting", BisectConstraint, True, id="fourcv-bisecting"
        ),
        pytest.param(
            fourcv, "fixed_chi", SampleConstraint, False, id="fourcv-fixed_chi"
        ),
        pytest.param(
            fourcv, "fixed_phi", SampleConstraint, False, id="fourcv-fixed_phi"
        ),
        pytest.param(
            fourcv,
            "constant_omega",
            SampleConstraint,
            False,
            id="fourcv-constant_omega",
        ),
        pytest.param(
            fourcv, "psi_constant", ReferenceConstraint, False, id="fourcv-psi_constant"
        ),
        pytest.param(
            fourcv,
            "double_diffraction",
            BisectConstraint,
            True,
            id="fourcv-double_diffraction",
        ),
        pytest.param(
            fourch, "bisecting", BisectConstraint, True, id="fourch-bisecting"
        ),
        pytest.param(
            fourch, "fixed_chi", SampleConstraint, False, id="fourch-fixed_chi"
        ),
        pytest.param(
            fourch, "fixed_phi", SampleConstraint, False, id="fourch-fixed_phi"
        ),
        pytest.param(
            fourch,
            "constant_omega",
            SampleConstraint,
            False,
            id="fourch-constant_omega",
        ),
        pytest.param(
            fourch, "psi_constant", ReferenceConstraint, False, id="fourch-psi_constant"
        ),
        pytest.param(
            fourch,
            "double_diffraction",
            BisectConstraint,
            True,
            id="fourch-double_diffraction",
        ),
    ],
)
def test_four_circle_mode_constraint_type(
    factory, mode_name, expected_type, expected_has_bisect
):
    """Each mode's leading constraint is of the expected type."""
    g = factory()
    cs = g.modes[mode_name]
    assert isinstance(cs, ConstraintSet)
    assert cs.has_bisect == expected_has_bisect
    assert any(isinstance(c, expected_type) for c in cs.constraints)


@pytest.mark.parametrize(
    "factory, mode_name, expected_implemented",
    [
        pytest.param(fourcv, "bisecting", True, id="fourcv-bisecting-implemented"),
        pytest.param(fourcv, "fixed_chi", True, id="fourcv-fixed_chi-implemented"),
        pytest.param(fourcv, "fixed_phi", True, id="fourcv-fixed_phi-implemented"),
        pytest.param(
            fourcv, "constant_omega", True, id="fourcv-constant_omega-implemented"
        ),
        pytest.param(fourcv, "psi_constant", False, id="fourcv-psi_constant-stub"),
        pytest.param(
            fourcv,
            "double_diffraction",
            True,
            id="fourcv-double_diffraction-implemented",
        ),
        pytest.param(fourch, "bisecting", True, id="fourch-bisecting-implemented"),
        pytest.param(fourch, "fixed_chi", True, id="fourch-fixed_chi-implemented"),
        pytest.param(fourch, "fixed_phi", True, id="fourch-fixed_phi-implemented"),
        pytest.param(
            fourch, "constant_omega", True, id="fourch-constant_omega-implemented"
        ),
        pytest.param(fourch, "psi_constant", False, id="fourch-psi_constant-stub"),
        pytest.param(
            fourch,
            "double_diffraction",
            True,
            id="fourch-double_diffraction-implemented",
        ),
    ],
)
def test_four_circle_mode_is_implemented(factory, mode_name, expected_implemented):
    """Implemented modes return True; stubs return False from is_implemented()."""
    g = factory()
    assert g.modes[mode_name].is_implemented(g) == expected_implemented


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(fourcv, id="fourcv"),
        pytest.param(fourch, id="fourch"),
    ],
)
def test_four_circle_default_mode_is_bisecting(factory):
    """Default mode for both four-circle geometries is 'bisecting'."""
    g = factory()
    assert g.mode_name == "bisecting"


@pytest.mark.parametrize(
    "factory, mode_name, expected_computed",
    [
        pytest.param(
            fourcv,
            "bisecting",
            ["omega", "chi", "phi", "ttheta"],
            id="fourcv-bisecting-computed",
        ),
        pytest.param(
            fourcv,
            "fixed_chi",
            ["omega", "phi", "ttheta"],
            id="fourcv-fixed_chi-computed",
        ),
        pytest.param(
            fourcv,
            "fixed_phi",
            ["omega", "chi", "ttheta"],
            id="fourcv-fixed_phi-computed",
        ),
        pytest.param(
            fourcv,
            "constant_omega",
            ["chi", "phi", "ttheta"],
            id="fourcv-constant_omega-computed",
        ),
    ],
)
def test_four_circle_computed_stages_declared(factory, mode_name, expected_computed):
    """computed field lists the correct stage names for each mode."""
    g = factory()
    cs = g.modes[mode_name]
    assert cs.computed == expected_computed


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(fourcv, id="fourcv"),
        pytest.param(fourch, id="fourch"),
    ],
)
def test_four_circle_free_dof(factory):
    """Four-circle geometries have free_dof_after_bragg == 1."""
    g = factory()
    assert g.free_dof_after_bragg == 1


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(fourcv, id="fourcv"),
        pytest.param(fourch, id="fourch"),
    ],
)
def test_four_circle_modes_round_trip_serialisation(factory):
    """Full to_dict / from_dict round-trip preserves all 6 modes."""
    import json

    g = factory()
    d = g.to_dict()
    assert json.dumps(d)  # must be JSON-serialisable
    assert set(d["modes"].keys()) == {
        "bisecting",
        "fixed_chi",
        "fixed_phi",
        "constant_omega",
        "psi_constant",
        "double_diffraction",
    }
    g2 = AdHocDiffractometer.from_dict(d)
    assert set(g2.modes.keys()) == set(g.modes.keys())
    assert g2.mode_name == "bisecting"


# ---------------------------------------------------------------------------
# Issue #154 — fivec mode structure
# ---------------------------------------------------------------------------

_FIVEC_MODES = {
    "bisecting_4c",
    "fixed_chi",
    "fixed_phi",
    "fixed_mu",
    "constant_omega_noncoplanar",
}


def test_fivec_factory_mode_names():
    """fivec exposes exactly the 5 declared mode names."""
    assert set(fivec().modes.keys()) == _FIVEC_MODES


def test_fivec_default_mode_is_bisecting_4c():
    """Default mode for fivec is 'bisecting_4c'."""
    assert fivec().mode_name == "bisecting_4c"


def test_fivec_free_dof():
    """fivec has free_dof_after_bragg == 2."""
    assert fivec().free_dof_after_bragg == 2


@pytest.mark.parametrize(
    "mode_name",
    [pytest.param(m, id=m) for m in _FIVEC_MODES],
)
def test_fivec_mode_is_constraint_set(mode_name):
    """Every fivec mode is a ConstraintSet with exactly 2 constraints."""
    g = fivec()
    cs = g.modes[mode_name]
    assert isinstance(cs, ConstraintSet)
    assert len(cs) == 2


@pytest.mark.parametrize(
    "mode_name, expected_implemented",
    [
        pytest.param("bisecting_4c", True, id="bisecting_4c-implemented"),
        pytest.param("fixed_chi", True, id="fixed_chi-implemented"),
        pytest.param("fixed_phi", True, id="fixed_phi-implemented"),
        pytest.param("fixed_mu", True, id="fixed_mu-implemented"),
        pytest.param(
            "constant_omega_noncoplanar",
            True,
            id="constant_omega_noncoplanar-implemented",
        ),
    ],
)
def test_fivec_mode_is_implemented(mode_name, expected_implemented):
    """All fivec modes are implemented by the generic solver."""
    g = fivec()
    assert g.modes[mode_name].is_implemented(g) == expected_implemented


@pytest.mark.parametrize(
    "mode_name, expected_has_bisect",
    [
        pytest.param("bisecting_4c", True, id="bisecting_4c-has-bisect"),
        pytest.param("fixed_chi", False, id="fixed_chi-no-bisect"),
        pytest.param("fixed_phi", False, id="fixed_phi-no-bisect"),
        pytest.param("fixed_mu", True, id="fixed_mu-has-bisect"),
        pytest.param(
            "constant_omega_noncoplanar", False, id="constant_omega-no-bisect"
        ),
    ],
)
def test_fivec_mode_has_bisect(mode_name, expected_has_bisect):
    """BisectConstraint presence matches expected for each fivec mode."""
    g = fivec()
    assert g.modes[mode_name].has_bisect == expected_has_bisect


@pytest.mark.parametrize(
    "mode_name, expected_computed",
    [
        pytest.param(
            "bisecting_4c", ["omega", "chi", "phi", "ttheta"], id="bisecting_4c"
        ),
        pytest.param("fixed_chi", ["omega", "phi", "ttheta"], id="fixed_chi"),
        pytest.param("fixed_phi", ["omega", "chi", "ttheta"], id="fixed_phi"),
        pytest.param(
            "constant_omega_noncoplanar",
            ["mu", "chi", "phi", "ttheta"],
            id="constant_omega_noncoplanar",
        ),
    ],
)
def test_fivec_computed_stages(mode_name, expected_computed):
    """computed field lists the correct stage names for each fivec mode."""
    cs = fivec().modes[mode_name]
    assert cs.computed == expected_computed


def test_fivec_modes_round_trip_serialisation():
    """Full to_dict / from_dict round-trip preserves all 5 fivec modes."""
    import json

    g = fivec()
    d = g.to_dict()
    assert json.dumps(d)
    assert set(d["modes"].keys()) == _FIVEC_MODES
    g2 = AdHocDiffractometer.from_dict(d)
    assert set(g2.modes.keys()) == _FIVEC_MODES
    assert g2.mode_name == "bisecting_4c"


# ---------------------------------------------------------------------------
# Issue #155 — sixc mode structure
# ---------------------------------------------------------------------------

_SIXC_MODES = {
    "bisecting_4c",
    "fixed_gamma_5c",
    "fixed_alpha_5c",
    "fixed_alpha_zaxis",
    "fixed_beta_zaxis",
    "alpha_eq_beta_zaxis",
}

_SIXC_IMPLEMENTED = {"bisecting_4c", "fixed_gamma_5c", "fixed_alpha_5c"}
_SIXC_STUBS = {"fixed_alpha_zaxis", "fixed_beta_zaxis", "alpha_eq_beta_zaxis"}


def test_sixc_factory_mode_names():
    """sixc exposes exactly the 6 declared mode names."""
    assert set(sixc().modes.keys()) == _SIXC_MODES


def test_sixc_default_mode_is_bisecting_4c():
    """Default mode for sixc is 'bisecting_4c'."""
    assert sixc().mode_name == "bisecting_4c"


def test_sixc_free_dof():
    """sixc has free_dof_after_bragg == 3."""
    assert sixc().free_dof_after_bragg == 3


@pytest.mark.parametrize(
    "mode_name",
    [pytest.param(m, id=m) for m in _SIXC_MODES],
)
def test_sixc_mode_is_constraint_set(mode_name):
    """Every sixc mode is a ConstraintSet with exactly 3 constraints."""
    g = sixc()
    cs = g.modes[mode_name]
    assert isinstance(cs, ConstraintSet)
    assert len(cs) == 3


@pytest.mark.parametrize(
    "mode_name, expected_implemented",
    [pytest.param(m, True, id=f"{m}-implemented") for m in _SIXC_IMPLEMENTED]
    + [pytest.param(m, False, id=f"{m}-stub") for m in _SIXC_STUBS],
)
def test_sixc_mode_is_implemented(mode_name, expected_implemented):
    """Implemented modes return True; zaxis stubs return False."""
    g = sixc()
    assert g.modes[mode_name].is_implemented(g) == expected_implemented


@pytest.mark.parametrize(
    "mode_name, expected_has_bisect",
    [
        pytest.param("bisecting_4c", True, id="bisecting_4c-has-bisect"),
        pytest.param("fixed_gamma_5c", True, id="fixed_gamma_5c-has-bisect"),
        pytest.param("fixed_alpha_5c", True, id="fixed_alpha_5c-has-bisect"),
        pytest.param("fixed_alpha_zaxis", False, id="fixed_alpha_zaxis-no-bisect"),
        pytest.param("fixed_beta_zaxis", False, id="fixed_beta_zaxis-no-bisect"),
        pytest.param("alpha_eq_beta_zaxis", False, id="alpha_eq_beta_zaxis-no-bisect"),
    ],
)
def test_sixc_mode_has_bisect(mode_name, expected_has_bisect):
    """BisectConstraint presence matches expected for each sixc mode."""
    assert sixc().modes[mode_name].has_bisect == expected_has_bisect


@pytest.mark.parametrize(
    "mode_name, extras_key, expected_value",
    [
        pytest.param(
            "fixed_alpha_zaxis", "n_hat", "REQUIRED", id="fixed_alpha_zaxis-n_hat"
        ),
        pytest.param(
            "fixed_alpha_zaxis", "alpha_i", None, id="fixed_alpha_zaxis-alpha_i-output"
        ),
        pytest.param(
            "fixed_beta_zaxis", "n_hat", "REQUIRED", id="fixed_beta_zaxis-n_hat"
        ),
        pytest.param(
            "fixed_beta_zaxis", "beta_out", None, id="fixed_beta_zaxis-beta_out-output"
        ),
        pytest.param(
            "alpha_eq_beta_zaxis", "n_hat", "REQUIRED", id="alpha_eq_beta_zaxis-n_hat"
        ),
    ],
)
def test_sixc_zaxis_extras_declared(mode_name, extras_key, expected_value):
    """zaxis mode extras carry the expected sentinel or output placeholder."""
    g = sixc()
    mode = g.modes[mode_name]
    actual = mode.extras.get(extras_key)
    if expected_value == "REQUIRED":
        assert actual is REQUIRED
    else:
        assert actual is expected_value


def test_sixc_four_circle_computed_stages():
    """four_circle mode computed field lists omega, chi, phi, delta."""
    cs = sixc().modes["bisecting_4c"]
    assert cs.computed == ["omega", "chi", "phi", "delta"]


def test_sixc_modes_round_trip_serialisation():
    """Full to_dict / from_dict round-trip preserves all 6 sixc modes."""
    import json

    g = sixc()
    d = g.to_dict()
    assert json.dumps(d)
    assert set(d["modes"].keys()) == _SIXC_MODES
    g2 = AdHocDiffractometer.from_dict(d)
    assert set(g2.modes.keys()) == _SIXC_MODES
    assert g2.mode_name == "bisecting_4c"


# ---------------------------------------------------------------------------
# Issue #156 — zaxis and s2d2 mode structure
# ---------------------------------------------------------------------------

_ZAXIS_MODES = {"zaxis", "reflectivity"}
_S2D2_MODES = {"mu_fixed", "reflectivity"}


@pytest.mark.parametrize(
    "factory, expected_modes",
    [
        pytest.param(zaxis, _ZAXIS_MODES, id="zaxis"),
        pytest.param(s2d2, _S2D2_MODES, id="s2d2"),
    ],
)
def test_zaxis_s2d2_factory_mode_names(factory, expected_modes):
    """zaxis and s2d2 expose exactly their declared mode names."""
    assert set(factory().modes.keys()) == expected_modes


@pytest.mark.parametrize(
    "factory, expected_modes",
    [
        pytest.param(zaxis, _ZAXIS_MODES, id="zaxis"),
        pytest.param(s2d2, _S2D2_MODES, id="s2d2"),
    ],
)
def test_zaxis_s2d2_free_dof(factory, expected_modes):
    """Both geometries have free_dof_after_bragg == 1."""
    assert factory().free_dof_after_bragg == 1


@pytest.mark.parametrize(
    "factory, mode_name",
    [
        pytest.param(zaxis, "zaxis", id="zaxis-zaxis"),
        pytest.param(zaxis, "reflectivity", id="zaxis-reflectivity"),
        pytest.param(s2d2, "reflectivity", id="s2d2-reflectivity"),
    ],
)
def test_zaxis_s2d2_reference_modes_are_stubs(factory, mode_name):
    """All reference-constraint modes return is_implemented=False."""
    g = factory()
    cs = g.modes[mode_name]
    assert isinstance(cs, ConstraintSet)
    assert len(cs) == 1
    assert cs.reference_constraint is not None
    assert cs.is_implemented(g) is False


def test_s2d2_mu_fixed_is_implemented():
    """s2d2 mu_fixed uses a SampleConstraint — is_implemented returns True."""
    g = s2d2()
    cs = g.modes["mu_fixed"]
    assert isinstance(cs, ConstraintSet)
    assert len(cs) == 1
    assert cs.is_implemented(g) is True


@pytest.mark.parametrize(
    "factory, mode_name, extras_key, expected_value",
    [
        pytest.param(zaxis, "zaxis", "n_hat", "REQUIRED", id="zaxis-n_hat"),
        pytest.param(zaxis, "zaxis", "alpha_i", None, id="zaxis-alpha_i-output"),
        pytest.param(zaxis, "zaxis", "beta_out", None, id="zaxis-beta_out-output"),
        pytest.param(zaxis, "reflectivity", "n_hat", "REQUIRED", id="zaxis-refl-n_hat"),
        pytest.param(s2d2, "reflectivity", "n_hat", "REQUIRED", id="s2d2-refl-n_hat"),
        pytest.param(s2d2, "reflectivity", "alpha_i", None, id="s2d2-alpha_i-output"),
    ],
)
def test_zaxis_s2d2_extras_declared(factory, mode_name, extras_key, expected_value):
    """Reference modes carry expected REQUIRED sentinels and output placeholders."""
    g = factory()
    mode = g.modes[mode_name]
    actual = mode.extras.get(extras_key)
    if expected_value == "REQUIRED":
        assert actual is REQUIRED
    else:
        assert actual is expected_value


@pytest.mark.parametrize(
    "factory, mode_name, expected_computed",
    [
        pytest.param(zaxis, "zaxis", ["Z", "delta", "gamma"], id="zaxis-computed"),
        pytest.param(
            zaxis,
            "reflectivity",
            ["Z", "delta", "alpha", "gamma"],
            id="zaxis-refl-computed",
        ),
        pytest.param(
            s2d2, "mu_fixed", ["Z", "nu", "delta"], id="s2d2-mu_fixed-computed"
        ),
        pytest.param(
            s2d2, "reflectivity", ["mu", "Z", "nu", "delta"], id="s2d2-refl-computed"
        ),
    ],
)
def test_zaxis_s2d2_computed_stages(factory, mode_name, expected_computed):
    """computed field lists the correct stage names."""
    cs = factory().modes[mode_name]
    assert cs.computed == expected_computed


@pytest.mark.parametrize(
    "factory, expected_modes",
    [
        pytest.param(zaxis, _ZAXIS_MODES, id="zaxis"),
        pytest.param(s2d2, _S2D2_MODES, id="s2d2"),
    ],
)
def test_zaxis_s2d2_modes_round_trip(factory, expected_modes):
    """Full to_dict / from_dict round-trip preserves all modes."""
    import json

    g = factory()
    d = g.to_dict()
    assert json.dumps(d)
    assert set(d["modes"].keys()) == expected_modes
    g2 = AdHocDiffractometer.from_dict(d)
    assert set(g2.modes.keys()) == expected_modes
    assert g2.mode_name is None


# ---------------------------------------------------------------------------
# Issue #151 — kappa4cv and kappa4ch mode structure
# ---------------------------------------------------------------------------

_KAPPA4_MODES = {
    "bisecting",
    "fixed_kphi",
    "constant_omega",
    "constant_chi",
    "constant_phi",
    "psi_constant",
}

_KAPPA4_IMPLEMENTED = {
    "bisecting",
    "fixed_kphi",
    "constant_omega",
    "constant_chi",
    "constant_phi",
}
_KAPPA4_STUBS = _KAPPA4_MODES - _KAPPA4_IMPLEMENTED


@pytest.mark.parametrize(
    "factory",
    [pytest.param(kappa4cv, id="kappa4cv"), pytest.param(kappa4ch, id="kappa4ch")],
)
def test_kappa4_factory_mode_names(factory):
    """kappa4cv and kappa4ch expose exactly the 6 declared mode names."""
    assert set(factory().modes.keys()) == _KAPPA4_MODES


@pytest.mark.parametrize(
    "factory",
    [pytest.param(kappa4cv, id="kappa4cv"), pytest.param(kappa4ch, id="kappa4ch")],
)
def test_kappa4_free_dof(factory):
    """kappa4cv and kappa4ch have free_dof_after_bragg == 1."""
    assert factory().free_dof_after_bragg == 1


@pytest.mark.parametrize(
    "factory",
    [pytest.param(kappa4cv, id="kappa4cv"), pytest.param(kappa4ch, id="kappa4ch")],
)
def test_kappa4_default_mode(factory):
    """Default mode for kappa4cv and kappa4ch is 'bisecting'."""
    assert factory().mode_name == "bisecting"


@pytest.mark.parametrize(
    "factory, mode_name, expected_implemented",
    [
        pytest.param(kappa4cv, m, True, id=f"kappa4cv-{m}-impl")
        for m in sorted(_KAPPA4_IMPLEMENTED)
    ]
    + [
        pytest.param(kappa4cv, m, False, id=f"kappa4cv-{m}-stub")
        for m in sorted(_KAPPA4_STUBS)
    ]
    + [
        pytest.param(kappa4ch, m, True, id=f"kappa4ch-{m}-impl")
        for m in sorted(_KAPPA4_IMPLEMENTED)
    ]
    + [
        pytest.param(kappa4ch, m, False, id=f"kappa4ch-{m}-stub")
        for m in sorted(_KAPPA4_STUBS)
    ],
)
def test_kappa4_mode_is_implemented(factory, mode_name, expected_implemented):
    """Implemented modes return True; virtual-angle stubs return False."""
    g = factory()
    assert g.modes[mode_name].is_implemented(g) == expected_implemented


@pytest.mark.parametrize(
    "factory, mode_name, expected_has_bisect",
    [
        pytest.param(kappa4cv, "bisecting", True, id="kappa4cv-bisecting-has-bisect"),
        pytest.param(kappa4cv, "fixed_kphi", False, id="kappa4cv-fixed_kphi-no-bisect"),
        pytest.param(
            kappa4cv, "constant_omega", False, id="kappa4cv-constant_omega-no-bisect"
        ),
        pytest.param(
            kappa4cv, "psi_constant", False, id="kappa4cv-psi_constant-no-bisect"
        ),
    ],
)
def test_kappa4_mode_has_bisect(factory, mode_name, expected_has_bisect):
    """BisectConstraint presence matches expected for kappa4cv modes."""
    assert factory().modes[mode_name].has_bisect == expected_has_bisect


@pytest.mark.parametrize(
    "mode_name, extras_key, expected_value",
    [
        pytest.param("psi_constant", "n_hat", "REQUIRED", id="psi_constant-n_hat"),
        pytest.param("psi_constant", "psi", None, id="psi_constant-psi-output"),
    ],
)
def test_kappa4_psi_constant_extras(mode_name, extras_key, expected_value):
    """psi_constant carries REQUIRED n_hat and None psi output extras."""
    for g in (kappa4cv(), kappa4ch()):
        actual = g.modes[mode_name].extras.get(extras_key)
        if expected_value == "REQUIRED":
            assert actual is REQUIRED
        else:
            assert actual is expected_value


@pytest.mark.parametrize(
    "factory, mode_name, expected_computed",
    [
        pytest.param(
            kappa4cv,
            "bisecting",
            ["komega", "kappa", "kphi", "ttheta"],
            id="kappa4cv-bisecting",
        ),
        pytest.param(
            kappa4cv,
            "fixed_kphi",
            ["komega", "kappa", "ttheta"],
            id="kappa4cv-fixed_kphi",
        ),
        pytest.param(
            kappa4cv,
            "constant_chi",
            ["komega", "kappa", "kphi", "ttheta"],
            id="kappa4cv-constant_chi",
        ),
    ],
)
def test_kappa4_computed_stages(factory, mode_name, expected_computed):
    """computed field lists the correct stage names."""
    assert factory().modes[mode_name].computed == expected_computed


@pytest.mark.parametrize(
    "factory",
    [pytest.param(kappa4cv, id="kappa4cv"), pytest.param(kappa4ch, id="kappa4ch")],
)
def test_kappa4_modes_round_trip(factory):
    """Full to_dict / from_dict round-trip preserves all 6 modes."""
    import json

    g = factory()
    d = g.to_dict()
    assert json.dumps(d)
    assert set(d["modes"].keys()) == _KAPPA4_MODES
    g2 = AdHocDiffractometer.from_dict(d)
    assert set(g2.modes.keys()) == _KAPPA4_MODES
    assert g2.mode_name == "bisecting"


# ---------------------------------------------------------------------------
# Issue #152 — kappa6c mode structure
# ---------------------------------------------------------------------------

_KAPPA6C_MODES = {
    "bisecting_vertical",
    "bisecting_horizontal",
    "fixed_kphi",
    "fixed_mu",
    "fixed_nu",
    "fixed_delta",
    "lifting_detector_mu",
    "lifting_detector_kphi",
    "psi_constant_vertical",
    "psi_constant_horizontal",
}

_KAPPA6C_IMPLEMENTED = {
    "bisecting_vertical",
    "bisecting_horizontal",
    "fixed_kphi",
    "fixed_mu",
    "fixed_nu",
    "fixed_delta",
}
_KAPPA6C_STUBS = (
    _KAPPA6C_MODES - _KAPPA6C_IMPLEMENTED
)  # psi_constant_*, lifting_detector_*


def test_kappa6c_factory_mode_names():
    """kappa6c exposes exactly the 7 declared mode names."""
    assert set(kappa6c().modes.keys()) == _KAPPA6C_MODES


def test_kappa6c_free_dof():
    """kappa6c has free_dof_after_bragg == 3."""
    assert kappa6c().free_dof_after_bragg == 3


def test_kappa6c_default_mode():
    """Default mode for kappa6c is 'bisecting_vertical'."""
    assert kappa6c().mode_name == "bisecting_vertical"


@pytest.mark.parametrize(
    "mode_name, expected_implemented",
    [pytest.param(m, True, id=f"{m}-impl") for m in sorted(_KAPPA6C_IMPLEMENTED)]
    + [pytest.param(m, False, id=f"{m}-stub") for m in sorted(_KAPPA6C_STUBS)],
)
def test_kappa6c_mode_is_implemented(mode_name, expected_implemented):
    """Implemented modes return True; stubs return False."""
    g = kappa6c()
    assert g.modes[mode_name].is_implemented(g) == expected_implemented


@pytest.mark.parametrize(
    "mode_name, expected_has_bisect",
    [
        pytest.param("bisecting_vertical", True, id="bisecting_vertical-has-bisect"),
        pytest.param(
            "bisecting_horizontal", True, id="bisecting_horizontal-has-bisect"
        ),
        pytest.param("fixed_kphi", False, id="fixed_kphi-no-bisect"),
        pytest.param("fixed_mu", True, id="fixed_mu-has-bisect"),
        pytest.param("fixed_nu", True, id="fixed_nu-has-bisect"),
        pytest.param("fixed_delta", True, id="fixed_delta-has-bisect"),
        pytest.param("lifting_detector_mu", False, id="lifting_detector_mu-no-bisect"),
        pytest.param(
            "psi_constant_vertical", True, id="psi_constant_vertical-has-bisect"
        ),
        pytest.param(
            "psi_constant_horizontal", True, id="psi_constant_horizontal-has-bisect"
        ),
    ],
)
def test_kappa6c_mode_has_bisect(mode_name, expected_has_bisect):
    """BisectConstraint presence matches expected."""
    assert kappa6c().modes[mode_name].has_bisect == expected_has_bisect


@pytest.mark.parametrize(
    "mode_name",
    [
        pytest.param("psi_constant_vertical", id="psi_constant_vertical"),
        pytest.param("psi_constant_horizontal", id="psi_constant_horizontal"),
    ],
)
def test_kappa6c_psi_constant_extras(mode_name):
    """Both psi_constant modes carry REQUIRED n_hat and None psi output."""
    cs = kappa6c().modes[mode_name]
    assert cs.extras.get("n_hat") is REQUIRED
    assert cs.extras.get("psi") is None


@pytest.mark.parametrize(
    "mode_name, expected_computed",
    [
        pytest.param(
            "bisecting_vertical",
            ["komega", "kappa", "kphi", "delta"],
            id="bisecting_vertical",
        ),
        pytest.param(
            "bisecting_horizontal",
            ["mu", "kappa", "kphi", "nu"],
            id="bisecting_horizontal",
        ),
        pytest.param("fixed_kphi", ["komega", "kappa", "delta"], id="fixed_kphi"),
        pytest.param("fixed_mu", ["komega", "kappa", "kphi", "delta"], id="fixed_mu"),
        pytest.param("fixed_nu", ["komega", "kappa", "kphi", "delta"], id="fixed_nu"),
        pytest.param("fixed_delta", ["mu", "kappa", "kphi", "nu"], id="fixed_delta"),
    ],
)
def test_kappa6c_computed_stages(mode_name, expected_computed):
    """computed field lists the correct stage names."""
    assert kappa6c().modes[mode_name].computed == expected_computed


def test_kappa6c_modes_round_trip():
    """Full to_dict / from_dict round-trip preserves all 7 kappa6c modes."""
    import json

    g = kappa6c()
    d = g.to_dict()
    assert json.dumps(d)
    assert set(d["modes"].keys()) == _KAPPA6C_MODES
    g2 = AdHocDiffractometer.from_dict(d)
    assert set(g2.modes.keys()) == _KAPPA6C_MODES
    assert g2.mode_name == "bisecting_vertical"
