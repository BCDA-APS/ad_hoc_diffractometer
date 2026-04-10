"""
Unit tests for ad_hoc_diffractometer.refinement.

Covers refine_lattice_bl1967 (#32) and refine_lattice_simplex (#33).

Both functions share the same interface and the same ``refine_all`` option:
  - refine_all=False (default): only the free parameters for the current
    crystal system are varied; symmetry constraints are maintained.
  - refine_all=True: all six cell parameters are treated as independent.

Tests
-----
- API: return dict keys, sample updated in-place, shape of arrays
- refine_all=False: correct number of free params per crystal system;
  constrained params stay equal (b=a for cubic, alpha=90 for orthorhombic, …)
- refine_all=True: all six parameters varied
- rms decreases from a perturbed UB (algorithmic correctness)
- Error cases: no parent, no UB, <3 reflections, bad refine flags, bad types
- _active_cell_params helper returns correct names per system
- _free_params_for_system helper
"""

import math
import re

import numpy as np
import pytest

from ad_hoc_diffractometer import Lattice
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer import psic
from ad_hoc_diffractometer import ub_from_three_reflections_bl1967
from ad_hoc_diffractometer import ub_from_two_reflections_bl1967
from ad_hoc_diffractometer.refinement import _active_cell_params
from ad_hoc_diffractometer.refinement import _free_params_for_system
from ad_hoc_diffractometer.refinement import refine_lattice_bl1967
from ad_hoc_diffractometer.refinement import refine_lattice_simplex
from ad_hoc_diffractometer.reflection import ReflectionList
from ad_hoc_diffractometer.sample import Sample
from ad_hoc_diffractometer.spec import g1_to_sample
from ad_hoc_diffractometer.spec import parse_fourc_g1

# ---------------------------------------------------------------------------
# Shared test geometry
# ---------------------------------------------------------------------------

_TWO_PI = 2.0 * math.pi

# psic with cubic a=1 (B=I), lambda=2*pi.  Three orthogonal reflections with
# exact motor angles.  UB from ub_from_three_reflections_bl1967 = I exactly.
_R1_ANG = {"mu": 0.0, "eta": 30.0, "chi": 0.0, "phi": 0.0, "nu": 0.0, "delta": 60.0}
_R2_ANG = {"mu": 0.0, "eta": 30.0, "chi": 0.0, "phi": 90.0, "nu": 0.0, "delta": 60.0}
_R3_ANG = {"mu": 0.0, "eta": 30.0, "chi": 90.0, "phi": 30.0, "nu": 0.0, "delta": 60.0}


@pytest.fixture
def cubic_geom():
    """
    psic geometry, cubic a=1 (B=I), lambda=2*pi.

    UB is set by ub_from_three_reflections_bl1967 (exact UB=I).
    """
    g = psic()
    g.wavelength = _TWO_PI
    g.sample.lattice = Lattice(a=1.0)
    g.add_reflection("r1", hkl=(1, 0, 0), angles=_R1_ANG)
    g.add_reflection("r2", hkl=(0, 1, 0), angles=_R2_ANG)
    g.add_reflection("r3", hkl=(0, 0, 1), angles=_R3_ANG)
    ub_from_three_reflections_bl1967(g.sample, "r1", "r2", "r3")
    return g


@pytest.fixture
def perturbed_geom(cubic_geom):
    """
    Same as cubic_geom but UB deliberately perturbed so initial rms > 0.
    Used to verify that the refinement actually reduces rms.
    """
    g = cubic_geom
    g.sample.UB = np.array([[1.1, 0.0, 0.0], [0.0, 0.9, 0.0], [0.0, 0.0, 1.0]])
    g.sample.U = g.sample.UB.copy()
    return g


# SPEC-based fourcv fixture (sapphire, real measured reflections)
_G1_LINE_A = (
    "#G1 4.785 4.785 12.991 90 90 120 "
    "1.516237713 1.516237713 0.483656786 90 90 60 "
    "0 0 6  1 0 0  "
    "41.94188 20.97 90 0  0 0  "
    "60 30 0 0  0 0  "
    "1.549802558 1.549802558  0 0"
)


@pytest.fixture
def sapphire_geom():
    """
    fourcv geometry with sapphire lattice from Align4Pete SPEC data.
    Three reflections: or1=(0,0,6), or2=(1,0,0), r3=(1,0,4).
    Initial UB set by ub_from_two_reflections_bl1967.
    """
    g = fourcv()
    g1_to_sample(parse_fourc_g1(_G1_LINE_A), g)
    g.add_reflection(
        "r3",
        hkl=(1, 0, 4),
        angles={
            "omega": 17.6428,
            "chi": 50.8925,
            "phi": 29.95,
            "two_theta": 35.392375,
        },
        wavelength=1.549802558,
    )
    ub_from_two_reflections_bl1967(g.sample)
    return g


# ---------------------------------------------------------------------------
# _free_params_for_system helper
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "system,expected",
    [
        pytest.param("cubic", ("a",), id="cubic"),
        pytest.param("tetragonal", ("a", "c"), id="tetragonal"),
        pytest.param("orthorhombic", ("a", "b", "c"), id="orthorhombic"),
        pytest.param("hexagonal", ("a", "c"), id="hexagonal"),
        pytest.param("trigonal", ("a", "alpha"), id="trigonal"),
        pytest.param("monoclinic", ("a", "b", "c", "beta"), id="monoclinic"),
        pytest.param(
            "triclinic",
            ("a", "b", "c", "alpha", "beta", "gamma"),
            id="triclinic",
        ),
    ],
)
def test_free_params_for_system(system, expected):
    assert _free_params_for_system(system) == expected


# ---------------------------------------------------------------------------
# _active_cell_params helper
# ---------------------------------------------------------------------------


def test_active_cell_params_cubic_refine_all_false():
    lat = Lattice(a=5.0)
    assert _active_cell_params(lat, refine_all=False) == ("a",)


def test_active_cell_params_cubic_refine_all_true():
    lat = Lattice(a=5.0)
    assert _active_cell_params(lat, refine_all=True) == (
        "a",
        "b",
        "c",
        "alpha",
        "beta",
        "gamma",
    )


def test_active_cell_params_hexagonal_refine_all_false():
    lat = Lattice(a=4.785, c=12.991, gamma=120)
    assert _active_cell_params(lat, refine_all=False) == ("a", "c")


def test_active_cell_params_monoclinic_refine_all_false():
    lat = Lattice(a=4.0, b=5.0, c=6.0, beta=95.0)
    assert _active_cell_params(lat, refine_all=False) == ("a", "b", "c", "beta")


# ---------------------------------------------------------------------------
# refine_lattice_bl1967 — API and result dict
# ---------------------------------------------------------------------------


class TestRefineBL1967API:
    def test_returns_dict(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        assert isinstance(result, dict)

    def test_result_keys(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        assert set(result.keys()) == {
            "lattice",
            "UB",
            "U",
            "residuals",
            "rms",
            "converged",
            "n_iter",
        }

    def test_result_lattice_is_Lattice(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        assert isinstance(result["lattice"], Lattice)

    def test_result_UB_shape(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        assert result["UB"].shape == (3, 3)

    def test_result_residuals_shape(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        assert result["residuals"].shape == (3, 3)

    def test_result_rms_is_float(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        assert isinstance(result["rms"], float)

    def test_result_converged_is_bool(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        assert isinstance(result["converged"], bool)

    def test_result_n_iter_is_int(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        assert isinstance(result["n_iter"], int)
        assert result["n_iter"] >= 1

    def test_sample_lattice_updated_in_place(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        assert perturbed_geom.sample.lattice is result["lattice"]

    def test_sample_UB_updated_in_place(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        np.testing.assert_array_equal(perturbed_geom.sample.UB, result["UB"])

    def test_sample_U_updated_in_place(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        np.testing.assert_array_equal(perturbed_geom.sample.U, result["U"])

    def test_accepts_reflection_objects(self, perturbed_geom):
        g = perturbed_geom
        r_objs = [g.sample.reflections[n] for n in ("r1", "r2", "r3")]
        result = refine_lattice_bl1967(g.sample, r_objs)
        assert result["UB"].shape == (3, 3)

    def test_accepts_mixed_str_and_objects(self, perturbed_geom):
        g = perturbed_geom
        r2_obj = g.sample.reflections["r2"]
        result = refine_lattice_bl1967(g.sample, ["r1", r2_obj, "r3"])
        assert result["UB"].shape == (3, 3)


# ---------------------------------------------------------------------------
# refine_lattice_bl1967 — rms reduction
# ---------------------------------------------------------------------------


class TestRefineBL1967Convergence:
    def test_rms_decreases_from_perturbed_UB(self, perturbed_geom):
        """Starting from a perturbed UB, refinement must reduce rms."""
        g = perturbed_geom
        # Compute initial rms manually
        UB0 = g.sample.UB.copy()
        from ad_hoc_diffractometer.orientation import angles_to_phi_vector

        rms0 = float(
            np.sqrt(
                np.mean(
                    [
                        np.sum(
                            (
                                angles_to_phi_vector(
                                    g, **g.sample.reflections[n].angles
                                )
                                - UB0 @ np.array(g.sample.reflections[n].hkl)
                            )
                            ** 2
                        )
                        for n in ("r1", "r2", "r3")
                    ]
                )
            )
        )
        result = refine_lattice_bl1967(g.sample, ["r1", "r2", "r3"])
        assert result["rms"] < rms0

    def test_rms_nonnegative(self, perturbed_geom):
        result = refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "r3"])
        assert result["rms"] >= 0.0


# ---------------------------------------------------------------------------
# refine_lattice_bl1967 — refine_all and crystal-system constraints
# ---------------------------------------------------------------------------


class TestRefineBL1967RefineAll:
    def test_refine_all_false_cubic_b_equals_a(self, perturbed_geom):
        """cubic: b must equal a after refinement (symmetry constraint)."""
        result = refine_lattice_bl1967(
            perturbed_geom.sample, ["r1", "r2", "r3"], refine_all=False
        )
        assert result["lattice"].b == pytest.approx(result["lattice"].a)

    def test_refine_all_false_cubic_alpha_is_90(self, perturbed_geom):
        """cubic: all angles must remain 90° after refinement."""
        result = refine_lattice_bl1967(
            perturbed_geom.sample, ["r1", "r2", "r3"], refine_all=False
        )
        assert result["lattice"].alpha == pytest.approx(90.0)
        assert result["lattice"].beta == pytest.approx(90.0)
        assert result["lattice"].gamma == pytest.approx(90.0)

    def test_refine_all_true_returns_six_independent(self, perturbed_geom):
        """refine_all=True does not crash and returns a valid lattice."""
        result = refine_lattice_bl1967(
            perturbed_geom.sample,
            ["r1", "r2", "r3"],
            refine_all=True,
        )
        assert isinstance(result["lattice"], Lattice)
        assert result["UB"].shape == (3, 3)

    def test_refine_all_false_hexagonal_gamma_is_120(self, sapphire_geom):
        """hexagonal: gamma must remain 120° when refine_all=False."""
        result = refine_lattice_bl1967(
            sapphire_geom.sample,
            ["or1", "or2", "r3"],
            refine_all=False,
        )
        assert result["lattice"].gamma == pytest.approx(120.0)

    def test_refine_all_false_hexagonal_b_equals_a(self, sapphire_geom):
        """hexagonal: b must equal a after refinement."""
        result = refine_lattice_bl1967(
            sapphire_geom.sample,
            ["or1", "or2", "r3"],
            refine_all=False,
        )
        assert result["lattice"].b == pytest.approx(result["lattice"].a)


# ---------------------------------------------------------------------------
# refine_lattice_bl1967 — error cases
# ---------------------------------------------------------------------------


class TestRefineBL1967Errors:
    def test_no_parent_raises(self):
        rl = ReflectionList(geometry_name="psic", valid_stages={"mu"})
        s = Sample(name="t", lattice=Lattice(a=1.0), reflections=rl)
        s.UB = np.eye(3)
        with pytest.raises(ValueError, match=re.escape("sample.parent")):
            refine_lattice_bl1967(s, [])

    def test_no_UB_raises(self, cubic_geom):
        cubic_geom.sample.UB = None
        with pytest.raises(ValueError, match=re.escape("sample.UB")):
            refine_lattice_bl1967(cubic_geom.sample, ["r1", "r2", "r3"])

    def test_fewer_than_3_reflections_raises(self, perturbed_geom):
        with pytest.raises(ValueError, match=re.escape("3 reflections")):
            refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2"])

    def test_neither_flag_raises(self, perturbed_geom):
        with pytest.raises(ValueError, match=re.escape("At least one")):
            refine_lattice_bl1967(
                perturbed_geom.sample,
                ["r1", "r2", "r3"],
                refine_cell=False,
                refine_orientation=False,
            )

    def test_bad_reflection_type_raises(self, perturbed_geom):
        with pytest.raises(TypeError, match=re.escape("Reflection object")):
            refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", 42])

    def test_unknown_reflection_name_raises(self, perturbed_geom):
        with pytest.raises(KeyError):
            refine_lattice_bl1967(perturbed_geom.sample, ["r1", "r2", "no_such"])


# ---------------------------------------------------------------------------
# refine_lattice_simplex — API and result dict (mirrors BL1967 tests)
# ---------------------------------------------------------------------------


class TestRefineSimplex:
    def test_returns_dict(self, perturbed_geom):
        result = refine_lattice_simplex(
            perturbed_geom.sample, ["r1", "r2", "r3"], max_iter=200
        )
        assert isinstance(result, dict)

    def test_result_keys(self, perturbed_geom):
        result = refine_lattice_simplex(
            perturbed_geom.sample, ["r1", "r2", "r3"], max_iter=200
        )
        assert set(result.keys()) == {
            "lattice",
            "UB",
            "U",
            "residuals",
            "rms",
            "converged",
            "n_iter",
        }

    def test_rms_nonnegative(self, perturbed_geom):
        result = refine_lattice_simplex(
            perturbed_geom.sample, ["r1", "r2", "r3"], max_iter=200
        )
        assert result["rms"] >= 0.0

    def test_UB_shape(self, perturbed_geom):
        result = refine_lattice_simplex(
            perturbed_geom.sample, ["r1", "r2", "r3"], max_iter=200
        )
        assert result["UB"].shape == (3, 3)

    def test_sample_updated_in_place(self, perturbed_geom):
        result = refine_lattice_simplex(
            perturbed_geom.sample, ["r1", "r2", "r3"], max_iter=200
        )
        assert perturbed_geom.sample.lattice is result["lattice"]
        np.testing.assert_array_equal(perturbed_geom.sample.UB, result["UB"])

    def test_refine_all_false_cubic_constraint(self, perturbed_geom):
        """Simplex with refine_all=False must maintain b=a for cubic."""
        result = refine_lattice_simplex(
            perturbed_geom.sample,
            ["r1", "r2", "r3"],
            refine_all=False,
            max_iter=500,
        )
        assert result["lattice"].b == pytest.approx(result["lattice"].a)

    def test_refine_all_false_cubic_angles_90(self, perturbed_geom):
        """Simplex: cubic angles remain 90° with refine_all=False."""
        result = refine_lattice_simplex(
            perturbed_geom.sample,
            ["r1", "r2", "r3"],
            refine_all=False,
            max_iter=500,
        )
        assert result["lattice"].alpha == pytest.approx(90.0)
        assert result["lattice"].beta == pytest.approx(90.0)
        assert result["lattice"].gamma == pytest.approx(90.0)

    def test_refine_all_false_hexagonal_b_equals_a(self, sapphire_geom):
        """Simplex: hexagonal b=a maintained with refine_all=False."""
        result = refine_lattice_simplex(
            sapphire_geom.sample,
            ["or1", "or2", "r3"],
            refine_all=False,
            max_iter=500,
        )
        assert result["lattice"].b == pytest.approx(result["lattice"].a)

    def test_refine_all_true_returns_valid(self, perturbed_geom):
        """Simplex refine_all=True does not crash."""
        result = refine_lattice_simplex(
            perturbed_geom.sample,
            ["r1", "r2", "r3"],
            refine_all=True,
            max_iter=200,
        )
        assert isinstance(result["lattice"], Lattice)

    def test_no_parent_raises(self):
        rl = ReflectionList(geometry_name="psic", valid_stages={"mu"})
        s = Sample(name="t", lattice=Lattice(a=1.0), reflections=rl)
        s.UB = np.eye(3)
        with pytest.raises(ValueError, match=re.escape("sample.parent")):
            refine_lattice_simplex(s, [])

    def test_no_UB_raises(self, cubic_geom):
        cubic_geom.sample.UB = None
        with pytest.raises(ValueError, match=re.escape("sample.UB")):
            refine_lattice_simplex(cubic_geom.sample, ["r1", "r2", "r3"])

    def test_fewer_than_3_reflections_raises(self, perturbed_geom):
        with pytest.raises(ValueError, match=re.escape("3 reflections")):
            refine_lattice_simplex(perturbed_geom.sample, ["r1", "r2"])

    def test_neither_flag_raises(self, perturbed_geom):
        with pytest.raises(ValueError, match=re.escape("At least one")):
            refine_lattice_simplex(
                perturbed_geom.sample,
                ["r1", "r2", "r3"],
                refine_cell=False,
                refine_orientation=False,
            )

    def test_bad_reflection_type_raises(self, perturbed_geom):
        with pytest.raises(TypeError, match=re.escape("Reflection object")):
            refine_lattice_simplex(perturbed_geom.sample, ["r1", "r2", 99])

    def test_unknown_reflection_name_raises(self, perturbed_geom):
        with pytest.raises(KeyError):
            refine_lattice_simplex(perturbed_geom.sample, ["r1", "r2", "no_such"])
