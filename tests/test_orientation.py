"""
Unit tests for ad_hoc_diffractometer.orientation.

Covers:
  - ub_identity()
  - ub_from_one_reflection(): Stage, str, None+parent, None+no-parent
  - Rodrigues edge cases: parallel, anti-parallel
  - U orthonormality
  - UB @ reference_hkl ≈ direction of reference_stage axis
  - sample.U and sample.UB updated in-place
  - reference_stage type handling (Stage, str, ndarray)
"""

import re

import numpy as np
import pytest

from ad_hoc_diffractometer import Lattice
from ad_hoc_diffractometer import psic
from ad_hoc_diffractometer import ub_from_one_reflection
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer.reflection import ReflectionList
from ad_hoc_diffractometer.sample import Sample

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_PSIC_ANGLES = {
    "mu": 0.0,
    "eta": 20.97,
    "chi": 90.0,
    "phi": 0.0,
    "nu": 0.0,
    "delta": 41.94,
}


@pytest.fixture
def sapphire_geom():
    """psic geometry with a sapphire sample and one (006) reflection."""
    g = psic()
    g.add_sample("sapphire", Lattice(a=4.758, c=12.991))
    g.sample = "sapphire"
    g.wavelength = 1.5406
    g.add_reflection(
        "r1",
        hkl=(0, 0, 6),
        angles=_PSIC_ANGLES,
    )
    g.sample.reflections.setor1("r1")
    return g


# ---------------------------------------------------------------------------
# ub_identity()
# ---------------------------------------------------------------------------


def test_ub_identity_sets_U_to_eye(psic_geom):
    ub_identity(psic_geom.sample)
    np.testing.assert_array_equal(psic_geom.sample.U, np.eye(3))


def test_ub_identity_sets_UB_to_B(psic_geom):
    ub_identity(psic_geom.sample)
    np.testing.assert_allclose(
        psic_geom.sample.UB, psic_geom.sample.lattice.B, atol=1e-12
    )


def test_ub_identity_returns_UB(psic_geom):
    UB = ub_identity(psic_geom.sample)
    np.testing.assert_array_equal(UB, psic_geom.sample.UB)


def test_ub_identity_updates_in_place(psic_geom):
    sample = psic_geom.sample
    UB = ub_identity(sample)
    assert sample.UB is UB


# ---------------------------------------------------------------------------
# ub_from_one_reflection() — reference_stage variants
# ---------------------------------------------------------------------------


def test_ub_one_refl_stage_object(sapphire_geom):
    """Pass a Stage object as reference_stage."""
    g = sapphire_geom
    UB = ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(0, 0, 1),
        reference_stage=g.stage("phi"),
    )
    assert g.sample.UB is UB
    assert g.sample.U is not None


def test_ub_one_refl_stage_string(sapphire_geom):
    """Pass a stage name string; resolved via sample.parent."""
    g = sapphire_geom
    UB = ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(0, 0, 1),
        reference_stage="phi",
    )
    assert UB.shape == (3, 3)


def test_ub_one_refl_none_with_parent_defaults_to_phi(sapphire_geom):
    """reference_stage=None with a parent geometry defaults to phi."""
    g = sapphire_geom
    UB_explicit = ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(0, 0, 1),
        reference_stage=g.stage("phi"),
    )
    # reset U/UB
    g.sample.U = None
    g.sample.UB = None
    UB_default = ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(0, 0, 1),
        reference_stage=None,
    )
    np.testing.assert_allclose(UB_explicit, UB_default, atol=1e-12)


def test_ub_one_refl_none_no_parent_raises():
    """reference_stage=None with no parent geometry raises ValueError."""
    from ad_hoc_diffractometer.reflection import ReflectionList

    rl = ReflectionList(
        geometry_name="psic", valid_stages={"mu", "eta", "chi", "phi", "nu", "delta"}
    )
    rl.add("r1", hkl=(0, 0, 6), angles={"mu": 0.0})
    sample = Sample(name="test", lattice=Lattice(a=1.0), reflections=rl)
    rl.setor1("r1")
    with pytest.raises(ValueError, match=re.escape("no parent geometry")):
        ub_from_one_reflection(sample, "r1", reference_stage=None)


def test_ub_one_refl_stage_string_no_parent_raises():
    """Stage string with no parent raises ValueError."""
    rl = ReflectionList(geometry_name="psic", valid_stages={"mu"})
    rl.add("r1", hkl=(0, 0, 1), angles={})
    sample = Sample(name="test", lattice=Lattice(a=1.0), reflections=rl)
    with pytest.raises(ValueError, match=re.escape("no parent geometry")):
        ub_from_one_reflection(sample, "r1", reference_stage="phi")


def test_ub_one_refl_ndarray_reference_stage(sapphire_geom):
    """Pass a raw ndarray as reference_stage."""
    g = sapphire_geom
    phi_axis = g.stage("phi").axis
    UB = ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(0, 0, 1),
        reference_stage=phi_axis,
    )
    assert UB.shape == (3, 3)


def test_ub_one_refl_reflection_object(sapphire_geom):
    """Pass a Reflection object directly."""
    g = sapphire_geom
    r = g.sample.reflections["r1"]
    UB = ub_from_one_reflection(
        g.sample,
        r,
        reference_hkl=(0, 0, 1),
        reference_stage=g.stage("phi"),
    )
    assert UB.shape == (3, 3)


def test_ub_one_refl_unknown_reflection_raises(sapphire_geom):
    g = sapphire_geom
    with pytest.raises(KeyError):
        ub_from_one_reflection(
            g.sample,
            "missing",
            reference_stage=g.stage("phi"),
        )


# ---------------------------------------------------------------------------
# ub_from_one_reflection() — mathematical correctness
# ---------------------------------------------------------------------------


def test_ub_one_refl_U_is_orthonormal(sapphire_geom):
    """U must be an orthonormal matrix (U.T @ U = I, det = 1)."""
    g = sapphire_geom
    ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(0, 0, 1),
        reference_stage=g.stage("phi"),
    )
    U = g.sample.U
    np.testing.assert_allclose(U.T @ U, np.eye(3), atol=1e-10)
    assert abs(np.linalg.det(U) - 1.0) < 1e-10


def test_ub_one_refl_UB_equals_U_at_B(sapphire_geom):
    """UB must equal U @ B."""
    g = sapphire_geom
    ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(0, 0, 1),
        reference_stage=g.stage("phi"),
    )
    np.testing.assert_allclose(g.sample.UB, g.sample.U @ g.sample.lattice.B, atol=1e-12)


def test_ub_one_refl_crystal_direction_maps_to_stage_axis(sapphire_geom):
    """UB @ reference_hkl should be parallel to reference_stage.axis."""
    g = sapphire_geom
    reference_hkl = np.array([0.0, 0.0, 1.0])
    phi_axis = g.stage("phi").axis
    ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=tuple(reference_hkl),
        reference_stage=g.stage("phi"),
    )
    # UB @ h must be parallel to phi_axis (same direction, any magnitude)
    q = g.sample.UB @ reference_hkl
    q_hat = q / np.linalg.norm(q)
    r_hat = phi_axis / np.linalg.norm(phi_axis)
    np.testing.assert_allclose(np.abs(np.dot(q_hat, r_hat)), 1.0, atol=1e-10)


def test_ub_one_refl_updates_sample_in_place(sapphire_geom):
    g = sapphire_geom
    assert g.sample.U is None
    assert g.sample.UB is None
    ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(0, 0, 1),
        reference_stage=g.stage("phi"),
    )
    assert g.sample.U is not None
    assert g.sample.UB is not None


# ---------------------------------------------------------------------------
# Rodrigues edge cases
# ---------------------------------------------------------------------------


def test_ub_one_refl_parallel_gives_identity_U(psic_geom):
    """When Bh_hat already points along r_hat, U should be identity."""
    g = psic_geom
    # cubic a=1: B = I; reference_hkl=(0,0,1) → Bh = (0,0,1)
    # phi axis of psic is -ZHAT = (0,0,-1) — not parallel
    # Use mu axis (+XHAT) and reference_hkl that maps to XHAT
    # B @ (1,0,0) = (1,0,0) = XHAT = mu.axis
    g.add_reflection(
        "r1",
        hkl=(1, 0, 0),
        angles={
            "mu": 0.0,
            "eta": 0.0,
            "chi": 0.0,
            "phi": 0.0,
            "nu": 0.0,
            "delta": 0.0,
        },
    )
    ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(1, 0, 0),
        reference_stage=g.stage("mu"),  # mu axis = +XHAT
    )
    np.testing.assert_allclose(g.sample.U, np.eye(3), atol=1e-10)


def test_ub_one_refl_antipara_U_is_rotation_matrix(psic_geom):
    """Anti-parallel case: U must still be a valid rotation matrix."""
    g = psic_geom
    # B @ (0,0,1) = XHAT (cubic a=1, c-axis maps to x in BL convention)
    # Actually default lattice is cubic a=1: B = I, so B@(0,0,1) = (0,0,1)
    # phi axis = -ZHAT = (0,0,-1) — anti-parallel to (0,0,1)
    g.add_reflection(
        "r1",
        hkl=(0, 0, 1),
        angles={
            "mu": 0.0,
            "eta": 0.0,
            "chi": 0.0,
            "phi": 0.0,
            "nu": 0.0,
            "delta": 0.0,
        },
    )
    ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(0, 0, 1),
        reference_stage=g.stage("phi"),  # phi = -ZHAT = (0,0,-1)
    )
    U = g.sample.U
    np.testing.assert_allclose(U.T @ U, np.eye(3), atol=1e-10)
    assert abs(np.linalg.det(U) - 1.0) < 1e-10


# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------


def test_ub_one_refl_zero_reference_hkl_raises(sapphire_geom):
    g = sapphire_geom
    with pytest.raises(ValueError, match=re.escape("zero vector")):
        ub_from_one_reflection(
            g.sample,
            "r1",
            reference_hkl=(0, 0, 0),
            reference_stage=g.stage("phi"),
        )


def test_ub_one_refl_bad_type_raises(sapphire_geom):
    g = sapphire_geom
    with pytest.raises(TypeError, match=re.escape("Reflection or a name string")):
        ub_from_one_reflection(
            g.sample,
            42,
            reference_stage=g.stage("phi"),
        )
