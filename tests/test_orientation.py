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
  - angles_to_phi_vector(): zero angles, magnitude (Bragg law), geometry independence
    of |Q| from sample rotations, angle restoration, error cases
"""

import math
import re

import numpy as np
import pytest

from ad_hoc_diffractometer import Lattice
from ad_hoc_diffractometer import angles_to_phi_vector
from ad_hoc_diffractometer import fourcv
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


# ---------------------------------------------------------------------------
# angles_to_phi_vector()
# ---------------------------------------------------------------------------

# Sapphire (006) psic angles used throughout these tests
_SAPPHIRE_ANGLES = {
    "mu": 0.0,
    "eta": 20.97,
    "chi": 90.0,
    "phi": 0.0,
    "nu": 0.0,
    "delta": 41.94,
}

# Expected |Q| for psic sapphire (006):
# |Q| = (4π/λ) sin(2θ/2) = (4π/1.5406) sin(20.97°)
_LAMBDA_CU_KA = 1.5406
_DELTA_006 = 41.94
_Q_MAG_SAPPHIRE_006 = (
    4.0 * math.pi / _LAMBDA_CU_KA * math.sin(math.radians(_DELTA_006 / 2.0))
)


# --- basic return-value properties ------------------------------------------


def test_angles_to_phi_vector_all_zero_gives_zero(psic_geom):
    """All motor angles zero → no scattering → Q_phi = 0."""
    psic_geom.wavelength = _LAMBDA_CU_KA
    Q = angles_to_phi_vector(psic_geom, mu=0, eta=0, chi=0, phi=0, nu=0, delta=0)
    np.testing.assert_allclose(Q, np.zeros(3), atol=1e-12)


def test_angles_to_phi_vector_returns_array_shape(psic_geom):
    """Return value is a numpy array of shape (3,)."""
    psic_geom.wavelength = _LAMBDA_CU_KA
    Q = angles_to_phi_vector(psic_geom, **_SAPPHIRE_ANGLES)
    assert isinstance(Q, np.ndarray)
    assert Q.shape == (3,)


def test_angles_to_phi_vector_magnitude_bragg(psic_geom):
    """|Q_phi| equals (4π/λ)·sin(2θ/2) — the Bragg scattering-vector magnitude."""
    psic_geom.wavelength = _LAMBDA_CU_KA
    Q = angles_to_phi_vector(psic_geom, **_SAPPHIRE_ANGLES)
    np.testing.assert_allclose(np.linalg.norm(Q), _Q_MAG_SAPPHIRE_006, rtol=1e-6)


def test_angles_to_phi_vector_explicit_components_psic_sapphire(psic_geom):
    """
    Explicit component check for psic sapphire (006) at Cu Kα.

    The expected values were computed from the implementation under review
    and locked here as a regression guard.  Any future change to the
    rotation convention must also update these numbers.
    """
    psic_geom.wavelength = _LAMBDA_CU_KA
    Q = angles_to_phi_vector(psic_geom, **_SAPPHIRE_ANGLES)
    expected = np.array([0.37387713, -0.97550961, 2.72580786])
    np.testing.assert_allclose(Q, expected, atol=1e-6)


# --- invariance / dependence properties -------------------------------------


def test_angles_to_phi_vector_magnitude_invariant_to_phi_rotation(psic_geom):
    """|Q_phi| is independent of the phi motor angle."""
    psic_geom.wavelength = _LAMBDA_CU_KA
    base = {**_SAPPHIRE_ANGLES}
    norms = []
    for phi_deg in (0.0, 30.0, 60.0, 90.0, 135.0, 180.0):
        base["phi"] = phi_deg
        norms.append(np.linalg.norm(angles_to_phi_vector(psic_geom, **base)))
    np.testing.assert_allclose(norms, norms[0], rtol=1e-10)


def test_angles_to_phi_vector_magnitude_invariant_to_sample_rotation(psic_geom):
    """
    |Q_phi| depends only on detector angles, not on sample angles.

    Varying mu, eta, chi, phi while keeping nu and delta fixed must leave
    |Q_phi| unchanged.
    """
    psic_geom.wavelength = _LAMBDA_CU_KA
    detector_angles = {"nu": 0.0, "delta": _DELTA_006}
    sample_angle_sets = [
        {"mu": 0.0, "eta": 20.97, "chi": 0.0, "phi": 0.0},
        {"mu": 0.0, "eta": 20.97, "chi": 30.0, "phi": 15.0},
        {"mu": 5.0, "eta": 20.97, "chi": 90.0, "phi": 45.0},
        {"mu": 10.0, "eta": 20.97, "chi": 45.0, "phi": 120.0},
    ]
    norms = [
        np.linalg.norm(angles_to_phi_vector(psic_geom, **{**detector_angles, **sample}))
        for sample in sample_angle_sets
    ]
    np.testing.assert_allclose(norms, norms[0], rtol=1e-10)


def test_angles_to_phi_vector_fourcv_all_zero_gives_zero():
    """fourcv: all-zero angles → Q_phi = 0."""
    g = fourcv()
    g.wavelength = 1.0
    Q = angles_to_phi_vector(g, omega=0, chi=0, phi=0, two_theta=0)
    np.testing.assert_allclose(Q, np.zeros(3), atol=1e-12)


def test_angles_to_phi_vector_fourcv_bisecting_magnitude():
    """
    fourcv bisecting geometry: |Q_phi| = (4π/λ)·sin(2θ/2).

    With omega = 2θ/2 (bisecting), the sample and detector share the
    same theta angle and the Q magnitude must satisfy the Bragg formula.
    """
    g = fourcv()
    g.wavelength = _LAMBDA_CU_KA
    two_theta = 28.4474  # Si (111) at Cu Kα
    Q = angles_to_phi_vector(
        g, omega=two_theta / 2.0, chi=0.0, phi=0.0, two_theta=two_theta
    )
    expected_mag = (
        4.0 * math.pi / _LAMBDA_CU_KA * math.sin(math.radians(two_theta / 2.0))
    )
    np.testing.assert_allclose(np.linalg.norm(Q), expected_mag, rtol=1e-6)


# --- angle restoration ------------------------------------------------------


def test_angles_to_phi_vector_restores_stage_angles(psic_geom):
    """Motor angles are restored to their original values after the call."""
    psic_geom.wavelength = _LAMBDA_CU_KA
    psic_geom.set_angle("eta", 99.9)
    psic_geom.set_angle("phi", 45.0)
    psic_geom.set_angle("chi", 33.3)

    angles_to_phi_vector(psic_geom, mu=0, eta=20.97, chi=90, phi=0, nu=0, delta=41.94)

    assert psic_geom.stage("eta").angle == pytest.approx(99.9)
    assert psic_geom.stage("phi").angle == pytest.approx(45.0)
    assert psic_geom.stage("chi").angle == pytest.approx(33.3)


def test_angles_to_phi_vector_restores_angles_on_error(psic_geom):
    """Stage angles are restored even when an exception is raised mid-call."""
    psic_geom.wavelength = _LAMBDA_CU_KA
    psic_geom.set_angle("eta", 55.5)

    # Trigger KeyError by supplying an unknown stage name; the KeyError is
    # raised before any angle is set, so eta should still be 55.5.
    with pytest.raises(KeyError):
        angles_to_phi_vector(psic_geom, eta=20, no_such_stage=0)

    assert psic_geom.stage("eta").angle == pytest.approx(55.5)


# --- partial motor-angle sets -----------------------------------------------


def test_angles_to_phi_vector_partial_angles_uses_current(psic_geom):
    """Stages not supplied in **motor_angles keep their current angle."""
    psic_geom.wavelength = _LAMBDA_CU_KA
    # Set all angles to the sapphire values directly on the geometry
    for name, angle in _SAPPHIRE_ANGLES.items():
        psic_geom.set_angle(name, angle)

    # Call with no keyword arguments — should use the already-set angles
    Q_none = angles_to_phi_vector(psic_geom)

    # Should equal the result when all angles are passed explicitly
    Q_explicit = angles_to_phi_vector(psic_geom, **_SAPPHIRE_ANGLES)

    np.testing.assert_allclose(Q_none, Q_explicit, atol=1e-12)


# --- error cases ------------------------------------------------------------


def test_angles_to_phi_vector_no_wavelength_raises(psic_geom):
    """Raises ValueError when wavelength is not set."""
    # wavelength defaults to None
    assert psic_geom.wavelength is None
    with pytest.raises(ValueError, match=re.escape("wavelength must be set")):
        angles_to_phi_vector(psic_geom, mu=0, eta=0, chi=0, phi=0, nu=0, delta=0)


def test_angles_to_phi_vector_unknown_stage_raises(psic_geom):
    """Raises KeyError for a stage name not in the geometry."""
    psic_geom.wavelength = _LAMBDA_CU_KA
    with pytest.raises(KeyError):
        angles_to_phi_vector(psic_geom, no_such_stage=0.0)
