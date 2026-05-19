# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
Unit tests for the Lattice class in lattice.py.

Covers:
  - Crystal system deduction from minimum parameter sets (all 7 systems)
  - Default parameter filling for each system
  - Default lattice (cubic, a=1)
  - Lazy computation and cache invalidation
  - Parameter setters re-deduce system and invalidate cache
  - cartesian_lattice_vectors, reciprocal_lattice_vectors, B matrix correctness
  - __str__ reports only free parameters for the crystal system
  - __repr__ reports all six parameters
  - Parameter validation: strictly positive lengths, angle range (1,179),
    non-degenerate unit cell volume, system-specific constraints
  - Invalid inputs: zero/negative lengths, extreme angles, degenerate cells,
    ambiguous parameter combinations

Convention for context parameters:
  Each parametrize set includes a 'context' entry that is either
    does_not_raise()          -- the call must succeed
    pytest.raises(Exc, ...)   -- the call must raise Exc
  The match string is embedded directly in pytest.raises(..., match=...) and
  is not passed as a separate parameter.
"""

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import Lattice

# ---------------------------------------------------------------------------
# Crystal system deduction and default parameter filling
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, expected_system, expected_params, context",
    [
        # ------------------------------------------------------------------
        # Cubic
        # ------------------------------------------------------------------
        pytest.param(
            {"a": 5.0},
            "cubic",
            {"a": 5.0, "b": 5.0, "c": 5.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
            does_not_raise(),
            id="cubic-a-only",
        ),
        pytest.param(
            {"a": 3.0, "b": 3.0, "c": 3.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
            "cubic",
            {"a": 3.0, "b": 3.0, "c": 3.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
            does_not_raise(),
            id="cubic-all-six",
        ),
        # ------------------------------------------------------------------
        # Tetragonal
        # ------------------------------------------------------------------
        pytest.param(
            {"a": 3.0, "c": 6.0},
            "tetragonal",
            {"a": 3.0, "b": 3.0, "c": 6.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
            does_not_raise(),
            id="tetragonal-a-c",
        ),
        pytest.param(
            {"a": 3.0, "b": 3.0, "c": 6.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
            "tetragonal",
            {"a": 3.0, "b": 3.0, "c": 6.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
            does_not_raise(),
            id="tetragonal-all-six",
        ),
        # ------------------------------------------------------------------
        # Hexagonal
        # ------------------------------------------------------------------
        pytest.param(
            {"a": 4.785, "c": 12.991, "gamma": 120.0},
            "hexagonal",
            {
                "a": 4.785,
                "b": 4.785,
                "c": 12.991,
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 120.0,
            },
            does_not_raise(),
            id="hexagonal-a-c-gamma120-explicit",
        ),
        pytest.param(
            {
                "a": 4.785,
                "b": 4.785,
                "c": 12.991,
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 120.0,
            },
            "hexagonal",
            {
                "a": 4.785,
                "b": 4.785,
                "c": 12.991,
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 120.0,
            },
            does_not_raise(),
            id="hexagonal-all-six",
        ),
        # ------------------------------------------------------------------
        # Trigonal
        # ------------------------------------------------------------------
        pytest.param(
            {"a": 5.0, "alpha": 60.0},
            "trigonal",
            {"a": 5.0, "b": 5.0, "c": 5.0, "alpha": 60.0, "beta": 60.0, "gamma": 60.0},
            does_not_raise(),
            id="trigonal-a-alpha",
        ),
        pytest.param(
            {"a": 5.0, "b": 5.0, "c": 5.0, "alpha": 60.0, "beta": 60.0, "gamma": 60.0},
            "trigonal",
            {"a": 5.0, "b": 5.0, "c": 5.0, "alpha": 60.0, "beta": 60.0, "gamma": 60.0},
            does_not_raise(),
            id="trigonal-all-six",
        ),
        # ------------------------------------------------------------------
        # Orthorhombic
        # ------------------------------------------------------------------
        pytest.param(
            {"a": 2.0, "b": 3.0, "c": 4.0},
            "orthorhombic",
            {"a": 2.0, "b": 3.0, "c": 4.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
            does_not_raise(),
            id="orthorhombic-a-b-c",
        ),
        pytest.param(
            {"a": 2.0, "b": 3.0, "c": 4.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
            "orthorhombic",
            {"a": 2.0, "b": 3.0, "c": 4.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
            does_not_raise(),
            id="orthorhombic-all-six",
        ),
        # ------------------------------------------------------------------
        # Monoclinic
        # ------------------------------------------------------------------
        pytest.param(
            {"a": 5.0, "b": 6.0, "c": 7.0, "beta": 110.0},
            "monoclinic",
            {"a": 5.0, "b": 6.0, "c": 7.0, "alpha": 90.0, "beta": 110.0, "gamma": 90.0},
            does_not_raise(),
            id="monoclinic-a-b-c-beta",
        ),
        pytest.param(
            {"a": 5.0, "b": 6.0, "c": 7.0, "alpha": 90.0, "beta": 110.0, "gamma": 90.0},
            "monoclinic",
            {"a": 5.0, "b": 6.0, "c": 7.0, "alpha": 90.0, "beta": 110.0, "gamma": 90.0},
            does_not_raise(),
            id="monoclinic-all-six",
        ),
        # ------------------------------------------------------------------
        # Triclinic
        # ------------------------------------------------------------------
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 80.0, "beta": 85.0, "gamma": 95.0},
            "triclinic",
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 80.0, "beta": 85.0, "gamma": 95.0},
            does_not_raise(),
            id="triclinic-all-six",
        ),
    ],
)
def test_lattice_system_deduction(kwargs, expected_system, expected_params, context):
    with context:
        lat = Lattice(**kwargs)
        assert lat.system == expected_system
        for param, expected_val in expected_params.items():
            np.testing.assert_allclose(
                getattr(lat, param),
                expected_val,
                atol=1e-12,
                err_msg=f"Mismatch in {param}",
            )


# ---------------------------------------------------------------------------
# Default lattice
# ---------------------------------------------------------------------------


def test_lattice_default():
    lat = Lattice()
    assert lat.system == "cubic"
    assert lat.a == 1.0
    assert lat.b == 1.0
    assert lat.c == 1.0
    assert lat.alpha == 90.0
    assert lat.beta == 90.0
    assert lat.gamma == 90.0


# ---------------------------------------------------------------------------
# Invalid inputs — missing or insufficient parameters
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, context",
    [
        # Lattice(b=3.0): a defaults to 1.0 but b without c is ambiguous
        pytest.param(
            {"b": 3.0},
            pytest.raises(ValueError, match=re.escape("Cannot deduce crystal system")),
            id="invalid-b-without-a-or-c",
        ),
        # Lattice(alpha=70, beta=80): two non-90 angles, no system matches
        pytest.param(
            {"a": 3.0, "alpha": 70.0, "beta": 80.0},
            pytest.raises(ValueError, match=re.escape("Cannot deduce crystal system")),
            id="invalid-ambiguous-alpha-and-beta",
        ),
        # a + b + beta but no c
        pytest.param(
            {"a": 3.0, "b": 4.0, "beta": 100.0},
            pytest.raises(ValueError, match=re.escape("Cannot deduce crystal system")),
            id="invalid-a-b-beta-no-c",
        ),
        # Trigonal requires alpha != 90
        pytest.param(
            {"a": 3.0, "alpha": 90.0},
            pytest.raises(ValueError, match=re.escape("alpha must differ from 90")),
            id="invalid-trigonal-alpha-exactly-90",
        ),
        # Hexagonal requires gamma == 120
        pytest.param(
            {"a": 3.0, "c": 6.0, "gamma": 110.0},
            pytest.raises(ValueError, match=re.escape("gamma must be 120")),
            id="invalid-hexagonal-gamma-not-120",
        ),
    ],
)
def test_lattice_invalid_params(kwargs, context):
    with context:
        Lattice(**kwargs)


@pytest.mark.parametrize(
    "args, context",
    [
        pytest.param(
            ("_UNSET",) * 6,
            pytest.raises(ValueError, match=re.escape("'a' must always be provided")),
            id="no-a-supplied",
        ),
    ],
)
def test_deduce_system_and_params_errors(args, context):
    """_deduce_system_and_params raises for invalid direct inputs."""
    from ad_hoc_diffractometer.lattice import _UNSET
    from ad_hoc_diffractometer.lattice import _deduce_system_and_params

    real_args = tuple(_UNSET if a == "_UNSET" else a for a in args)
    with context:
        _deduce_system_and_params(*real_args)


# ---------------------------------------------------------------------------
# Invalid inputs — out-of-range values
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, context",
    [
        # Strictly positive lengths
        pytest.param(
            {"a": 0.0},
            pytest.raises(ValueError, match=re.escape("strictly positive")),
            id="invalid-a-zero",
        ),
        pytest.param(
            {"a": -1.0},
            pytest.raises(ValueError, match=re.escape("strictly positive")),
            id="invalid-a-negative",
        ),
        pytest.param(
            {"a": 3.0, "b": 0.0, "c": 4.0},
            pytest.raises(ValueError, match=re.escape("strictly positive")),
            id="invalid-b-zero",
        ),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": -2.0},
            pytest.raises(ValueError, match=re.escape("strictly positive")),
            id="invalid-c-negative",
        ),
        # Angle at or below minimum (1.0 degrees, open interval)
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 0.0, "beta": 90.0, "gamma": 90.0},
            pytest.raises(ValueError, match=re.escape("must be in")),
            id="invalid-alpha-zero",
        ),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 90.0, "beta": -10.0, "gamma": 90.0},
            pytest.raises(ValueError, match=re.escape("must be in")),
            id="invalid-beta-negative",
        ),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 90.0, "beta": 90.0, "gamma": 1.0},
            pytest.raises(ValueError, match=re.escape("must be in")),
            id="invalid-gamma-at-boundary-low",
        ),
        # Angle at or above maximum (179.0 degrees, open interval)
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 179.0, "beta": 90.0, "gamma": 90.0},
            pytest.raises(ValueError, match=re.escape("must be in")),
            id="invalid-alpha-at-boundary-high",
        ),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 90.0, "beta": 200.0, "gamma": 90.0},
            pytest.raises(ValueError, match=re.escape("must be in")),
            id="invalid-beta-exceeds-180",
        ),
        # Degenerate cell: volume factor <= 0
        pytest.param(
            {
                "a": 3.0,
                "b": 4.0,
                "c": 5.0,
                "alpha": 170.0,
                "beta": 170.0,
                "gamma": 170.0,
            },
            pytest.raises(ValueError, match=re.escape("degenerate")),
            id="invalid-degenerate-all-angles-170",
        ),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 5.0, "beta": 5.0, "gamma": 175.0},
            pytest.raises(ValueError, match=re.escape("degenerate")),
            id="invalid-degenerate-extreme-mix",
        ),
    ],
)
def test_lattice_invalid_values(kwargs, context):
    with context:
        Lattice(**kwargs)


# ---------------------------------------------------------------------------
# Invalid inputs via setters
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "initial_kwargs, setter, bad_value, context",
    [
        pytest.param(
            {"a": 5.0},
            "a",
            0.0,
            pytest.raises(ValueError, match=re.escape("strictly positive")),
            id="setter-a-zero",
        ),
        pytest.param(
            {"a": 5.0},
            "a",
            -3.0,
            pytest.raises(ValueError, match=re.escape("strictly positive")),
            id="setter-a-negative",
        ),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0},
            "alpha",
            0.0,
            pytest.raises(ValueError, match=re.escape("must be in")),
            id="setter-alpha-zero",
        ),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0},
            "beta",
            180.0,
            pytest.raises(ValueError, match=re.escape("must be in")),
            id="setter-beta-180",
        ),
        pytest.param(
            {
                "a": 3.0,
                "b": 4.0,
                "c": 5.0,
                "alpha": 160.0,
                "beta": 160.0,
                "gamma": 30.0,
            },
            "alpha",
            175.0,
            pytest.raises(ValueError, match=re.escape("degenerate")),
            id="setter-alpha-makes-cell-degenerate",
        ),
    ],
)
def test_lattice_setter_invalid(initial_kwargs, setter, bad_value, context):
    with context:
        lat = Lattice(**initial_kwargs)
        setattr(lat, setter, bad_value)


# ---------------------------------------------------------------------------
# Setter preserves state on failure
# ---------------------------------------------------------------------------


def test_lattice_setter_preserves_state_on_failure():
    lat = Lattice(a=5.0)
    original_a = lat.a
    with pytest.raises(ValueError):
        lat.a = -1.0
    # State must be unchanged after failed setter
    assert lat.a == original_a
    assert lat.system == "cubic"


# ---------------------------------------------------------------------------
# Lazy computation and cache invalidation
# ---------------------------------------------------------------------------


def test_lattice_lazy_not_computed_at_construction():
    lat = Lattice(a=5.0)
    assert lat._cartesian_lattice_vectors is None
    assert lat._reciprocal_lattice_vectors is None
    assert lat._B is None


def test_lattice_lazy_computed_after_access():
    """Each lazy property populates its own cache when accessed.

    Under issue #280, the three properties (``cartesian_lattice_vectors``,
    ``reciprocal_lattice_vectors``, ``B``) are computed independently
    from the lattice parameters via :func:`b_matrix_bl1967` and the
    BL1967 reciprocal-of-reciprocal identity, rather than chained
    through one another.  Accessing ``B`` therefore only populates
    ``_B``; the other caches populate on their own first access.
    """
    lat = Lattice(a=5.0)
    _ = lat.B
    assert lat._B is not None
    _ = lat.cartesian_lattice_vectors
    assert lat._cartesian_lattice_vectors is not None
    _ = lat.reciprocal_lattice_vectors
    assert lat._reciprocal_lattice_vectors is not None


@pytest.mark.parametrize(
    "initial_kwargs, param, new_value, context",
    [
        pytest.param(
            {"a": 5.0},
            "a",
            6.0,
            does_not_raise(),
            id="cache-invalidated-on-a",
        ),
        pytest.param(
            {"a": 5.0, "c": 10.0, "gamma": 120.0},
            "c",
            8.0,
            does_not_raise(),
            id="cache-invalidated-on-c",
        ),
        pytest.param(
            {"a": 5.0, "c": 10.0, "gamma": 120.0},
            "gamma",
            120.0,
            does_not_raise(),
            id="cache-invalidated-on-gamma",
        ),
    ],
)
def test_lattice_cache_invalidated_on_setter(initial_kwargs, param, new_value, context):
    with context:
        lat = Lattice(**initial_kwargs)
        _ = lat.B
        assert lat._B is not None
        setattr(lat, param, new_value)
        assert lat._cartesian_lattice_vectors is None
        assert lat._reciprocal_lattice_vectors is None
        assert lat._B is None


def test_lattice_cache_reused():
    lat = Lattice(a=5.0)
    v1 = lat.cartesian_lattice_vectors
    v2 = lat.cartesian_lattice_vectors
    assert v1 is v2


# ---------------------------------------------------------------------------
# Correctness of computed properties
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, expected_a1, context",
    [
        pytest.param(
            {"a": 5.0},
            np.array([5.0, 0.0, 0.0]),
            does_not_raise(),
            id="cubic-a1-along-x",
        ),
        pytest.param(
            {"a": 2.0, "b": 3.0, "c": 4.0},
            np.array([2.0, 0.0, 0.0]),
            does_not_raise(),
            id="orthorhombic-a1-along-x",
        ),
    ],
)
def test_lattice_cartesian_a1_along_x(kwargs, expected_a1, context):
    with context:
        lat = Lattice(**kwargs)
        a1, _, _ = lat.cartesian_lattice_vectors
        np.testing.assert_allclose(a1, expected_a1, atol=1e-12)


@pytest.mark.parametrize(
    "kwargs, context",
    [
        pytest.param({"a": 5.0}, does_not_raise(), id="cubic"),
        pytest.param(
            {"a": 2.0, "b": 3.0, "c": 4.0}, does_not_raise(), id="orthorhombic"
        ),
        pytest.param(
            {"a": 4.785, "c": 12.991, "gamma": 120.0}, does_not_raise(), id="hexagonal"
        ),
        pytest.param(
            {"a": 5.0, "b": 6.0, "c": 7.0, "beta": 110.0},
            does_not_raise(),
            id="monoclinic",
        ),
        pytest.param({"a": 5.0, "alpha": 60.0}, does_not_raise(), id="trigonal"),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 80.0, "beta": 85.0, "gamma": 95.0},
            does_not_raise(),
            id="triclinic",
        ),
    ],
)
def test_lattice_reciprocal_orthogonality(kwargs, context):
    """b_i . a_j = 2*pi * delta_ij for all i, j."""
    with context:
        lat = Lattice(**kwargs)
        a1, a2, a3 = lat.cartesian_lattice_vectors
        b1, b2, b3 = lat.reciprocal_lattice_vectors
        twopi = 2 * np.pi
        np.testing.assert_allclose(np.dot(b1, a1), twopi, atol=1e-10)
        np.testing.assert_allclose(np.dot(b2, a2), twopi, atol=1e-10)
        np.testing.assert_allclose(np.dot(b3, a3), twopi, atol=1e-10)
        np.testing.assert_allclose(np.dot(b1, a2), 0.0, atol=1e-10)
        np.testing.assert_allclose(np.dot(b1, a3), 0.0, atol=1e-10)
        np.testing.assert_allclose(np.dot(b2, a1), 0.0, atol=1e-10)
        np.testing.assert_allclose(np.dot(b2, a3), 0.0, atol=1e-10)
        np.testing.assert_allclose(np.dot(b3, a1), 0.0, atol=1e-10)
        np.testing.assert_allclose(np.dot(b3, a2), 0.0, atol=1e-10)


@pytest.mark.parametrize(
    "a, context",
    [
        pytest.param(5.0, does_not_raise(), id="cubic-a5"),
        pytest.param(3.0, does_not_raise(), id="cubic-a3"),
        pytest.param(1.0, does_not_raise(), id="cubic-a1-default"),
    ],
)
def test_lattice_B_cubic_diagonal(a, context):
    """For cubic, B is diagonal with 2π/a on the diagonal (BL1967 convention)."""
    with context:
        lat = Lattice(a=a)
        B = lat.B
        twopi_over_a = 2 * np.pi / a
        np.testing.assert_allclose(
            np.diag(B), [twopi_over_a, twopi_over_a, twopi_over_a], atol=1e-12
        )
        off = B - np.diag(np.diag(B))
        np.testing.assert_allclose(off, np.zeros((3, 3)), atol=1e-12)


@pytest.mark.parametrize(
    "kwargs, context",
    [
        pytest.param({"a": 5.0}, does_not_raise(), id="cubic"),
        pytest.param(
            {"a": 4.785, "c": 12.991, "gamma": 120.0}, does_not_raise(), id="hexagonal"
        ),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 80.0, "beta": 85.0, "gamma": 95.0},
            does_not_raise(),
            id="triclinic",
        ),
    ],
)
def test_lattice_B_bl1967_convention(kwargs, context):
    """Verify (b1, b2, b3) are the columns of B (BL1967 / SPEC convention).

    With 2π included in the reciprocal vectors, the columns of B are
    b1, b2, b3, so ``column_stack([b1, b2, b3]) == B`` and
    ``B @ h == h*b1 + k*b2 + l*b3``.
    """
    with context:
        lat = Lattice(**kwargs)
        b1, b2, b3 = lat.reciprocal_lattice_vectors
        B = lat.B
        rec_matrix = np.column_stack([b1, b2, b3])
        np.testing.assert_allclose(rec_matrix, B, atol=1e-10)
        # Verify the BL1967 action: B @ h == h*b1 + k*b2 + l*b3
        for hkl in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1]):
            expected = hkl[0] * b1 + hkl[1] * b2 + hkl[2] * b3
            np.testing.assert_allclose(B @ np.array(hkl), expected, atol=1e-10)


# ---------------------------------------------------------------------------
# Lower-symmetry B-matrix correctness (issue #237)
#
# These tests exercise non-orthogonal crystal systems (hexagonal, trigonal,
# monoclinic, triclinic) to confirm that ``Q = B @ h`` returns the correct
# scattering vector, with magnitude equal to ``2π/d_hkl`` computed
# independently from the reciprocal metric tensor.
#
# This is the regression coverage for issue #237, where ``b_matrix()`` had a
# spurious ``.T`` that placed the reciprocal vectors as rows of B instead of
# columns.  The bug was silent for cubic/tetragonal/orthorhombic cells (whose
# reciprocal vectors are mutually orthogonal Cartesian axes, so ``M`` and
# ``M.T`` are equal) but produced wrong results for any ``(h, k, 0)``-type
# reflection in a cell with non-orthogonal reciprocal vectors.
# ---------------------------------------------------------------------------


def _d_spacing_from_metric(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
    h: int,
    k: int,
    l: int,  # noqa: E741
) -> float:
    """
    Compute ``d_hkl`` from the reciprocal metric tensor, independent of B.

    Uses the standard crystallographic formula ``1/d² = h G* h^T``, where
    ``G* = (G)⁻¹`` and ``G`` is the real-space metric tensor:

        G_{ij} = a_i · a_j

    All angles in degrees.  Returns d in Angstroms.
    """
    ar, br_, gr = np.deg2rad([alpha, beta, gamma])
    G = np.array(
        [
            [a * a, a * b * np.cos(gr), a * c * np.cos(br_)],
            [a * b * np.cos(gr), b * b, b * c * np.cos(ar)],
            [a * c * np.cos(br_), b * c * np.cos(ar), c * c],
        ]
    )
    Gstar = np.linalg.inv(G)
    h_vec = np.array([h, k, l], dtype=float)
    inv_d_sq = float(h_vec @ Gstar @ h_vec)
    return 1.0 / np.sqrt(inv_d_sq)


_LOW_SYMMETRY_LATTICES = [
    pytest.param(
        # hexagonal sapphire (α-Al₂O₃)
        dict(a=4.785, c=12.991, gamma=120.0),
        dict(a=4.785, b=4.785, c=12.991, alpha=90.0, beta=90.0, gamma=120.0),
        id="hexagonal-sapphire",
    ),
    pytest.param(
        # trigonal (rhombohedral setting), e.g. dolomite-like
        dict(a=4.81, alpha=47.0),
        dict(a=4.81, b=4.81, c=4.81, alpha=47.0, beta=47.0, gamma=47.0),
        id="trigonal-rhombohedral",
    ),
    pytest.param(
        # monoclinic gypsum (CaSO₄·2H₂O)
        dict(a=6.284, b=15.200, c=5.678, beta=114.09),
        dict(a=6.284, b=15.200, c=5.678, alpha=90.0, beta=114.09, gamma=90.0),
        id="monoclinic-gypsum",
    ),
    pytest.param(
        # triclinic K₂Cr₂O₇-like cell (all six parameters distinct from 90°)
        dict(a=7.45, b=7.38, c=13.40, alpha=98.0, beta=96.5, gamma=91.5),
        dict(a=7.45, b=7.38, c=13.40, alpha=98.0, beta=96.5, gamma=91.5),
        id="triclinic-K2Cr2O7-like",
    ),
]


@pytest.mark.parametrize(
    "lattice_kwargs, full_params",
    [(p.values[0], p.values[1]) for p in _LOW_SYMMETRY_LATTICES],
    ids=[p.id for p in _LOW_SYMMETRY_LATTICES],
)
@pytest.mark.parametrize(
    "h, k, l",
    [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 1, 0),  # mixes b1 and b2 — exposes B-matrix layout bug
        (1, 0, 1),  # mixes b1 and b3
        (0, 1, 1),  # mixes b2 and b3
        (1, 1, 1),  # mixes all three
        (2, -1, 0),  # negative index
        (3, 2, 1),  # generic
    ],
)
def test_lattice_B_lower_symmetry_d_spacing(lattice_kwargs, full_params, h, k, l):  # noqa: E741
    """
    For non-orthogonal cells, ``|B @ h| == 2π / d_hkl``.

    ``d_hkl`` is computed independently from the reciprocal metric tensor
    (no use of B), so this test verifies the B-matrix layout end-to-end.

    Regression test for issue #237.
    """
    lat = Lattice(**lattice_kwargs)
    B = lat.B
    hkl = np.array([h, k, l], dtype=float)
    Q = B @ hkl
    Q_mag = float(np.linalg.norm(Q))

    d_expected = _d_spacing_from_metric(**full_params, h=h, k=k, l=l)
    Q_expected_mag = 2.0 * np.pi / d_expected

    np.testing.assert_allclose(Q_mag, Q_expected_mag, rtol=1e-10)


@pytest.mark.parametrize(
    "lattice_kwargs",
    [
        pytest.param(dict(a=4.785, c=12.991, gamma=120.0), id="hexagonal-sapphire"),
        pytest.param(dict(a=4.81, alpha=47.0), id="trigonal-rhombohedral"),
        pytest.param(
            dict(a=6.284, b=15.200, c=5.678, beta=114.09), id="monoclinic-gypsum"
        ),
        pytest.param(
            dict(a=7.45, b=7.38, c=13.40, alpha=98.0, beta=96.5, gamma=91.5),
            id="triclinic-K2Cr2O7-like",
        ),
    ],
)
def test_lattice_B_lower_symmetry_columns_are_reciprocal(lattice_kwargs):
    """
    For non-orthogonal cells, the columns of B are the reciprocal lattice vectors.

    This is the layout assertion that caught issue #237: ``B @ (1, 0, 0) == b1``
    only holds when columns of B are b1, b2, b3 (BL1967 / SPEC convention).
    The pre-fix implementation used ``B = column_stack([b1, b2, b3]).T`` which
    placed the reciprocal vectors as rows, causing ``B @ (1, 0, 0)`` to return
    ``(b1·x̂, b2·x̂, b3·x̂)`` instead of ``b1``.

    Regression test for issue #237.
    """
    lat = Lattice(**lattice_kwargs)
    b1, b2, b3 = lat.reciprocal_lattice_vectors
    B = lat.B

    np.testing.assert_allclose(B[:, 0], b1, atol=1e-12)
    np.testing.assert_allclose(B[:, 1], b2, atol=1e-12)
    np.testing.assert_allclose(B[:, 2], b3, atol=1e-12)

    # B @ unit-h returns the corresponding reciprocal vector.
    np.testing.assert_allclose(B @ np.array([1.0, 0.0, 0.0]), b1, atol=1e-12)
    np.testing.assert_allclose(B @ np.array([0.0, 1.0, 0.0]), b2, atol=1e-12)
    np.testing.assert_allclose(B @ np.array([0.0, 0.0, 1.0]), b3, atol=1e-12)


# ---------------------------------------------------------------------------
# Setter re-deduction of crystal system
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "initial_kwargs, setter, new_value, expected_system, context",
    [
        pytest.param(
            {"a": 5.0},
            "c",
            8.0,
            "tetragonal",
            does_not_raise(),
            id="cubic-becomes-tetragonal-on-c-set",
        ),
        pytest.param(
            {"a": 5.0},
            "b",
            7.0,
            "orthorhombic",
            does_not_raise(),
            id="cubic-becomes-orthorhombic-on-b-set",
        ),
        pytest.param(
            {"a": 5.0},
            "a",
            3.0,
            "cubic",
            does_not_raise(),
            id="cubic-stays-cubic-on-a-change",
        ),
    ],
)
def test_lattice_setter_rededuces_system(
    initial_kwargs, setter, new_value, expected_system, context
):
    with context:
        lat = Lattice(**initial_kwargs)
        setattr(lat, setter, new_value)
        assert lat.system == expected_system


# ---------------------------------------------------------------------------
# __str__ reports only free parameters for the crystal system
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, expected_system, expected_in_str, expected_not_in_str, context",
    [
        # Expected fragments use default precision=6 (package default)
        pytest.param(
            {"a": 5.0},
            "cubic",
            ["cubic", "a=5.000000"],
            ["b=", "c=", "alpha=", "beta=", "gamma="],
            does_not_raise(),
            id="str-cubic",
        ),
        pytest.param(
            {"a": 3.0, "c": 6.0},
            "tetragonal",
            ["tetragonal", "a=3.000000", "c=6.000000"],
            ["b=", "alpha=", "beta=", "gamma="],
            does_not_raise(),
            id="str-tetragonal",
        ),
        pytest.param(
            {"a": 4.785, "c": 12.991, "gamma": 120.0},
            "hexagonal",
            ["hexagonal", "a=4.785000", "c=12.991000"],
            ["b=", "alpha=", "beta=", "gamma="],
            does_not_raise(),
            id="str-hexagonal",
        ),
        pytest.param(
            {"a": 5.0, "alpha": 60.0},
            "trigonal",
            ["trigonal", "a=5.000000", "alpha=60.000000"],
            ["b=", "c=", "beta=", "gamma="],
            does_not_raise(),
            id="str-trigonal",
        ),
        pytest.param(
            {"a": 2.0, "b": 3.0, "c": 4.0},
            "orthorhombic",
            ["orthorhombic", "a=2.000000", "b=3.000000", "c=4.000000"],
            ["alpha=", "beta=", "gamma="],
            does_not_raise(),
            id="str-orthorhombic",
        ),
        pytest.param(
            {"a": 5.0, "b": 6.0, "c": 7.0, "beta": 110.0},
            "monoclinic",
            ["monoclinic", "a=5.000000", "b=6.000000", "c=7.000000", "beta=110.000000"],
            ["alpha=", "gamma="],
            does_not_raise(),
            id="str-monoclinic",
        ),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 80.0, "beta": 85.0, "gamma": 95.0},
            "triclinic",
            [
                "triclinic",
                "a=3.000000",
                "b=4.000000",
                "c=5.000000",
                "alpha=80.000000",
                "beta=85.000000",
                "gamma=95.000000",
            ],
            [],
            does_not_raise(),
            id="str-triclinic",
        ),
    ],
)
def test_lattice_str(
    kwargs, expected_system, expected_in_str, expected_not_in_str, context
):
    with context:
        lat = Lattice(**kwargs)
        s = str(lat)
        assert expected_system in s
        for fragment in expected_in_str:
            assert fragment in s, f"Expected {fragment!r} in str: {s!r}"
        for fragment in expected_not_in_str:
            assert fragment not in s, f"Unexpected {fragment!r} in str: {s!r}"


# ---------------------------------------------------------------------------
# __repr__ reports all six parameters
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, context",
    [
        pytest.param({"a": 5.0}, does_not_raise(), id="repr-cubic"),
        pytest.param(
            {"a": 2.0, "b": 3.0, "c": 4.0}, does_not_raise(), id="repr-orthorhombic"
        ),
        pytest.param(
            {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 80.0, "beta": 85.0, "gamma": 95.0},
            does_not_raise(),
            id="repr-triclinic",
        ),
    ],
)
def test_lattice_repr(kwargs, context):
    with context:
        lat = Lattice(**kwargs)
        r = repr(lat)
        assert "Lattice" in r
        assert lat.system in r
        for param in ("a", "b", "c", "alpha", "beta", "gamma"):
            assert f"{param}=" in r


# ---------------------------------------------------------------------------
# Display precision
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "precision, expected_fragment, context",
    [
        pytest.param(
            3,
            "a=5.000",
            does_not_raise(),
            id="precision-3-decimal-places",
        ),
        pytest.param(
            0,
            "a=5",
            does_not_raise(),
            id="precision-0-decimal-places",
        ),
        pytest.param(
            10,
            "a=5.0000000000",
            does_not_raise(),
            id="precision-10-decimal-places",
        ),
        pytest.param(
            None,
            "a=5.000000",
            does_not_raise(),
            id="precision-none-uses-package-default",
        ),
    ],
)
def test_lattice_instance_precision(precision, expected_fragment, context):
    """Per-instance precision controls __str__ output."""
    with context:
        lat = Lattice(a=5.0, precision=precision)
        assert expected_fragment in str(lat)


@pytest.mark.parametrize(
    "set_digits, expected_fragment, context",
    [
        pytest.param(2, "a=5.00", does_not_raise(), id="package-precision-2"),
        pytest.param(4, "a=5.0000", does_not_raise(), id="package-precision-4"),
    ],
)
def test_package_precision(set_digits, expected_fragment, context):
    """Package-level set_precision() controls __str__ for instances with precision=None."""
    with context:
        original = ahd.display.get_precision()
        try:
            ahd.display.set_precision(set_digits)
            lat = Lattice(a=5.0)  # precision=None -> uses package default
            assert expected_fragment in str(lat)
        finally:
            ahd.display.set_precision(original)  # always restore


@pytest.mark.parametrize(
    "digits, context",
    [
        pytest.param(0, does_not_raise(), id="set-precision-zero"),
        pytest.param(6, does_not_raise(), id="set-precision-six"),
        pytest.param(15, does_not_raise(), id="set-precision-fifteen"),
        pytest.param(
            -1,
            pytest.raises(ValueError, match=re.escape("non-negative")),
            id="invalid-set-precision-negative",
        ),
        pytest.param(
            "six",
            pytest.raises(ValueError, match=re.escape("non-negative integer")),
            id="invalid-set-precision-string",
        ),
        pytest.param(
            3.0,
            pytest.raises(ValueError, match=re.escape("non-negative integer")),
            id="invalid-set-precision-float",
        ),
    ],
)
def test_set_precision_validation(digits, context):
    with context:
        original = ahd.display.get_precision()
        try:
            ahd.display.set_precision(digits)
        finally:
            ahd.display.set_precision(original)


def test_instance_precision_overrides_package():
    """Instance precision takes precedence over the package-level default."""
    original = ahd.display.get_precision()
    try:
        ahd.display.set_precision(6)
        lat = Lattice(a=5.0, precision=2)
        assert "a=5.00 " in str(lat)
        assert "a=5.000000" not in str(lat)
    finally:
        ahd.display.set_precision(original)


# ---------------------------------------------------------------------------
# Lattice.__eq__ with tolerance
# ---------------------------------------------------------------------------


def test_lattice_eq_identical():
    assert Lattice(a=5.431) == Lattice(a=5.431)


def test_lattice_eq_within_default_tolerance():
    """Values differing by less than half a unit in the 6th decimal place."""
    a = 5.4310000
    b = 5.4310003  # differs by 3e-7, within atol=5e-7
    assert Lattice(a=a) == Lattice(a=b)


def test_lattice_eq_outside_default_tolerance():
    a = 5.4310000
    b = 5.4320000  # differs by 1e-3, outside atol=5e-7
    assert Lattice(a=a) != Lattice(a=b)


def test_lattice_eq_explicit_atol():
    assert Lattice(a=5.431).__eq__(Lattice(a=5.432), atol=0.01) is True
    assert Lattice(a=5.431).__eq__(Lattice(a=5.451), atol=0.01) is False


def test_lattice_eq_not_implemented_for_non_lattice():
    assert Lattice(a=5.431).__eq__("not a lattice") is NotImplemented


def test_lattice_eq_different_parameters():
    assert Lattice(a=5.431) != Lattice(a=5.431, c=10.0)  # tetragonal vs cubic


# ---------------------------------------------------------------------------
# Lattice.to_dict() / from_dict()
# ---------------------------------------------------------------------------


import json  # noqa: E402  (needed for json.dumps check)


def test_lattice_to_dict_structure():
    """to_dict() returns a JSON-serialisable dict with exactly the six keys."""
    d = Lattice(a=5.431).to_dict()
    assert isinstance(d, dict)
    assert set(d.keys()) == {"a", "b", "c", "alpha", "beta", "gamma"}
    assert json.dumps(d)  # must not raise


@pytest.mark.parametrize(
    "lattice, key, expected, context",
    [
        pytest.param(Lattice(a=5.431), "a", 5.431, does_not_raise(), id="cubic-a"),
        pytest.param(
            Lattice(a=5.431), "b", 5.431, does_not_raise(), id="cubic-b-equals-a"
        ),
        pytest.param(
            Lattice(a=5.431), "alpha", 90.0, does_not_raise(), id="cubic-alpha"
        ),
        pytest.param(
            Lattice(a=4.785, c=12.991, gamma=120),
            "gamma",
            120.0,
            does_not_raise(),
            id="hex-gamma",
        ),
        pytest.param(
            Lattice(a=4.785, c=12.991, gamma=120),
            "a",
            4.785,
            does_not_raise(),
            id="hex-a",
        ),
        pytest.param(
            Lattice(a=4.785, c=12.991, gamma=120),
            "c",
            12.991,
            does_not_raise(),
            id="hex-c",
        ),
    ],
)
def test_lattice_to_dict_values(lattice, key, expected, context):
    """to_dict() stores the correct value for each lattice parameter."""
    with context:
        assert lattice.to_dict()[key] == pytest.approx(expected)


@pytest.mark.parametrize(
    "lattice, context",
    [
        pytest.param(Lattice(a=5.431), does_not_raise(), id="cubic"),
        pytest.param(
            Lattice(a=4.785, c=12.991, gamma=120), does_not_raise(), id="hexagonal"
        ),
        pytest.param(
            Lattice(a=5.0, b=6.0, c=7.0, alpha=80.0, beta=90.0, gamma=100.0),
            does_not_raise(),
            id="triclinic",
        ),
    ],
)
def test_lattice_from_dict_roundtrip(lattice, context):
    """from_dict(to_dict()) reproduces the original Lattice exactly."""
    with context:
        assert Lattice.from_dict(lattice.to_dict()) == lattice
