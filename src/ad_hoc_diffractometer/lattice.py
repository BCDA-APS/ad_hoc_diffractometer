"""
lattice.py — crystallographic lattice calculations.

Provides functions for computing direct lattice vectors, reciprocal lattice
vectors, and the B matrix from crystal lattice parameters.

Based on:
  Busing & Levy, Acta Cryst. 22, 457-464 (1967)
  I16 Diffractometer Rotation Matrix document
  Lecture-2-Reciprocal-lattice-notes
"""

import numpy as np


def lattice_vectors(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the three Cartesian direct lattice vectors from crystal lattice
    parameters.

    The convention places a1 along xHat, a2 in the xHat-yHat plane, and a3
    determined by the right-hand rule.

    Parameters
    ----------
    a, b, c : float
        Lattice parameters in Angstroms.
    alpha, beta, gamma : float
        Lattice angles in degrees.
        alpha = angle between b and c axes
        beta  = angle between a and c axes
        gamma = angle between a and b axes

    Returns
    -------
    a1, a2, a3 : numpy.ndarray, shape (3,)
        Cartesian direct lattice vectors in Angstroms.
    """
    alpha_r = np.deg2rad(alpha)
    beta_r = np.deg2rad(beta)
    gamma_r = np.deg2rad(gamma)

    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([b * np.cos(gamma_r), b * np.sin(gamma_r), 0.0])

    a3x = c * np.cos(beta_r)
    a3y = c * (np.cos(alpha_r) - np.cos(beta_r) * np.cos(gamma_r)) / np.sin(gamma_r)
    a3z = np.sqrt(max(c**2 - a3x**2 - a3y**2, 0.0))
    a3 = np.array([a3x, a3y, a3z])

    return a1, a2, a3


def reciprocal_vectors(
    a1: np.ndarray,
    a2: np.ndarray,
    a3: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the three reciprocal lattice vectors from Cartesian direct lattice
    vectors.

    Uses the standard crystallographic definition:

        b1 = 2*pi * (a2 x a3) / (a1 . (a2 x a3))
        b2 = 2*pi * (a3 x a1) / (a1 . (a2 x a3))
        b3 = 2*pi * (a1 x a2) / (a1 . (a2 x a3))

    where the denominator is the unit cell volume Vc = a1 . (a2 x a3).

    The orthogonality condition is satisfied:

        b_i . a_j = 2*pi * delta_ij

    Parameters
    ----------
    a1, a2, a3 : numpy.ndarray, shape (3,)
        Cartesian direct lattice vectors in Angstroms.

    Returns
    -------
    b1, b2, b3 : numpy.ndarray, shape (3,)
        Reciprocal lattice vectors in inverse Angstroms (with 2*pi factor).
    """
    Vc = np.dot(a1, np.cross(a2, a3))

    b1 = 2 * np.pi * np.cross(a2, a3) / Vc
    b2 = 2 * np.pi * np.cross(a3, a1) / Vc
    b3 = 2 * np.pi * np.cross(a1, a2) / Vc

    return b1, b2, b3


def b_matrix(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
) -> np.ndarray:
    """
    Compute the B matrix from the reciprocal lattice vectors.

    The B matrix transforms Miller indices h = (h, k, l) to Cartesian
    coordinates in the crystal frame (Busing & Levy, 1967, eq. 3):

        hc = B . h

    Following the I16 diffractometer convention:

        (b1, b2, b3) = 2*pi * B.T

    so the columns of B.T are the reciprocal lattice vectors divided by 2*pi,
    equivalently:

        2*pi * B.T @ h = h[0]*b1 + h[1]*b2 + h[2]*b3

    B is not in general orthonormal (the crystal need not be cubic).

    Parameters
    ----------
    b1, b2, b3 : numpy.ndarray, shape (3,)
        Reciprocal lattice vectors in inverse Angstroms (with 2*pi factor).

    Returns
    -------
    B : numpy.ndarray, shape (3, 3)
        B matrix in inverse Angstroms (no 2*pi factor).
    """
    return np.column_stack([b1, b2, b3]).T / (2 * np.pi)
