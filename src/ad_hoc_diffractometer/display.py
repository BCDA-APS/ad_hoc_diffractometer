"""
display.py — display precision settings and numeric comparison helper.

Provides a package-level default number of decimal places for displaying
floating-point values.  Any class or function that formats numbers for
display should use get_precision() to honour the current setting.

The default precision can be changed at the package level:

    import ad_hoc_diffractometer as ahd
    ahd.set_precision(4)

Or per-object where the class supports it (e.g. Lattice):

    lat = Lattice(a=4.785, c=12.991, gamma=120.0)
    lat.precision = 4
    print(lat)   # shows 4 decimal places

Internal representation is always full floating-point precision.  Only
display (str, repr-like summaries) is affected.

``isclose(a, b)`` provides a tolerance-aware numeric comparison whose
absolute tolerance is derived from the current display precision:
``atol = 0.5 * 10 ** (-get_precision())``.  This means two values are
considered equal when they agree to within half a unit in the last
*displayed* decimal place — the same resolution the user sees.  An
explicit tolerance can be supplied to override this default.
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Package-level default
# ---------------------------------------------------------------------------

_DEFAULT_PRECISION: int = 6  # decimal places


def get_precision() -> int:
    """Return the current package-level display precision (decimal places)."""
    return _DEFAULT_PRECISION


def set_precision(digits: int) -> None:
    """
    Set the package-level display precision.

    Parameters
    ----------
    digits : int
        Number of decimal places for displaying floating-point values.
        Must be a non-negative integer.

    Raises
    ------
    ValueError
        If digits is negative or not an integer.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> ahd.set_precision(4)
    >>> ahd.get_precision()
    4
    """
    global _DEFAULT_PRECISION
    if not isinstance(digits, int) or isinstance(digits, bool):
        raise ValueError(f"Precision must be a non-negative integer; got {digits!r}.")
    if digits < 0:
        raise ValueError(f"Precision must be non-negative; got {digits}.")
    _DEFAULT_PRECISION = digits


def precision_atol(digits: int | None = None) -> float:
    """
    Return an absolute tolerance matching the current display precision.

    ``atol = 0.5 * 10 ** (-digits)`` — half a unit in the last displayed
    decimal place, i.e. the smallest difference that would show up in the
    formatted output.

    Parameters
    ----------
    digits : int or None
        Number of decimal places.  If None, uses ``get_precision()``.

    Returns
    -------
    float
    """
    if digits is None:
        digits = get_precision()
    return 0.5 * 10 ** (-digits)


def allclose(
    a,
    b,
    atol: float | None = None,
    digits: int | None = None,
) -> bool:
    """
    Tolerance-aware comparison for scalars or array-like sequences.

    Uses ``np.allclose`` with ``rtol=0`` and an absolute tolerance derived
    from the display precision (or an explicit value).

    Parameters
    ----------
    a, b : scalar or array-like
        Values to compare.
    atol : float or None
        Explicit absolute tolerance.  If None, derived from ``digits``.
    digits : int or None
        Display precision to derive tolerance from.  If None, uses
        ``get_precision()``.  Ignored when ``atol`` is given explicitly.

    Returns
    -------
    bool

    Examples
    --------
    >>> allclose(1.0000001, 1.0000002)          # within default 6-digit tol
    True
    >>> allclose(1.001, 1.002, digits=3)        # within 3-digit tol (0.0005)
    True
    >>> allclose(1.001, 1.002, digits=4)        # outside 4-digit tol (0.00005)
    False
    >>> allclose([1.0, 2.0], [1.0, 2.0])
    True
    """
    if atol is None:
        atol = precision_atol(digits)
    return bool(np.allclose(a, b, atol=atol, rtol=0))


def fmt(value: float, digits: int | None = None) -> str:
    """
    Format a float for display with the given number of decimal places.

    Parameters
    ----------
    value : float
        The value to format.
    digits : int or None
        Number of decimal places.  If None, uses the package-level default
        from get_precision().


    Returns
    -------
    str
        Formatted string, e.g. '4.785000' for fmt(4.785, 6).
    """
    if digits is None:
        digits = get_precision()
    return f"{value:.{digits}f}"
