"""
display.py — display precision settings.

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
"""

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
