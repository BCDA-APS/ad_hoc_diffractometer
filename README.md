# *Ad hoc* diffractometer

[![Documentation](https://img.shields.io/badge/docs-latest-blue)](https://prjemian.github.io/ad_hoc_diffractometer/latest/)
[![License: CC-BY-4.0](https://img.shields.io/badge/License-CC--BY--4.0-brightgreen)](https://creativecommons.org/licenses/by/4.0/)

Multi-circle diffractometer geometry and related calculations.

`import ad_hoc_diffractometer as ahd`

## References

Papers describing various diffractometer geometries.

- 1967 (fourc), Acta  Cryst., 22, 457-464, W.R. Busing and H.A. Levy.
- 1985 (zaxis), J. Appl. Cryst. 18, 33-36, J.M. Bloch.
- 1987 (fivec), J. Appl. Cryst. 20, 330-337, E. Vlieg et al.
- 1993 (sixc), J. Appl. Cryst., 26, 706, M. Lohmeier and E. Vlieg.
- 1995 (s2d2), J. Appl. Cryst. 28, 318-326, K.W. Evans-Lutterodt & M.-T. Tang.
- 1999 (psic), J. Appl. Cryst. 32, 614-623, H. You.
- 2006 (kappa), International Tables for Crystallography, Vol. C, Section 2.2.6. DOI:
    10.1107/97809553602060000577.
- 2016 (review), "Multicircle Diffractometry Methods," pp. 1-10, D.A. Walko,
  in "Reference Module in Materials Science and Materials Engineering,"
  Oxford: Elsevier,  Saleem Hashmi (editor-in-chief), ISBN: 978-0-12-803581-8.

Also note:

- https://en.wikipedia.org/wiki/Orientation_(geometry)
- https://www.cuemath.com/algebra/rotation-matrix/
- https://en.wikipedia.org/wiki/Rotation_matrix
- https://en.wikipedia.org/wiki/Quaternion
- https://quaternion.readthedocs.io/en/latest/ (Python)
- https://sot.github.io/Quaternion/ (Python)

## INSTALL

### conda/pip

```bash
conda create -y -n ad_hoc_diffractometer python
conda activate ad_hoc_diffractometer
pip install -e .[dev]
```

### uv

[uv](https://docs.astral.sh/uv/) creates and manages virtual environments automatically:

```bash
uv sync --extra dev
```

Run commands inside the managed environment with `uv run`, e.g.:

```bash
uv run pytest
```

### hatch

[hatch](https://hatch.pypa.io/) uses the `default` environment defined in
`pyproject.toml`, which already includes the `dev` dependencies:

```bash
pip install hatch          # install hatch once, globally
hatch env create           # create the default environment
hatch run pytest           # run commands inside the environment
```
