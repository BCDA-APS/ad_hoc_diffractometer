# *Ad hoc* diffractometer

Multi-circle diffractometer geometry and related calculations.

## References

Papers describing various diffractometer geometries.

- 1967 (fourc), Acta  Cryst., 22, 457-464, W.R. Busing and H.A. Levy.
- 1999 (psic), J. Appl. Cryst. 32, 614-623, H. You.
- 1993 (sixc), J. Appl. Cryst., 26, 706, M. Lohmeier and E. Vlieg.
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
