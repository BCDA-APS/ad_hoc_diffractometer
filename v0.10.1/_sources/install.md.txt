# Installation

## Requirements

- Python ≥ 3.10
- [NumPy](https://numpy.org/)

## Install from PyPI

```bash
pip install ad_hoc_diffractometer
```

## Install for development

```bash
git clone https://github.com/BCDA-APS/ad_hoc_diffractometer
cd ad_hoc_diffractometer
pip install -e .[dev]
```

## Install documentation dependencies

```bash
pip install -e .[doc]
```

## Install all optional dependencies

```bash
pip install -e .[all]
```

## Build the documentation locally

```bash
cd docs
make html
```

The built HTML is written to `docs/build/html/`.
Open `docs/build/html/index.html` in a browser to view it.

## Verify the installation

After installing, run the short example on the
{ref}`Demonstration Geometries Quick start <quick-start-demo>` page to
confirm the package is working end-to-end.
