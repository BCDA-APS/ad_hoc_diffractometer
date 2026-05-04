# *Ad hoc* diffractometer

[![PyPI](https://img.shields.io/pypi/v/ad_hoc_diffractometer)](https://pypi.org/project/ad_hoc_diffractometer/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue)](https://prjemian.github.io/ad_hoc_diffractometer/latest/)
[![License: CC-BY-4.0](https://img.shields.io/badge/License-CC--BY--4.0-brightgreen)](https://creativecommons.org/licenses/by/4.0/)

`ad_hoc_diffractometer` is a Python package that lets you describe any
multi-circle diffractometer geometry and perform X-ray/neutron
crystallography calculations.

See the [**Quick Start
guide**](https://prjemian.github.io/ad_hoc_diffractometer/latest/quick_start.html)
for a step-by-step walkthrough building an Eulerian four-circle geometry
— choosing a coordinate basis, stage stacking, diffraction mode
definition(s), and running a forward calculation. Common geometries are
provided as examples.

## Features

`ad_hoc_diffractometer` handles the core calculations you need for diffractometer work:

- **Geometry setup**: Describe your diffractometer using observable
  physical directions (vertical, longitudinal, lattransverseeral).
- **Orientation calculations**: Compute orientation matrices from
  reflections, refine crystal lattice.
- **Reciprocal space mapping**: Convert from rotation axes to reciprocal
  space coordinates.
- **Diffractometer control**: Convert from reciprocal space coordinates
  to rotation axes.
- **Mode Definitions**: You define which axes are free, fixed, or
  coupled when solving kinematics.

## It's your diffractometer

You get **full control over your setup** — whether you're using a
standard four-circle geometry or something custom, the package adapts to
you. No hard-coded configurations mean new geometries require zero
changes to the code.

## Minimal Requirements

Only [**Python**](https://python.org) with its Standard Library and
[**NumPy**](https://numpy.org). No scipy, sympy, or other scientific
dependencies required.

## Use Cases

- Simulating diffractometer behavior.
- Real-time operations in reciprocal space during beamtime.
- Planning experiments and trajectories before you run them.
- Backend support for diffractometer control systems.
- Creating visualizations of diffractometer geometry.

**Important**: The package assumes **monochromatic radiation**
throughout — all calculations are at a fixed wavelength.
