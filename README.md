# *Ad hoc* diffractometer

[![PyPI](https://img.shields.io/pypi/v/ad_hoc_diffractometer)](https://pypi.org/project/ad_hoc_diffractometer/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue)](https://prjemian.github.io/ad_hoc_diffractometer/latest/)
[![License: CC-BY-4.0](https://img.shields.io/badge/License-CC--BY--4.0-brightgreen)](https://creativecommons.org/licenses/by/4.0/)

`ad_hoc_diffractometer` is a **pure-Python** package for operating
multi-circle diffractometers in reciprocal space for X-ray and neutron
crystallography.  It is built around a key design principle: **any
multi-circle diffractometer geometry can be fully described by the caller**
— no geometry is hard-coded, and new geometries require no changes to the
package itself.

Its only runtime dependency beyond the Python Standard Library is
[NumPy](https://numpy.org) — no scipy, sympy, or other scientific
libraries are required.

> **Note:** The package assumes **monochromatic radiation** throughout.
> All diffraction calculations are performed at a fixed wavelength.

## Features

- A class-based description of diffractometer stages (rotary axes) and
  their stacking order
- Predefined factory functions for standard synchrotron and laboratory
  diffractometer geometries (psic, fourcv, fourch, sixc, kappa families,
  zaxis, s2d2, fivec)
- Crystallographic lattice calculations (B matrix, reciprocal lattice)
- U and UB matrix computation from orienting reflections
- Forward diffraction calculations (hkl → motor angles), with
  diffraction modes controlling which stages are free, fixed, or coupled
- Reciprocal-space trajectory planning

## Quick start

```python
import ad_hoc_diffractometer as ahd

# Create a six-circle psic geometry and set the wavelength
g = ahd.psic()
g.wavelength = 1.0  # Å

# Define the sample lattice (cubic silicon)
g.sample.lattice = ahd.Lattice(a=5.431)

# Show a summary of the geometry
print(g.summary())
```

## Install

```bash
pip install ad_hoc_diffractometer
```

For development:

```bash
git clone https://github.com/prjemian/ad_hoc_diffractometer
cd ad_hoc_diffractometer
pip install -e .[dev]
```

See the [documentation](https://prjemian.github.io/ad_hoc_diffractometer/latest/)
for full installation options (conda, uv, hatch) and usage guides.

## References

Papers describing the diffractometer geometries provided:

- Busing & Levy (1967) — fourc. *Acta Cryst.* 22, 457–464.
- Bloch (1985) — zaxis. *J. Appl. Cryst.* 18, 33–36.
- Vlieg et al. (1987) — fivec. *J. Appl. Cryst.* 20, 330–337.
- Lohmeier & Vlieg (1993) — sixc. *J. Appl. Cryst.* 26, 706–716.
- Evans-Lutterodt & Tang (1995) — s2d2. *J. Appl. Cryst.* 28, 318–326.
- You (1999) — psic. *J. Appl. Cryst.* 32, 614–623.
  DOI: [10.1107/S0021889899001223](https://doi.org/10.1107/S0021889899001223)
- ITC Vol. C §2.2.6 (2006) — kappa.
  DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)
- Walko (2016) — geometry survey. *Reference Module in Materials Science
  and Materials Engineering*, Elsevier.
