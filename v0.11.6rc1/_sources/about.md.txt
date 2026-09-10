(about)=
# About `ad_hoc_diffractometer`

`ad_hoc_diffractometer` handles the core calculations you need for
diffractometer work.

## Features

- **Diffractometer description** — describe your diffractometer using
  observable physical directions (vertical, longitudinal, transverse).
- **Orientation calculations** — compute orientation matrices from
  reflections; refine the crystal lattice.
- **Reciprocal-space mapping** — convert from rotation axes to
  reciprocal-space coordinates.
- **Diffractometer control** — convert from reciprocal-space coordinates
  to rotation axes.
- **Mode definitions** — you define which axes are free, fixed, or
  coupled when solving the kinematics.

## It's your diffractometer

You get full control over your setup — whether you're using a standard
four-circle diffractometer or something custom, the package adapts to
you.  No hard-coded diffractometers: a new instrument requires zero
changes to the package code.

## Minimal requirements

Only [Python](https://python.org) (with its Standard Library),
[NumPy](https://numpy.org), and [YAML](https://pyyaml.org/).  No SciPy, SymPy,
or other scientific dependencies required.

## Use cases

- Simulating diffractometer behavior.
- Real-time operations in reciprocal space during beamtime.
- Planning experiments and trajectories before you run them.
- Backend support for diffractometer control systems.
- Creating visualizations of diffractometer geometry.

```{note}
The package assumes **monochromatic radiation** throughout — all
calculations are at a fixed wavelength.
```

## See also

- {doc}`quick_start` — build a four-circle diffractometer step by step
- {ref}`Demonstration Geometries Quick start <quick-start-demo>` — a
  short example using a demo geometry
- {doc}`concepts` — coordinate conventions, axis sign convention,
  B/U/UB matrices, diffraction modes, and the ψ angle
- {doc}`howto/index` — task-oriented how-to guides
