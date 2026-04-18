(howto-serialize)=
# Save and Restore a Diffractometer Configuration

This guide shows how to save the complete state of a diffractometer to a
file and restore it later.  The state includes the geometry, wavelength,
sample lattice, orienting reflections, UB matrix, and all other parameters.

## What is serialized

`to_dict()` captures the full diffractometer state:

- Geometry name, description, and basis
- Wavelength, source type
- All motor stage names and angles
- Active diffraction mode and cut points
- All samples: lattice parameters, reflections (with angles and wavelength),
  orienting reflection designations (or0/or1), and UB matrix
- Azimuthal reference, surface normal, inclination matrix
- Detector geometry parameters (distance, tilt, offset)

The resulting dict contains only JSON-compatible types (str, float, int,
list, dict, None) — it can be round-tripped through JSON or YAML without
loss.

## Basic round-trip

```python
import ad_hoc_diffractometer as ahd

# Set up a geometry
g = ahd.presets.fourcv()
g.wavelength = 1.5406                       # Å  (Cu Kα)
g.sample.lattice = ahd.Lattice(a=5.431)     # cubic silicon
ahd.ub_identity(g.sample)

g.sample.reflections.add(
    "r1",
    hkl=(0, 0, 4),
    angles={"omega": 34.56, "chi": 90.0, "phi": 0.0, "ttheta": 69.13},
)
g.sample.reflections.setor0("r1")

# Serialize to a dict
d = g.to_dict()

# Restore from the dict
g2 = ahd.AdHocDiffractometer.from_dict(d)
print(g2.name)          # fourcv
print(g2.wavelength)    # 1.5406
print(g2.sample.lattice.a)  # 5.431
```

## Save to a JSON file

JSON is the simplest option — it requires only the Python standard library:

```python
import json

# Save
with open("diffractometer.json", "w") as f:
    json.dump(g.to_dict(), f, indent=2)

# Restore
with open("diffractometer.json") as f:
    g2 = ahd.AdHocDiffractometer.from_dict(json.load(f))
```

## Save to a YAML file

YAML produces a more human-readable file.  It requires `pyyaml`
(not installed automatically — run `pip install pyyaml`):

```python
try:
    import yaml
except ImportError:
    raise ImportError("Install pyyaml: pip install pyyaml")

# Save
with open("diffractometer.yaml", "w") as f:
    yaml.dump(g.to_dict(), f, default_flow_style=False, sort_keys=False)

# Restore
with open("diffractometer.yaml") as f:
    g2 = ahd.AdHocDiffractometer.from_dict(yaml.safe_load(f))
```

## Verify a round-trip

Check that the restored geometry matches the original:

```python
import numpy as np

assert g2.name == g.name
assert g2.wavelength == g.wavelength
assert np.isclose(g2.sample.lattice.a, g.sample.lattice.a)
assert list(g2.sample.reflections.keys()) == list(g.sample.reflections.keys())
assert np.allclose(g2.sample.UB, g.sample.UB)
print("Round-trip verified.")
```

## Serialize individual objects

All major classes support `to_dict()` / `from_dict()` independently:

```python
from ad_hoc_diffractometer import Lattice, Sample

# Lattice
lat_dict = g.sample.lattice.to_dict()
lat2 = Lattice.from_dict(lat_dict)

# Sample (including reflections and UB)
sample_dict = g.sample.to_dict()
# sample2 = Sample.from_dict(sample_dict, parent=g)  # requires parent geometry
```

## ConstraintSet round-trip

Diffraction modes ({class}`~ad_hoc_diffractometer.mode.ConstraintSet`) are
serialised as part of the geometry dict.  They can also be inspected and
round-tripped independently:

```python
import json
from ad_hoc_diffractometer import ConstraintSet, BisectConstraint
from ad_hoc_diffractometer import SampleConstraint, DetectorConstraint
from ad_hoc_diffractometer import REQUIRED

# A custom psic mode
cs = ConstraintSet(
    [
        BisectConstraint("eta", "delta"),
        SampleConstraint("mu", 0.0),
        DetectorConstraint("nu", 0.0),
    ],
    computed=["eta", "chi", "phi", "delta"],
)

# Serialise
d = cs.to_dict()
print(json.dumps(d, indent=2))

# Restore
cs2 = ConstraintSet.from_dict(d)
assert cs2 == cs

# REQUIRED / OPTIONAL sentinels survive the round-trip
cs_extras = ConstraintSet(
    [BisectConstraint("eta", "delta"),
     SampleConstraint("mu", 0.0),
     SampleConstraint("chi", 0.0)],
    extras={"n_hat": REQUIRED, "psi": None},
)
d2 = cs_extras.to_dict()
cs_extras2 = ConstraintSet.from_dict(d2)
from ad_hoc_diffractometer import REQUIRED as REQ
assert cs_extras2.extras["n_hat"] is REQ  # sentinel restored
```

The geometry-level round-trip preserves all modes, the active mode name,
and cut points:

```python
import ad_hoc_diffractometer as ahd

g = ahd.presets.psic()
g.mode_name = "bisecting_vertical"

d = g.to_dict()
g2 = ahd.AdHocDiffractometer.from_dict(d)

assert set(g2.modes.keys()) == set(g.modes.keys())
assert g2.mode_name == "bisecting_vertical"
```

## Typical workflow

A common pattern is to save the configuration after alignment and restore
it at the start of the next session:

```python
import json
import ad_hoc_diffractometer as ahd

SESSION_FILE = "session_config.json"

def save_session(geometry):
    with open(SESSION_FILE, "w") as f:
        json.dump(geometry.to_dict(), f, indent=2)
    print(f"Session saved to {SESSION_FILE}")

def restore_session():
    with open(SESSION_FILE) as f:
        return ahd.AdHocDiffractometer.from_dict(json.load(f))

# End of session: save
g = ahd.presets.fourcv()
# ... do alignment ...
save_session(g)

# Start of next session: restore
g = restore_session()
ahd.pa(g)   # verify the restored configuration
```

## See also

- {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.to_dict`
- {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.from_dict`
- {class}`~ad_hoc_diffractometer.lattice.Lattice`
- {class}`~ad_hoc_diffractometer.sample.Sample`
- {func}`~ad_hoc_diffractometer.geometry.pa`
- {doc}`orient`
