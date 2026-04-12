(howto-refine-lattice)=
# Refine the Lattice Constants

This guide shows how to refine unit-cell parameters from a set of measured
Bragg peak positions.  Two methods are available:

| Method | Function | Algorithm | When to use |
|---|---|---|---|
| Busing & Levy (1967) | `refine_lattice_bl1967()` | Iterative least-squares with analytic Jacobian | Good starting point; can refine cell *and* orientation simultaneously |
| Nelder-Mead simplex | `refine_lattice_simplex()` | Derivative-free simplex | Far from solution; irregular residual surface; orientation held fixed |

Both methods support constrained refinement (crystal-system symmetry enforced)
and unconstrained refinement (all six parameters free).

## Prerequisites

- A geometry with a wavelength set
- A UB matrix set on the sample (see {doc}`orient`)
- At least three measured reflections with accurate motor angles

More reflections give a more reliable result; six or more are recommended.

## Set up the geometry and collect reflections

```python
import ad_hoc_diffractometer as ahd

g = ahd.fourcv()
g.wavelength = 1.5498   # Å
g.sample.lattice = ahd.Lattice(a=4.785, c=12.991, gamma=120.0)  # sapphire

# Set UB from two orienting reflections (see howto/orient)
ahd.ub_from_two_reflections_bl1967(
    g.sample,
    r1=g.reflections["r1"],
    r2=g.reflections["r2"],
)

# Collect additional measured reflections for refinement
measured = [
    ("r1",  (0, 0,  6), {"ttheta": 41.937, "omega": 20.931, "chi": 90.032, "phi":  0.0}),
    ("r2",  (1, 0,  4), {"ttheta": 33.508, "omega": 16.714, "chi": 76.074, "phi":  0.0}),
    ("r3",  (0, 1,  2), {"ttheta": 27.823, "omega": 13.882, "chi": 57.587, "phi": 90.0}),
    ("r4",  (1, 1,  0), {"ttheta": 36.874, "omega": 18.449, "chi": 76.989, "phi": 45.0}),
    ("r5",  (0, 0, 12), {"ttheta": 91.371, "omega": 45.657, "chi": 90.061, "phi":  0.0}),
    ("r6",  (2, 0,  0), {"ttheta": 40.181, "omega": 20.053, "chi": 90.034, "phi": 30.0}),
]
for name, hkl, angles in measured:
    g.reflections.add(name, hkl=hkl, angles=angles)

rlist = [g.reflections[name] for name, *_ in measured]
```

---

## Method 1 — Busing & Levy least-squares

Iteratively solves the linearised normal equations to minimise the RMS
misfit between observed and calculated phi-frame scattering vectors.
Can refine cell parameters and orientation (U matrix) simultaneously.

### Step 1 — Constrained (crystal-system symmetry enforced)

```python
result1 = ahd.refine_lattice_bl1967(
    g.sample,
    rlist,
    refine_cell=True,
    refine_orientation=True,  # also refine U
    refine_all=False,         # enforce hexagonal constraints
)

lat1 = result1["lattice"]
print(f"a = {lat1.a:.6f} Å   c = {lat1.c:.6f} Å   gamma = {lat1.gamma:.4f}°")
print(f"RMS misfit : {result1['rms']:.3e}")
print(f"Converged  : {result1['converged']}  ({result1['n_iter']} iterations)")
```

### Step 2 — Unconstrained (all six parameters free)

```python
g.sample.lattice = result1["lattice"]

result2 = ahd.refine_lattice_bl1967(
    g.sample,
    rlist,
    refine_all=True,   # all six parameters independent
)

lat2 = result2["lattice"]
print(f"a = {lat2.a:.6f}  b = {lat2.b:.6f}  c = {lat2.c:.6f}")
print(f"α = {lat2.alpha:.4f}°  β = {lat2.beta:.4f}°  γ = {lat2.gamma:.4f}°")
print(f"RMS misfit : {result2['rms']:.3e}")
```

---

## Method 2 — Nelder-Mead simplex

Derivative-free global minimisation.  More robust when the starting point
is far from the solution, but slower and does not refine U by default.
Use this as a first pass when the BL1967 method fails to converge.

### Step 1 — Constrained

```python
result1 = ahd.refine_lattice_simplex(
    g.sample,
    rlist,
    refine_all=False,   # enforce crystal-system constraints
)

lat1 = result1["lattice"]
print(f"a = {lat1.a:.6f} Å   c = {lat1.c:.6f} Å   gamma = {lat1.gamma:.4f}°")
print(f"RMS misfit : {result1['rms']:.3e}")
print(f"Converged  : {result1['converged']}  ({result1['n_iter']} evaluations)")
```

### Step 2 — Unconstrained

```python
g.sample.lattice = result1["lattice"]

result2 = ahd.refine_lattice_simplex(
    g.sample,
    rlist,
    refine_all=True,
)

lat2 = result2["lattice"]
print(f"a = {lat2.a:.6f}  b = {lat2.b:.6f}  c = {lat2.c:.6f}")
print(f"α = {lat2.alpha:.4f}°  β = {lat2.beta:.4f}°  γ = {lat2.gamma:.4f}°")
print(f"RMS misfit : {result2['rms']:.3e}")
```

---

## Interpreting the results

| Quantity | Meaning |
|---|---|
| `rms` | RMS misfit between observed and predicted Q-vectors (Å⁻¹) |
| `converged` | `True` if the method converged within `max_iter` |
| `n_iter` | Iterations (BL1967) or function evaluations (simplex) used |
| `residuals` | Per-reflection misfit vector |

**Constrained vs unconstrained:**
Always start with constrained refinement (`refine_all=False`).  It uses
fewer free parameters, converges faster, and enforces the known crystal
symmetry.  Use unconstrained refinement (`refine_all=True`) as a
diagnostic: if b deviates significantly from a, or γ from 120°, this
may indicate a real distortion, an indexing error, or miscalibrated angles.

**Choosing a method:**

- Start with **BL1967** — it is faster and also refines the orientation.
- If BL1967 does not converge (poor starting lattice, highly distorted
  cell), use **simplex** first to get close, then hand off to BL1967 for
  final polishing.

## See also

- {func}`~ad_hoc_diffractometer.refine_lattice_bl1967`
- {func}`~ad_hoc_diffractometer.refine_lattice_simplex`
- {class}`~ad_hoc_diffractometer.lattice.Lattice`
- {doc}`orient`
- {doc}`lattice`
