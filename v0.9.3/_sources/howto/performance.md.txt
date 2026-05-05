(howto-performance)=
# Measure Forward/Inverse Performance

The {mod}`~ad_hoc_diffractometer.benchmark` module measures the throughput
(operations per second) and round-trip accuracy of
{meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.forward` and
{meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.inverse` across
all preset geometries and their declared diffraction modes.

## Run the benchmark

From the command line:

```bash
python -m ad_hoc_diffractometer.benchmark
```

From Python:

```python
from ad_hoc_diffractometer.benchmark import benchmark_all, benchmark_geometry

# All geometries, all modes
results = benchmark_all()

# Single geometry
results = benchmark_geometry("fourcv")
```

The output is a formatted table:

```text
geometry   mode                     status       fwd ops/s  inv ops/s  fwd/inv  round-trip err  solns
fourcv     bisecting                ok               3,241     18,502   0.1752        4.70e-13     11
fourcv     fixed_chi                ok               2,890     18,102   0.1596        3.40e-14     21
fourcv     fixed_psi                no_solutions      8,780          -       -               -      0
```

## Output columns

| Column | Description |
|--------|-------------|
| **geometry** | Preset geometry name |
| **mode** | Diffraction mode name |
| **status** | `ok`, `no_solutions`, `not_implemented`, or `error` |
| **fwd ops/s** | `forward()` operations per second |
| **inv ops/s** | `inverse()` operations per second |
| **fwd/inv** | `forward_ops / inverse_ops` — workstation-independent efficiency ratio. **Higher is better.** Since `inverse()` is a single direct computation, its speed characterizes the machine; the ratio measures algorithmic efficiency of the forward solver. |
| **round-trip err** | Maximum ‖hkl − inverse(forward(hkl))‖∞ |
| **solns** | Total solutions returned across all test reflections |

## Status values

- **ok** — `forward()` returned at least one solution; `inverse()` round-tripped
  successfully.
- **no_solutions** — `forward()` returned no solutions for any of the test
  reflections.  This is expected for some geometry+mode+reflection combinations
  (e.g. psi-constant modes where the target ψ does not match the natural value).
- **not_implemented** — the mode's solver is not yet implemented for this
  geometry.  `forward()` raises `NotImplementedError`.
- **error** — an unexpected exception occurred during setup or computation.
  The error message is included in the result dict.

## Controlling the benchmark

```python
from ad_hoc_diffractometer.benchmark import benchmark_all

# Fewer iterations (faster, noisier timing)
results = benchmark_all(n_iter=10)

# Custom reflections
results = benchmark_all(reflections=[(1, 0, 0), (1, 1, 1)])

# Suppress table output
results = benchmark_all(verbose=False)
```

## Factors that affect performance

**Number of solutions.**
Some geometry+mode combinations produce many valid solutions (e.g.
`double_diffraction` on a six-circle geometry can return dozens).  Each
solution requires inverse verification, so modes with more solutions are
slower per `forward()` call.

**Geometry complexity.**
Six-circle geometries have three free degrees of freedom after the Bragg
condition, requiring a three-dimensional search.  Four-circle geometries
have one free DOF and are correspondingly faster.

**Kappa geometries.**
The virtual-angle solver for kappa geometries uses a Newton–Raphson
iteration that adds overhead compared to direct Eulerian solutions.

**`inverse()` is fast.**
`inverse()` performs a single matrix solve (`UB⁻¹ @ Q_phi`) regardless
of geometry complexity — it does not enumerate solutions.  Typical
throughput is 5,000–20,000 ops/sec.

## Result dict structure

Each result is a Python dict with the following keys:

```python
{
    "geometry": "fourcv",
    "mode": "bisecting",
    "status": "ok",
    "forward_ops_per_sec": 3241.0,
    "inverse_ops_per_sec": 18502.0,
    "forward_inverse_ratio": 0.1752,
    "round_trip_max_error": 4.7e-13,
    "n_reflections": 5,
    "n_solutions": 11,
    "error_message": None,
}
```

## See also

- {func}`~ad_hoc_diffractometer.benchmark.benchmark_all`
- {func}`~ad_hoc_diffractometer.benchmark.benchmark_geometry`
- {func}`~ad_hoc_diffractometer.benchmark.benchmark_mode`
- {doc}`forward` — how forward calculations work
- {doc}`modes` — how to switch diffraction modes
