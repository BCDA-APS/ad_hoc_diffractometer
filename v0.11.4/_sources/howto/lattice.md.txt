(howto-lattice)=
# Define the Sample Lattice

This guide shows how to specify unit-cell parameters for any crystal
system and inspect the resulting B matrix and reciprocal lattice.

## Create a Lattice

Pass the minimum parameters needed for your crystal system:

```python
import ad_hoc_diffractometer as ahd

# Cubic: a only
silicon = ahd.Lattice(a=5.431)

# Tetragonal: a and c
tio2 = ahd.Lattice(a=4.594, c=2.959)

# Hexagonal: a, c, and gamma=120 (must be explicit)
sapphire = ahd.Lattice(a=4.785, c=12.991, gamma=120)

# Trigonal: a and alpha ≠ 90
calcite = ahd.Lattice(a=6.375, alpha=46.1)

# Orthorhombic: a, b, c
aragonite = ahd.Lattice(a=4.962, b=7.967, c=5.741)

# Monoclinic: a, b, c, beta
gypsum = ahd.Lattice(a=5.68, b=15.18, c=6.29, beta=118.4)

# Triclinic: a, b, c, alpha, beta, gamma
albite = ahd.Lattice(a=8.144, b=12.787, c=7.160,
                     alpha=94.33, beta=116.57, gamma=87.65)
```

The crystal system is deduced automatically from the supplied parameters:

```python
print(silicon.crystal_system)  # 'cubic'
print(sapphire.crystal_system) # 'hexagonal'
```

## Crystal system summary

| System | Minimum parameters required |
|---|---|
| Cubic | `a` |
| Tetragonal | `a`, `c` |
| Hexagonal | `a`, `c`, `gamma=120` (must be explicit) |
| Trigonal | `a`, `alpha` (alpha ≠ 90) |
| Orthorhombic | `a`, `b`, `c` |
| Monoclinic | `a`, `b`, `c`, `beta` |
| Triclinic | `a`, `b`, `c`, `alpha`, `beta`, `gamma` |

## Assign the lattice to a sample

```python
g = ahd.make_geometry("fourcv")
g.sample.lattice = silicon
```

## Inspect the B matrix

```python
import numpy as np

B = g.sample.lattice.B
print(B)
```

The B matrix satisfies $|\mathbf{B}\,\mathbf{h}| = 2\pi / d_{hkl}$:

```python
h = [0, 0, 1]
Q_mag = np.linalg.norm(B @ h)
d_hkl = 2 * np.pi / Q_mag
print(f"d(001) = {d_hkl:.4f} Å")
```

## Inspect the reciprocal lattice

```python
b1, b2, b3 = g.sample.lattice.reciprocal_lattice_vectors
print(f"b1 = {b1}")
```

## Compute d-spacing and Q for a reflection

```python
d   = ahd.hkl_to_d(g, 1, 1, 0)
Q   = ahd.hkl_to_Q(g, 1, 1, 0)
tth = ahd.conversions.hkl_to_two_theta(g, 1, 1, 0)
print(f"d(110) = {d:.4f} Å, 2θ = {tth:.3f}°")
```

## See also

- {class}`~ad_hoc_diffractometer.lattice.Lattice`
- {func}`~ad_hoc_diffractometer.conversions.hkl_to_d`
- {func}`~ad_hoc_diffractometer.conversions.hkl_to_Q`
- {func}`~ad_hoc_diffractometer.conversions.hkl_to_two_theta`
