(howto-wavelength)=
# Set Wavelength / Energy

This guide shows how to set the radiation wavelength on a geometry and
convert between wavelength, energy, d-spacing, and Q-vector magnitude.

## Set the wavelength

```python
import ad_hoc_diffractometer as ahd

g = ahd.psic()
g.wavelength = 1.0  # Å
print(g.wavelength)  # 1.0
```

The wavelength must be positive (in Å).  Setting `None` clears it:

```python
g.wavelength = None  # unset
```

## Convert wavelength ↔ energy

```python
# Wavelength (Å) → energy (keV)
E = ahd.wavelength_to_energy(1.5406)   # → 8.048 keV  (Cu Kα)

# Energy (keV) → wavelength (Å)
lam = ahd.energy_to_wavelength(12.398)  # → 1.0 Å
```

## Convert wavelength ↔ wavenumber

```python
k = ahd.wavelength_to_wavenumber(1.0)   # → 2π Å⁻¹
lam = ahd.wavenumber_to_wavelength(k)   # → 1.0 Å
```

## Convert d-spacing ↔ Q magnitude

```python
Q = ahd.d_to_Q_mag(2.0)        # d = 2.0 Å → Q = π Å⁻¹
d = ahd.Q_to_d(Q)              # back to d
```

## Convert d-spacing / Q ↔ two-theta

```python
lam = 1.5406  # Å, Cu Kα
tth = ahd.d_to_two_theta(2.0, lam)           # d → 2θ (degrees)
tth = ahd.Q_mag_to_two_theta(ahd.d_to_Q_mag(2.0), lam)  # Q → 2θ
d   = ahd.two_theta_to_d(tth, lam)           # 2θ → d
Q   = ahd.two_theta_to_Q_mag(tth, lam)       # 2θ → Q
```

## Neutron sources

Neutron energy is expressed in meV and wavelength in Å:

```python
lam = ahd.neutron_energy_to_wavelength(25.3)  # meV → Å  (thermal neutron)
E   = ahd.neutron_wavelength_to_energy(1.8)   # Å → meV
```

## Common X-ray laboratory lines

A dictionary of standard X-ray lines is available:

```python
print(ahd.XRAY_LINES["Cu Ka1"])   # 1.5406 Å
print(ahd.XRAY_LINES["Mo Ka1"])   # 0.7093 Å
```

## See also

- {func}`~ad_hoc_diffractometer.wavelength_to_energy`
- {func}`~ad_hoc_diffractometer.energy_to_wavelength`
- {func}`~ad_hoc_diffractometer.d_to_Q_mag`
- {func}`~ad_hoc_diffractometer.Q_to_d`
- {mod}`~ad_hoc_diffractometer.radiation`
- {mod}`~ad_hoc_diffractometer.conversions`
