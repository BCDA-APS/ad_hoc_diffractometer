(references)=
# References

Literature citations used throughout `ad_hoc_diffractometer`.
Geometries, algorithms, and conventions are traced to their primary sources.

---

## Diffractometer geometries

**Busing & Levy (1967)**
: W.R. Busing and H.A. Levy.
  *Angle calculations for 3- and 4-circle X-ray and neutron diffractometers.*
  Acta Crystallographica **22**, 457–464 (1967).
  DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970)

  Foundational reference for the four-circle geometry, B matrix, U matrix,
  and UB matrix.  Defines the orientation refinement least-squares procedure.
  Used by: {func}`~ad_hoc_diffractometer.fourcv`,
  {func}`~ad_hoc_diffractometer.fourch`,
  {func}`~ad_hoc_diffractometer.kappa4cv`,
  {func}`~ad_hoc_diffractometer.kappa4ch`.

**Bloch (1985)**
: J.M. Bloch.
  *Angle and distance calculations for X-ray diffraction with the Z-axis
  geometry.*
  Journal of Applied Crystallography **18**, 33–36 (1985).
  DOI: [10.1107/S0021889885009858](https://doi.org/10.1107/S0021889885009858)

  Defines the Z-axis diffractometer geometry.
  Used by: {func}`~ad_hoc_diffractometer.zaxis`.

**Vlieg et al. (1987)**
: E. Vlieg, A.E.M.J. Fischer, J.F. van der Veen, B.N. Dev, and G. Materlik.
  *Surface X-ray diffraction: a study of relaxation in the Cu(110) system.*
  Journal of Applied Crystallography **20**, 330–337 (1987).
  DOI: [10.1107/S0021889887087266](https://doi.org/10.1107/S0021889887087266)

  Defines the five-circle geometry.
  Used by: {func}`~ad_hoc_diffractometer.fivec`.

**Lohmeier & Vlieg (1993)**
: M. Lohmeier and E. Vlieg.
  *Angle calculations for a six-circle surface X-ray diffractometer.*
  Journal of Applied Crystallography **26**, 706–716 (1993).
  DOI: [10.1107/S0021889893006198](https://doi.org/10.1107/S0021889893006198)

  Defines the six-circle surface diffractometer geometry.
  Used by: {func}`~ad_hoc_diffractometer.sixc`.

**Evans-Lutterodt & Tang (1995)**
: K.W. Evans-Lutterodt and M.-T. Tang.
  *Angle calculations for a '2+2' surface X-ray diffractometer.*
  Journal of Applied Crystallography **28**, 318–326 (1995).
  DOI: [10.1107/S0021889895001063](https://doi.org/10.1107/S0021889895001063)

  Defines the S2D2 (2+2) diffractometer geometry.
  Used by: {func}`~ad_hoc_diffractometer.s2d2`.

**You (1999)**
: H. You.
  *Angle calculations for a '4S+2D' six-circle diffractometer.*
  Journal of Applied Crystallography **32**, 614–623 (1999).
  DOI: [10.1107/S0021889899001223](https://doi.org/10.1107/S0021889899001223)

  Defines the psic (4S+2D) six-circle geometry; axis sign conventions
  (mixed handedness); ψ angle definitions (eqs. 10–11).
  Used by: {func}`~ad_hoc_diffractometer.psic`,
  {func}`~ad_hoc_diffractometer.kappa6c`.

**ITC Vol. C §2.2.6 (2006)**
: International Tables for Crystallography, Volume C, Section 2.2.6.
  *Single-crystal X-ray techniques.*
  DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)

  Confirms the kappa 50° tilt convention and normal-beam equatorial geometry.
  Used by: {func}`~ad_hoc_diffractometer.kappa4cv`,
  {func}`~ad_hoc_diffractometer.kappa4ch`,
  {func}`~ad_hoc_diffractometer.kappa6c`.

**Walko (2016)**
: D.A. Walko.
  *X-ray diffractometers.*
  Reference Module in Materials Science and Materials Engineering,
  Elsevier (2016).

  Comprehensive survey of diffractometer geometry designations (S/D system);
  kappa convention; zaxis, s2d2 geometries.
  Used throughout the geometry factory descriptions.

---

## Physical constants

**CODATA 2022 / 2019 SI**
: NIST CODATA 2022 recommended values.
  BIPM SI Brochure, 9th edition (2019).

  - $hc = 12.398\,419\,843\,320\,026\,\text{keV·Å}$ — exact (h and c are
    defined constants since the 2019 SI redefinition).
  - $h^2/(2m_n) = 81.804\,210\,235\,2\,\text{meV·Å}^2$ — from CODATA 2022
    neutron mass $m_n$.

  See {data}`~ad_hoc_diffractometer.HC_KEV_ANGSTROM` and
  {data}`~ad_hoc_diffractometer.NEUTRON_MEV_ANGSTROM2`.

---

## Numerical methods

**Nelder & Mead (1965)**
: J.A. Nelder and R. Mead.
  *A simplex method for function minimization.*
  The Computer Journal **7**(4), 308–313 (1965).
  DOI: [10.1093/comjnl/7.4.308](https://doi.org/10.1093/comjnl/7.4.308)

  Derivative-free simplex optimisation algorithm used by
  {func}`~ad_hoc_diffractometer.refine_lattice_simplex`.

---

## APS alignment session

**Walko (2020)**
: D.A. Walko, private communication (December 2020).
  Crystal alignment session at APS beamline 7-ID-C using a sapphire sample.
  Documented in {doc}`howto/fourcv_alignment_howto`.
