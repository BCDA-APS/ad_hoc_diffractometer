(references)=
# References

Literature citations used throughout `ad_hoc_diffractometer`.
Geometries, algorithms, and conventions are traced to their primary sources.

---

## Diffractometer geometries

**Arndt & Willis (1966)**
: U.W. Arndt and B.T.M. Willis.
  *Single Crystal Diffractometry.*
  Cambridge University Press (1966); online edition 21 May 2010,
  ISBN 9780511735622.
  DOI: [10.1017/CBO9780511735622](https://doi.org/10.1017/CBO9780511735622)

  Early monograph on the mechanical design and use of single-crystal
  diffractometers, including the chapter
  *[Design of diffractometers](https://www.cambridge.org/core/books/abs/single-crystal-diffractometry/design-of-diffractometers/C522144B857706F4C973B46260FFB3D2)*.

**Busing & Levy (1967)**
: W.R. Busing and H.A. Levy.
  *Angle calculations for 3- and 4-circle X-ray and neutron diffractometers.*
  Acta Crystallographica **22**, 457–464 (1967).
  DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970)

  Foundational reference for the four-circle geometry, B matrix, U matrix,
  and UB matrix.  Defines the orientation refinement least-squares procedure.
  Used by: {func}`~ad_hoc_diffractometer.presets.fourcv`,
  {func}`~ad_hoc_diffractometer.presets.fourch`,
  {func}`~ad_hoc_diffractometer.presets.kappa4cv`,
  {func}`~ad_hoc_diffractometer.presets.kappa4ch`.

**Wyckoff (1985)**
: H.W. Wyckoff.
  *Diffractometry.*
  Methods in Enzymology **114**, 330–386 (1985).
  DOI: [10.1016/0076-6879(85)14026-7](https://doi.org/10.1016/0076-6879(85)14026-7)

  Figure 2(b) on p. 334 — the canonical Enraf-Nonius kappa
  diffractometer with explicit X/Y/Z coordinate axes and κ-rotation
  NEG/POS sense; the schematic cited by ITC Vol. C §2.2.6 as the
  reference for the kappa goniostat geometry.
  Used by: {func}`~ad_hoc_diffractometer.presets.kappa4ch`.

**Bloch (1985)**
: J.M. Bloch.
  *Angle and distance calculations for X-ray diffraction with the Z-axis
  geometry.*
  Journal of Applied Crystallography **18**, 33–36 (1985).
  DOI: [10.1107/S0021889885009858](https://doi.org/10.1107/S0021889885009858)

  Defines the Z-axis diffractometer geometry.
  Used by: {func}`~ad_hoc_diffractometer.presets.zaxis`.

**Vlieg et al. (1987)**
: E. Vlieg, A.E.M.J. Fischer, J.F. van der Veen, B.N. Dev, and G. Materlik.
  *Surface X-ray diffraction: a study of relaxation in the Cu(110) system.*
  Journal of Applied Crystallography **20**, 330–337 (1987).
  DOI: [10.1107/S0021889887087266](https://doi.org/10.1107/S0021889887087266)

  Defines the five-circle geometry.
  Used by: {func}`~ad_hoc_diffractometer.presets.fivec`.

**Lohmeier & Vlieg (1993)**
: M. Lohmeier and E. Vlieg.
  *Angle calculations for a six-circle surface X-ray diffractometer.*
  Journal of Applied Crystallography **26**, 706–716 (1993).
  DOI: [10.1107/S0021889893006198](https://doi.org/10.1107/S0021889893006198)

  Defines the six-circle surface diffractometer geometry.
  Used by: {func}`~ad_hoc_diffractometer.presets.sixc`.

**Evans-Lutterodt & Tang (1995)**
: K.W. Evans-Lutterodt and M.-T. Tang.
  *Angle calculations for a '2+2' surface X-ray diffractometer.*
  Journal of Applied Crystallography **28**, 318–326 (1995).
  DOI: [10.1107/S0021889895001063](https://doi.org/10.1107/S0021889895001063)

  Defines the S2D2 (2+2) diffractometer geometry.
  Used by: {func}`~ad_hoc_diffractometer.presets.s2d2`.

**Paciorek, Meyer & Chapuis (1999)**
: W.A. Paciorek, M. Meyer, and G. Chapuis.
  *On the geometry of a modern imaging diffractometer.*
  Acta Crystallographica A **55**, 543–557 (1999).
  DOI: [10.1107/S0108767399000951](https://doi.org/10.1107/S0108767399000951)

  Mathematical description of a four-circle κ-goniometer
  (KM4CCD) with explicit coordinate frame (Fig. 1) and the
  κ-axis tilt parameters χ, α (Fig. 2).  Cited by Sønsteby
  *et al.* (2013) as the basic mathematics of the
  six-axis κ instrument.

**You (1999)**
: H. You.
  *Angle calculations for a '4S+2D' six-circle diffractometer.*
  Journal of Applied Crystallography **32**, 614–623 (1999).
  DOI: [10.1107/S0021889899001223](https://doi.org/10.1107/S0021889899001223)

  Defines the psic (4S+2D) six-circle geometry; axis sign conventions
  (mixed handedness); ψ angle definitions (eqs. 10–11).
   Used by: {func}`~ad_hoc_diffractometer.presets.psic`,
   {func}`~ad_hoc_diffractometer.presets.kappa6c`.

**ITC Vol. C §2.2.6 (2006)**
: International Tables for Crystallography, Volume C, Section 2.2.6,
  p. 36.
  *Single-crystal X-ray techniques.*
  DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)

  Confirms the kappa 50° tilt convention; cites Wyckoff (1985, p. 334)
  for the schematic picture of the kappa goniostat.  §2.2.6.2
  documents the standard sign convention (right-handed for ω/χ/φ;
  left-handed for 2θ in Hamilton's choice) — note that the presets
  shipped here follow Walko's left-handed convention for ω/φ/2θ; see
  the {doc}`/concepts` page for the handedness discussion.
  Used by: {func}`~ad_hoc_diffractometer.presets.kappa4cv`,
  {func}`~ad_hoc_diffractometer.presets.kappa4ch`,
  {func}`~ad_hoc_diffractometer.presets.kappa6c`.

**Thorkildsen, Larsen & Beukes (2006)**
: G. Thorkildsen, H.B. Larsen, and J.A. Beukes.
  *Angle calculations for a three-circle goniostat.*
  Journal of Applied Crystallography **39**, 151–157 (2006).
  DOI: [10.1107/S0021889805041877](https://doi.org/10.1107/S0021889805041877)

  Vector formulation of the diffractometer angle-calculation problem,
  applicable to arbitrary goniostat designs (Eulerian, κ, and
  generalisations).  Table 1 + equation (3) give the canonical
  κ-goniostat axes used by ``kappa4cv``, ``kappa4ch``, and
  ``kappa6c``; §3 last paragraph explicitly notes the extension to
  additional rotation axes.
  Used by: {func}`~ad_hoc_diffractometer.presets.kappa4cv`,
  {func}`~ad_hoc_diffractometer.presets.kappa4ch`,
  {func}`~ad_hoc_diffractometer.presets.kappa6c`.

**Sønsteby et al. (2013)**
: H.H. Sønsteby, D. Chernyshov, M. Getz, O. Nilsen, and H. Fjellvåg.
  *On the application of a single-crystal κ-diffractometer and a CCD
  area detector for studies of thin films.*
  Journal of Synchrotron Radiation **20**, 644–647 (2013).
  DOI: [10.1107/S0909049513009102](https://doi.org/10.1107/S0909049513009102)

  Six-axis κ-diffractometer (KUMA6) at Swiss-Norwegian Beam Lines
  BM01A at the European Synchrotron Radiation Facility.  Cites
  Paciorek (1999) for the basic mathematics and Thorkildsen
  *et al.* (2006) for the angular-calculation method used by the
  CrysAlis control software.  The reference instrument for the
  ``kappa6c`` preset.
  Used by: {func}`~ad_hoc_diffractometer.presets.kappa6c`.

**Walko (2016)**
: D.A. Walko.
  *Multicircle Diffractometry Methods.*
  Reference Module in Materials Science and Materials Engineering,
  Elsevier (2016).
  DOI: [10.1016/B978-0-12-803581-8.01215-7](https://doi.org/10.1016/B978-0-12-803581-8.01215-7)

  Comprehensive survey of diffractometer geometry designations (S/D
  system); kappa convention; zaxis, s2d2, fivec geometries.  Figure 3
  shows the kappa diffractometer in vertical-scattering layout (the
  ``kappa4cv`` reference).
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

  See {data}`~ad_hoc_diffractometer.radiation.HC_KEV_ANGSTROM` and
  {data}`~ad_hoc_diffractometer.radiation.NEUTRON_MEV_ANGSTROM2`.

---

## Numerical methods

**Nelder & Mead (1965)**
: J.A. Nelder and R. Mead.
  *A simplex method for function minimization.*
  The Computer Journal **7**(4), 308–313 (1965).
  DOI: [10.1093/comjnl/7.4.308](https://doi.org/10.1093/comjnl/7.4.308)

  Derivative-free simplex optimisation algorithm used by
  {func}`~ad_hoc_diffractometer.refinement.refine_lattice_simplex`.

---

## APS alignment session

**Walko (2020)**
: D.A. Walko, private communication (December 2020).
  Crystal alignment session at APS beamline 7-ID-C using a sapphire sample.
  Documented in {doc}`howto/fourcv_alignment_howto`.
