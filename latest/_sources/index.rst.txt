.. _index:

============================
*Ad hoc* diffractometer
============================

**Version:** |release|

``ad_hoc_diffractometer`` is a **pure-Python** package for operating
multi-circle diffractometers in reciprocal space for X-ray and neutron
crystallography.  It is built around a key design principle: **any
multi-circle diffractometer geometry can be fully described by the caller**
— no geometry is hard-coded, and new geometries require no changes to the
package itself.

Its only runtime dependency beyond the Python Standard Library is
`NumPy <https://numpy.org>`_ — no scipy, sympy, or other scientific
libraries are required.

.. note::

   The package assumes **monochromatic radiation** throughout.  All
   diffraction calculations (Bragg angles, Q-vector magnitudes,
   forward/inverse problems) are performed at a fixed wavelength.

It provides:

- A class-based description of diffractometer stages (rotary axes) and
  their stacking order
- Predefined factory functions for standard synchrotron and laboratory
  diffractometer geometries (psic, fourcv, fourch, sixc, kappa families,
  zaxis, s2d2, fivec)
- Crystallographic lattice calculations (B matrix, reciprocal lattice)
- U and UB matrix computation from orienting reflections
- Forward diffraction calculations (hkl → motor angles), with
  diffraction modes controlling which stages are free, fixed, or coupled
- Wavelength, energy, d-spacing, and Q-vector conversions for X-ray
  and neutron sources
- Reciprocal-space operations: Q-vector magnitude, d-spacing, two-theta,
  and trajectory planning along arbitrary paths in reciprocal space

.. toctree::
   :hidden:

   install
   user_guide
   api
   changes
   direct-lattice
   problem1
   problem2

.. icons: https://fonts.google.com/icons

.. grid:: 2

   .. grid-item-card:: :material-regular:`install_desktop;3em` Get started
      :link: install
      :link-type: doc

      Install the package and verify the installation.

   .. grid-item-card:: :material-outlined:`menu_book;3em` User Guide
      :link: user_guide
      :link-type: doc

      Concepts, how-to guides, and a geometry reference.

   .. grid-item-card:: :material-regular:`api;3em` API Reference
      :link: api
      :link-type: doc

      Complete auto-generated reference for every public class,
      function, and constant.

   .. grid-item-card:: :material-outlined:`history;3em` Release Notes
      :link: changes
      :link-type: doc

      What changed in each release.

**Quick start**

.. code-block:: python

   import ad_hoc_diffractometer as ahd

   # Create a six-circle psic geometry and set the wavelength
   g = ahd.psic()
   g.wavelength = 1.0  # Å

   # Define the sample lattice (cubic silicon)
   g.sample.lattice = ahd.Lattice(a=5.431)

   # Show a summary of the geometry
   print(g.summary())

Background
----------

.. grid:: 3

   .. grid-item-card:: :material-outlined:`functions;3em` Direct Lattice
      :link: direct-lattice
      :link-type: doc

      Vector mathematics and lattice vector conventions in crystallography.

   .. grid-item-card:: :material-outlined:`science;3em` Case Study
      :link: problem1
      :link-type: doc

      The diffractometer problem that started this project.

   .. grid-item-card:: :material-outlined:`calculate;3em` You (1999) Geometry
      :link: problem2
      :link-type: doc

      Coordinate conventions and B/U/UB matrix derivation following
      You (1999) and Busing & Levy (1967).

About
-----

.. list-table::
   :stub-columns: 1

   * - Home
     - https://prjemian.github.io/ad_hoc_diffractometer/
   * - Source
     - https://github.com/prjemian/ad_hoc_diffractometer
   * - Version
     - |release|
   * - Published
     - |today|
   * - License
     - `CC-BY-4.0 <https://creativecommons.org/licenses/by/4.0/>`_
   * - Index
     - :ref:`genindex`
