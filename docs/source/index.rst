.. _index:

============================
*Ad hoc* diffractometer
============================

**Version:** |release|

``ad_hoc_diffractometer`` is a Python package that lets you describe any
multi-circle diffractometer geometry and perform X-ray/neutron
crystallography calculations.

See the :doc:`Quick Start guide <quick_start>` for a step-by-step
walkthrough building an Eulerian four-circle geometry — choosing a
coordinate basis, stage stacking, diffraction mode definition(s), and
running a forward calculation. Common geometries are provided as
examples.

.. less narrative on this page, it's the docs' home page
   TODO: This content could be moved to the user guide, perhaps.

   Features
   ------------

   ``ad_hoc_diffractometer`` handles the core calculations you need for diffractometer work:

   - **Geometry setup**: Describe your diffractometer using observable
   physical directions (vertical, longitudinal, transverse).
   - **Orientation calculations**: Compute orientation matrices from
   reflections, refine crystal lattice.
   - **Reciprocal space mapping**: Convert from rotation axes to reciprocal
   space coordinates.
   - **Diffractometer control**: Convert from reciprocal space coordinates
   to rotation axes.
   - **Mode Definitions**: You define which axes are free, fixed, or
   coupled when solving kinematics.

   It's your diffractometer
   --------------------------

   You get **full control over your setup** — whether you're using a
   standard four-circle geometry or something custom, the package adapts to
   you. No hard-coded configurations mean new geometries require zero
   changes to the code.

   Minimal Requirements
   -------------------------

   Only [**Python**](https://python.org) with its Standard Library and
   [**NumPy**](https://numpy.org). No scipy, sympy, or other scientific
   dependencies required.

   Use Cases
   ----------------

   - Simulating diffractometer behavior.
   - Real-time operations in reciprocal space during beamtime.
   - Planning experiments and trajectories before you run them.
   - Backend support for diffractometer control systems.
   - Creating visualizations of diffractometer geometry.

   .. note:: The package assumes **monochromatic radiation** throughout
      — all calculations are at a fixed wavelength.


.. toctree::
   :hidden:

   install
   user_guide
   api
   changes

.. icons: https://fonts.google.com/icons

.. grid:: 2

   .. grid-item-card:: :material-regular:`install_desktop;3em` Get started
      :link: install
      :link-type: doc

      Install the package and verify the installation.

   .. grid-item-card:: :material-outlined:`menu_book;3em` User Guide
      :link: user_guide
      :link-type: doc

      Concepts, how-to guides, and demonstration geometries.

   .. grid-item-card:: :material-regular:`api;3em` API Reference
      :link: api
      :link-type: doc

      Complete auto-generated reference for every public class,
      function, and constant.

   .. grid-item-card:: :material-outlined:`library_books;3em` References
      :link: references
      :link-type: doc

      All literature citations — geometry papers, physical constants,
      and numerical methods.

.. TODO: too many cards below
   Background
   ----------

   .. grid:: 3

      .. grid-item-card:: :material-outlined:`functions;3em` Direct Lattice
         :link: direct-lattice.html
         :link-type: url

         Vector mathematics and lattice vector conventions in crystallography.

      .. grid-item-card:: :material-outlined:`science;3em` Case Study
         :link: problem1.html
         :link-type: url

         The diffractometer problem that started this project.

      .. grid-item-card:: :material-outlined:`calculate;3em` Coordinate Convention & UB Matrix
         :link: problem2.html
         :link-type: url

         How a basis vector assignment leads to the B, U, and UB matrices, with a
         worked example using the convention described by H. You (1999.

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
