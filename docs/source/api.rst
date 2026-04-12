.. _api:

=============
API Reference
=============

The complete public API is auto-generated from docstrings by
`sphinx-autoapi <https://sphinx-autoapi.readthedocs.io/>`_.
The package is **pure Python** with a single runtime dependency:
`NumPy <https://numpy.org>`_.

.. icons: https://fonts.google.com/icons

.. grid:: 2

   .. grid-item-card:: :material-outlined:`place;3em` Foundation
      :link: autoapi/ad_hoc_diffractometer/constants/index
      :link-type: doc

      Core constants, rotation mathematics, and display formatting.

      :mod:`~ad_hoc_diffractometer.constants` ·
      :mod:`~ad_hoc_diffractometer.rotation` ·
      :mod:`~ad_hoc_diffractometer.display`

   .. grid-item-card:: :material-outlined:`settings_input_component;3em` Geometry Primitives
      :link: autoapi/ad_hoc_diffractometer/stage/index
      :link-type: doc

      Axis parsing and individual rotation stages.

      :mod:`~ad_hoc_diffractometer.axes` ·
      :mod:`~ad_hoc_diffractometer.stage`

   .. grid-item-card:: :material-outlined:`grain;3em` Crystallography
      :link: autoapi/ad_hoc_diffractometer/lattice/index
      :link-type: doc

      Lattice, reflections, sample, and radiation source.

      :mod:`~ad_hoc_diffractometer.lattice` ·
      :mod:`~ad_hoc_diffractometer.reflection` ·
      :mod:`~ad_hoc_diffractometer.sample` ·
      :mod:`~ad_hoc_diffractometer.radiation`

   .. grid-item-card:: :material-outlined:`precision_manufacturing;3em` Diffractometer
      :link: autoapi/ad_hoc_diffractometer/geometry/index
      :link-type: doc

      Geometry class, diffraction modes, and predefined geometry factories.

      :mod:`~ad_hoc_diffractometer.geometry` ·
      :mod:`~ad_hoc_diffractometer.mode` ·
      :mod:`~ad_hoc_diffractometer.factories`

   .. grid-item-card:: :material-outlined:`calculate;3em` Calculations
      :link: autoapi/ad_hoc_diffractometer/orientation/index
      :link-type: doc

      Orientation, forward problem, Q/energy conversions, trajectories,
      surface diffraction, and lattice refinement.

      :mod:`~ad_hoc_diffractometer.orientation` ·
      :mod:`~ad_hoc_diffractometer.forward` ·
      :mod:`~ad_hoc_diffractometer.conversions` ·
      :mod:`~ad_hoc_diffractometer.surface` ·
      :mod:`~ad_hoc_diffractometer.scan` ·
      :mod:`~ad_hoc_diffractometer.refinement`

   .. grid-item-card:: :material-regular:`api;3em` Full Module Index
      :link: autoapi/ad_hoc_diffractometer/index
      :link-type: doc

      Every public module, class, function, and constant,
      auto-generated from source docstrings.

Module Dependency Graph
-----------------------

The diagram below shows how the source modules depend on each other.
Arrows point from a module to the modules it imports.
Modules are grouped by layer.

.. graphviz::
   :caption: ad_hoc_diffractometer module dependencies
   :align: center

   digraph deps {
       graph [rankdir=TB, fontname="sans-serif", fontsize=11,
              bgcolor="transparent", splines=ortho, nodesep=0.4, ranksep=0.6]
       node  [shape=box, style="rounded,filled", fontname="sans-serif",
              fontsize=10, margin="0.15,0.08"]
       edge  [color="#555555", arrowsize=0.7]

       // ── Foundation (grey) ──────────────────────────────────
       subgraph cluster_foundation {
           label="Foundation"
           style=filled fillcolor="#f0f0f0" color="#aaaaaa"
           fontname="sans-serif" fontsize=11
           constants [label="constants" fillcolor="#e0e0e0"]
           rotation  [label="rotation"  fillcolor="#e0e0e0"]
           display   [label="display"   fillcolor="#e0e0e0"]
       }

       // ── Geometry primitives (blue) ──────────────────────────
       subgraph cluster_primitives {
           label="Geometry Primitives"
           style=filled fillcolor="#ddeeff" color="#6699cc"
           fontname="sans-serif" fontsize=11
           axes  [label="axes"  fillcolor="#c8dff8"]
           stage [label="stage" fillcolor="#c8dff8"]
       }

       // ── Crystallography (green) ─────────────────────────────
       subgraph cluster_crystal {
           label="Crystallography"
           style=filled fillcolor="#ddf0dd" color="#66aa66"
           fontname="sans-serif" fontsize=11
           lattice    [label="lattice"    fillcolor="#c0e8c0"]
           reflection [label="reflection" fillcolor="#c0e8c0"]
           sample     [label="sample"     fillcolor="#c0e8c0"]
           radiation  [label="radiation"  fillcolor="#c0e8c0"]
       }

       // ── Diffractometer (amber) ──────────────────────────────
       subgraph cluster_diffractometer {
           label="Diffractometer"
           style=filled fillcolor="#fff3cc" color="#ccaa00"
           fontname="sans-serif" fontsize=11
           mode      [label="mode"      fillcolor="#ffe999"]
           geometry  [label="geometry"  fillcolor="#ffe999"]
           factories [label="factories" fillcolor="#ffe999"]
       }

       // ── Calculations (orange) ───────────────────────────────
       subgraph cluster_calculations {
           label="Calculations"
           style=filled fillcolor="#ffe8d6" color="#cc7733"
           fontname="sans-serif" fontsize=11
           orientation [label="orientation" fillcolor="#ffd0b0"]
           forward     [label="forward"     fillcolor="#ffd0b0"]
           conversions [label="conversions" fillcolor="#ffd0b0"]
           surface     [label="surface"     fillcolor="#ffd0b0"]
           scan        [label="scan"        fillcolor="#ffd0b0"]
           refinement  [label="refinement"  fillcolor="#ffd0b0"]
       }

       // ── Edges ───────────────────────────────────────────────
       // Foundation
       axes -> constants
       stage -> rotation

       // Crystallography
       lattice    -> display
       reflection -> display
       sample     -> lattice
       sample     -> reflection

       // Diffractometer
       geometry -> axes         [style=dashed]
       geometry -> constants
       geometry -> display      [style=dashed]
       geometry -> forward      [style=dashed]
       geometry -> lattice      [style=dashed]
       geometry -> mode         [style=dashed]
       geometry -> orientation  [style=dashed]
       geometry -> radiation    [style=dashed]
       geometry -> reflection
       geometry -> rotation     [style=dashed]
       geometry -> sample       [style=dashed]
       geometry -> stage        [style=dashed]
       geometry -> surface      [style=dashed]
       factories -> axes
       factories -> constants
       factories -> geometry
       factories -> mode
       factories -> stage

       // Calculations
       orientation -> reflection
       orientation -> rotation
       orientation -> stage
       forward     -> geometry  [style=dashed]
       forward     -> mode
       forward     -> orientation
        conversions -> sample
       surface     -> geometry  [style=dashed]
       scan        -> geometry  [style=dashed]
       scan        -> rotation
       refinement  -> lattice
       refinement  -> orientation
       refinement  -> reflection
       refinement  -> rotation
   }

.. note::

   Dashed arrows indicate that the import is deferred (inside a function
   body or guarded by ``TYPE_CHECKING``) to avoid circular imports.

Full Module Reference
---------------------

.. toctree::
   :maxdepth: 1

   autoapi/ad_hoc_diffractometer/index
