(geometries)=
Geometry Reference
==================

One page per predefined diffractometer geometry.  Each factory function
returns a fully configured
:class:`~ad_hoc_diffractometer.geometry.AdHocDiffractometer` instance.

.. toctree::
   :hidden:

   fourcv
   fourch
   psic
   sixc
   kappa4cv
   kappa4ch
   kappa6c
   zaxis
   s2d2
   fivec

.. icons: https://fonts.google.com/icons

.. grid:: 2

   .. grid-item-card:: :material-outlined:`rotate_right;3em` fourcv
      :link: fourcv
      :link-type: doc

      Four-circle Eulerian, vertical scattering plane (synchrotron).

   .. grid-item-card:: :material-outlined:`rotate_right;3em` fourch
      :link: fourch
      :link-type: doc

      Four-circle Eulerian, horizontal scattering plane (laboratory).

   .. grid-item-card:: :material-outlined:`rotate_right;3em` psic
      :link: psic
      :link-type: doc

      Six-circle 4S+2D (You 1999), vertical scattering plane.

   .. grid-item-card:: :material-outlined:`rotate_right;3em` sixc
      :link: sixc
      :link-type: doc

      Six-circle surface diffractometer (Lohmeier & Vlieg 1993).

   .. grid-item-card:: :material-outlined:`rotate_right;3em` kappa4cv
      :link: kappa4cv
      :link-type: doc

      Kappa four-circle, vertical scattering plane (synchrotron).

   .. grid-item-card:: :material-outlined:`rotate_right;3em` kappa4ch
      :link: kappa4ch
      :link-type: doc

      Kappa four-circle, horizontal scattering plane (laboratory).

   .. grid-item-card:: :material-outlined:`rotate_right;3em` kappa6c
      :link: kappa6c
      :link-type: doc

      Six-circle kappa, psic-style outer axes (synchrotron).

   .. grid-item-card:: :material-outlined:`rotate_right;3em` zaxis
      :link: zaxis
      :link-type: doc

      Z-axis four-circle, surface diffraction (Bloch 1985).

   .. grid-item-card:: :material-outlined:`rotate_right;3em` s2d2
      :link: s2d2
      :link-type: doc

      S2D2 four-circle, fully decoupled sample and detector axes.

   .. grid-item-card:: :material-outlined:`rotate_right;3em` fivec
      :link: fivec
      :link-type: doc

      Five-circle: fourcv on a vertical mu base (Vlieg et al. 1987).
