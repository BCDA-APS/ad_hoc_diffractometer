# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Declarative geometry definitions (issue #267).

This subpackage contains the **demonstration** YAML files that ship with
``ad_hoc_diffractometer``.  Each ``*.yml`` file is a declarative
description of one diffractometer geometry, parsed at import time by
:mod:`ad_hoc_diffractometer.geometry_loader` and registered with the
geometry registry in :mod:`ad_hoc_diffractometer.factories`.

The files here serve **two** purposes:

1. They populate the registry returned by
   :func:`ad_hoc_diffractometer.list_geometries` so callers can use
   ``ahd.make_geometry("fourcv")`` etc. out of the box.
2. They are **templates** that real beamlines should copy and adapt to
   describe their own instruments.  These shipped files are
   *demonstrations* of the schema — they are not authoritative
   descriptions of any production diffractometer.

The schema reference page (``docs/source/reference/declarative_geometry_schema.md``)
describes every accepted key, axis form, constraint type, and error
message.  The companion how-to (``docs/source/howto/custom_geometry.md``)
walks through writing a new geometry from one of these templates.

This ``__init__.py`` exists so that ``importlib.resources`` can locate
the directory as a regular Python package.  No symbols are exported
from this module.
"""
