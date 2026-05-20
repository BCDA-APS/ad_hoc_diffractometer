# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
Regression tests for issue #285.

Issue #285 reported that three horizontal-bisecting modes returned
zero forward solutions for the trigonal-rhombohedral reflection
``(1, 1, 0)`` on a single-parameter rhombohedral cell
(``a = b = c = 7.0``, ``alpha = beta = gamma = 72°``):

* ``fourch  / bisecting``
* ``psic    / bisecting_vertical``
* ``kappa6c / bisecting_horizontal``

Investigation summary
---------------------

Replaying the issue's self-contained reproducer (verbatim) in a
clean Python 3.14 / numpy 2.4 environment across released versions
shows that two of the three claims do not reproduce on any version
of the package, and the third was already fixed before the issue
was filed:

+------------------+--------------------+---------------------------+--------------------------------+
|     version      | fourch / bisecting | psic / bisecting_vertical | kappa6c / bisecting_horizontal |
+==================+====================+===========================+================================+
| v0.10.0          | 2 valid solutions  | 2 valid solutions         | 0 solutions (real failure)     |
+------------------+--------------------+---------------------------+--------------------------------+
| v0.10.1          | 2 valid solutions  | 2 valid solutions         | 0 solutions (real failure)     |
+------------------+--------------------+---------------------------+--------------------------------+
| v0.11.0 (PyPI)   | 2 valid solutions  | 2 valid solutions         | 2 valid solutions              |
+------------------+--------------------+---------------------------+--------------------------------+
| main (post-#284) | 2 valid solutions  | 2 valid solutions         | 2 valid solutions              |
+------------------+--------------------+---------------------------+--------------------------------+

Every returned solution round-trips correctly through ``inverse()``
back to ``(1, 1, 0)``.

The ``kappa6c bisecting_horizontal`` zero-solution failure on
v0.10.x was eliminated by the issue #280 composition-order overhaul
(merged in v0.11.0) — most likely the B-matrix correction in phase 3,
which changed the reciprocal-lattice geometry for non-orthogonal
cells.  The ``fourch bisecting`` and ``psic bisecting_vertical``
claims in the issue body could not be reproduced on any released
version; the most plausible explanation is that the issue's
"self-contained reproducer" was synthesized from a downstream
cross-validation harness that uses a *bootstrapped* UB (built from
``calc_UB`` of two seed reflections at chosen positions) rather than
``ub_identity``, and the simplification to ``ub_identity`` lost the
specific UB orientation that produced the zero-solution branch
selection.

No code change was required to close issue #285.  These tests pin
down the working behavior so any future regression on the rhombo-
hedral path will be caught immediately.
"""

from __future__ import annotations

from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd

# Trigonal-rhombohedral single-parameter cell from the issue's
# reproducer block.
TRIGONAL_RHOMBOHEDRAL = dict(a=7.0, b=7.0, c=7.0, alpha=72.0, beta=72.0, gamma=72.0)
WAVELENGTH = 1.54


# ---------------------------------------------------------------------------
# Issue #285 reproducer: the three (geometry, mode) cases solve.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "geom_name, mode_name, context",
    [
        pytest.param("fourch", "bisecting", does_not_raise(), id="fourch-bisecting"),
        pytest.param(
            "psic",
            "bisecting_vertical",
            does_not_raise(),
            id="psic-bisecting_vertical",
        ),
        pytest.param(
            "kappa6c",
            "bisecting_horizontal",
            does_not_raise(),
            id="kappa6c-bisecting_horizontal",
        ),
    ],
)
def test_horizontal_bisecting_solves_trigonal_rhombohedral_110(
    geom_name, mode_name, context
):
    """``forward(1, 1, 0)`` returns at least one solution that round-trips.

    Direct copy of the reproducer block in the issue body.  Pre-fix
    (v0.10.x) ``kappa6c bisecting_horizontal`` returned ``n_sols=0``
    here; the other two cases returned valid solutions on every
    released version.
    """
    g = ahd.make_geometry(geom_name)
    g.sample.lattice = ahd.Lattice(**TRIGONAL_RHOMBOHEDRAL)
    g.wavelength = WAVELENGTH
    ahd.ub_identity(g.sample)
    g.mode_name = mode_name
    with context:
        sols = g.forward(1, 1, 0)
        assert len(sols) >= 1, (
            f"{geom_name} / {mode_name}: expected at least one solution "
            f"for trigonal-rhombohedral (1, 1, 0), got {len(sols)}."
        )
        for sol in sols:
            hkl_back = g.inverse(sol)
            np.testing.assert_allclose(hkl_back, (1.0, 1.0, 0.0), atol=1e-6)
