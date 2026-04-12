# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
constants.py — shared Cartesian basis vectors.

These are the internal numpy representations used throughout the package.
The caller-facing notation (+x, -z, etc.) is handled in axes.py.

Default coordinate convention (You 1999):

- XHAT (+x) — vertical (out of the floor)
- YHAT (+y) — longitudinal (along the beam, toward equipment)
- ZHAT (+z) — lateral (to our left when facing equipment)
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)

logger = logging.getLogger(__name__)

#: Unit vector along the +x axis.
#: In the You (1999) convention this is the **vertical** direction (out of the floor).
#: In the Busing & Levy (1967) convention this is the **lateral** direction.
XHAT = np.array([1.0, 0.0, 0.0])

#: Unit vector along the +y axis.
#: **Longitudinal** direction (along the beam, toward the equipment) in both
#: the You (1999) and Busing & Levy (1967) conventions.
YHAT = np.array([0.0, 1.0, 0.0])

#: Unit vector along the +z axis.
#: In the You (1999) convention this is the **lateral** direction (to our left
#: when facing the equipment).
#: In the Busing & Levy (1967) convention this is the **vertical** direction.
ZHAT = np.array([0.0, 0.0, 1.0])
