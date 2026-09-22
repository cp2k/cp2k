"""Direct Python access to CP2K, without starting a CP2K subprocess.

The native API uses atomic units (bohr, hartree). Importing this package does
not load libcp2k or initialize MPI. See :class:`CP2K` for runtime ownership.
"""

# SPDX-License-Identifier: GPL-2.0-or-later

from .input import input_to_string
from .library import CP2K, CalculationResult, ForceEnvironment, SCFConvergenceError

__all__ = [
    "CP2K",
    "CalculationResult",
    "ForceEnvironment",
    "SCFConvergenceError",
    "input_to_string",
]
__version__ = "0.1.0"
