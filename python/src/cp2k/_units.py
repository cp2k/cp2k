"""CODATA-2022 conversions maintained alongside src/common/physcon.F.

Keep the table values and derived expressions synchronized with physcon.
Source: https://doi.org/10.1103/RevModPhys.97.025002
"""

# SPDX-License-Identifier: GPL-2.0-or-later

CODATA_YEAR = 2022

# Exact SI defining constants, matching physcon names and SI units.
_C_LIGHT = 299792458.0  # c_light [m/s]
_H_PLANCK = 6.62607015e-34  # h_planck [J s]
_E_CHARGE = 1.602176634e-19  # e_charge [C]
_N_AVOGADRO = 6.02214076e23  # n_avogadro [1/mol]

# Recommended table entries, with physcon names and SI units.
_RYDBERG = 10973731.568157  # rydberg [1/m]
_A_BOHR = 5.29177210544e-11  # a_bohr [m]

# Derived conversions, using the same expressions as physcon.
BOHR_TO_ANGSTROM = 1.0e10 * _A_BOHR
BOHR_TO_NM = BOHR_TO_ANGSTROM * 0.1
HARTREE_TO_JOULE = 2.0 * _RYDBERG * _H_PLANCK * _C_LIGHT
HARTREE_TO_EV = HARTREE_TO_JOULE / _E_CHARGE
HARTREE_TO_KJMOL = 0.001 * HARTREE_TO_JOULE * _N_AVOGADRO
HARTREE_TO_KCALMOL = HARTREE_TO_KJMOL / 4.184
