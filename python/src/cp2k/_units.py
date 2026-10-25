"""CODATA-2006 conversions maintained alongside src/common/physcon.F.

Keep the table values and derived expressions synchronized with physcon.
Source: https://doi.org/10.1103/RevModPhys.80.633
"""

# SPDX-License-Identifier: GPL-2.0-or-later

CODATA_YEAR = 2006

# Exact SI value for this CODATA adjustment, matching physcon::c_light [m/s].
_C_LIGHT = 299792458.0

# Recommended table entries, with physcon names and SI units.
_H_PLANCK = 6.62606896e-34  # h_planck [J s]
_E_CHARGE = 1.602176487e-19  # e_charge [C]
_RYDBERG = 10973731.568527  # rydberg [1/m]
_N_AVOGADRO = 6.02214179e23  # n_avogadro [1/mol]
_A_BOHR = 0.52917720859e-10  # a_bohr [m]

# Derived conversions, using the same expressions as physcon.
BOHR_TO_ANGSTROM = 1.0e10 * _A_BOHR
BOHR_TO_NM = BOHR_TO_ANGSTROM * 0.1
HARTREE_TO_JOULE = 2.0 * _RYDBERG * _H_PLANCK * _C_LIGHT
HARTREE_TO_EV = HARTREE_TO_JOULE / _E_CHARGE
HARTREE_TO_KJMOL = 0.001 * HARTREE_TO_JOULE * _N_AVOGADRO
HARTREE_TO_KCALMOL = HARTREE_TO_KJMOL / 4.184
