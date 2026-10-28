"""
CODATA-2022 constants and conversions, maintained alongside src/common/physcon.F.

Source: https://doi.org/10.1103/RevModPhys.97.025002
"""

# Exact SI defining constants.
_C_LIGHT = 299792458.0
_H_PLANCK = 6.62607015e-34
_E_CHARGE = 1.602176634e-19
_BOLTZMANN = 1.380649e-23

# Recommended CODATA table entries.
_E_MASS = 9.1093837139e-31
_A_MASS = 1.66053906892e-27
_RYDBERG = 10973731.568157
_A_BOHR = 5.29177210544e-11

# Derived conversions, using the same expressions as physcon.
E_MASS = _E_MASS / _A_MASS  # Electron mass in u.
AU_TO_JOULES = 2.0 * _RYDBERG * _H_PLANCK * _C_LIGHT
AU_TO_EV = AU_TO_JOULES / _E_CHARGE
EV_TO_AU = 1 / AU_TO_EV
WAVENUMBER_TO_AU = 1 / (0.02 * _RYDBERG)
kB = _BOLTZMANN / AU_TO_JOULES
BOHR_TO_ANGSTROM = 1.0e10 * _A_BOHR

# Atomic masses in amu (Dalton)
# Source: www.angelo.edu/faculty/kboudrea/periodic/structure_mass.htm
ATOMIC_MASSES = {
    "H": 1.00797,
    "He": 4.00260,
    "Li": 6.941,
    "Be": 9.01218,
    "B": 10.81,
    "C": 12.011,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.998403,
    "Ne": 20.179,
    "Na": 22.98977,
    "Mg": 24.305,
    "Al": 26.98154,
    "Si": 28.0855,
    "P": 30.97376,
    "S": 32.06,
    "Cl": 35.453,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.08,
    "Sc": 44.9559,
    "Ti": 47.90,
    "V": 50.9415,
    "Cr": 51.996,
    "Mn": 54.9380,
    "Fe": 55.847,
    "Co": 58.9332,
    "Ni": 58.70,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.72,
    "Ge": 72.59,
    "As": 74.9216,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.80,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.9059,
    "Zr": 91.22,
    "Nb": 92.9064,
    "Mo": 95.94,
    "Ru": 101.07,
    "Rh": 102.9055,
    "Pd": 106.4,
    "Ag": 107.868,
    "Cd": 112.41,
    "In": 114.82,
    "Sn": 118.69,
    "Sb": 121.75,
    "I": 126.9045,
    "Te": 127.60,
    "Xe": 131.30,
    "Cs": 132.9054,
    "Ba": 137.33,
    "La": 138.9055,
    "Ce": 140.12,
    "Pr": 140.9077,
    "Nd": 144.24,
    "Sm": 150.4,
    "Eu": 151.96,
    "Gd": 157.25,
    "Tb": 158.9254,
    "Dy": 162.50,
    "Ho": 164.9304,
    "Er": 167.26,
    "Tm": 168.9342,
    "Yb": 173.04,
    "Lu": 174.967,
    "Hf": 178.49,
    "Ta": 180.9479,
    "W": 183.85,
    "Re": 186.207,
    "Os": 190.2,
    "Ir": 192.22,
    "Pt": 195.09,
    "Au": 196.9665,
    "Hg": 200.59,
    "Tl": 204.37,
    "Pb": 207.2,
    "Bi": 208.9804,
}
