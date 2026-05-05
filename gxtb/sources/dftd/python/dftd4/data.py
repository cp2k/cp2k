# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

from math import sqrt

import numpy as np

try:
    from os.path import dirname, join

    from ase.collections import Collection

    # using a collection will remove the data, but we get at least the structures
    references = Collection("references")
    # need to patch the collection immediately
    references.filename = join(dirname(__file__), references.name + ".json")

except ModuleNotFoundError:
    references = None


# covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197),
# values for metals decreased by 10 %
covalent_radii = np.array(
    [
        0.00,  # None
        0.32,  # H
        0.46,  # He
        1.20,  # Li (2nd)
        0.94,  # Be
        0.77,  # B
        0.75,  # C
        0.71,  # N
        0.63,  # O
        0.64,  # F
        0.67,  # Ne
        1.40,  # Na (3rd)
        1.25,  # Mg
        1.13,  # Al
        1.04,  # Si
        1.10,  # P
        1.02,  # S
        0.99,  # Cl
        0.96,  # Ar
        1.76,  # K  (4th)
        1.54,  # Ca
        1.33,  # Sc
        1.22,  # Ti
        1.21,  # V
        1.10,  # Cr
        1.07,  # Mn
        1.04,  # Fe
        1.00,  # Co
        0.99,  # Ni
        1.01,  # Cu
        1.09,  # Zn
        1.12,  # Ga
        1.09,  # Ge
        1.15,  # As
        1.10,  # Se
        1.14,  # Br
        1.17,  # Kr
        1.89,  # Rb (5th)
        1.67,  # Sr
        1.47,  # Y
        1.39,  # Zr
        1.32,  # Nb
        1.24,  # Mo
        1.15,  # Tc
        1.13,  # Ru
        1.13,  # Rh
        1.08,  # Pd
        1.15,  # Ag
        1.23,  # Cd
        1.28,  # In
        1.26,  # Sn
        1.26,  # Sb
        1.23,  # Te
        1.32,  # I
        1.31,  # Xe
        2.09,  # Cs (6th)
        1.76,  # Ba
        1.62,  # La
        1.47,  # Ce
        1.58,  # Pr
        1.57,  # Nd
        1.56,  # Pm
        1.55,  # Sm
        1.51,  # Eu
        1.52,  # Gd
        1.51,  # Tb
        1.50,  # Dy
        1.49,  # Ho
        1.49,  # Er
        1.48,  # Tm
        1.53,  # Yb
        1.46,  # Lu
        1.37,  # Hf
        1.31,  # Ta
        1.23,  # W
        1.18,  # Re
        1.16,  # Os
        1.11,  # Ir
        1.12,  # Pt
        1.13,  # Au
        1.32,  # Hg
        1.30,  # Tl
        1.30,  # Pb
        1.36,  # Bi
        1.31,  # Po
        1.38,  # At
        1.42,  # Rn
        2.01,  # Fr (7th)
        1.81,  # Ra
        1.67,  # Ac
        1.58,  # Th
        1.52,  # Pa
        1.53,  # U
        1.54,  # Np
        1.55,  # Pu
        1.49,  # Am
        1.49,  # Cm
        1.51,  # Bk
        1.51,  # Cf
        1.48,  # Es
        1.50,  # Fm
        1.56,  # Md
        1.58,  # No
        1.45,  # Lr
        1.41,  # Rf
        1.34,  # Db
        1.29,  # Sg
        1.27,  # Bh
        1.21,  # Hs
        1.16,  # Mt
        1.15,  # Ds
        1.09,  # Rg
        1.22,  # Cn
        1.36,  # Nh
        1.43,  # Fl
        1.46,  # Lv
        1.58,  # Mc
        1.48,  # Ts
        1.57,  # Og
    ]
)

pauling_en = np.array(
    [
        1.5,  # None
        2.20,  # H
        3.00,  # He
        0.98,  # Li (2nd)
        1.57,  # Be
        2.04,  # B
        2.55,  # C
        3.04,  # N
        3.44,  # O
        3.98,  # F
        4.50,  # Ne
        0.93,  # Na (3rd)
        1.31,  # Mg
        1.61,  # Al
        1.90,  # Si
        2.19,  # P
        2.58,  # S
        3.16,  # Cl
        3.50,  # Ar
        0.82,  # K  (4th)
        1.00,  # Ca
        1.36,  # Sc
        1.54,  # Ti
        1.63,  # V
        1.66,  # Cr
        1.55,  # Mn
        1.83,  # Fe
        1.88,  # Co
        1.91,  # Ni
        1.90,  # Cu
        1.65,  # Zn
        1.81,  # Ga
        2.01,  # Ge
        2.18,  # As
        2.55,  # Se
        2.96,  # Br
        3.00,  # Kr
        0.82,  # Rb (5th)
        0.95,  # Sr
        1.22,  # Y
        1.33,  # Zr
        1.60,  # Nb
        2.16,  # Mo
        1.90,  # Tc
        2.20,  # Ru
        2.28,  # Rh
        2.20,  # Pd
        1.93,  # Ag
        1.69,  # Cd
        1.78,  # In
        1.96,  # Sn
        2.05,  # Sb
        2.10,  # Te
        2.66,  # I
        2.60,  # Xe
        0.79,  # Cs (6th)
        0.89,  # Ba
        1.10,  # La
        1.12,  # Ce
        1.13,  # Pr
        1.14,  # Nd
        1.15,  # Pm
        1.17,  # Sm
        1.18,  # Eu
        1.20,  # Gd
        1.21,  # Tb
        1.22,  # Dy
        1.23,  # Ho
        1.24,  # Er
        1.25,  # Tm
        1.26,  # Yb
        1.27,  # Lu
        1.30,  # Hf
        1.50,  # Ta
        2.36,  # W
        1.90,  # Re
        2.20,  # Os
        2.20,  # Ir
        2.28,  # Pt
        2.54,  # Au
        2.00,  # Hg
        1.62,  # Tl
        2.33,  # Pb
        2.02,  # Bi
        2.00,  # Po
        2.20,  # At
        2.20,  # Rn
        0.79,  # Fr (7th)
        0.90,  # Ra
        1.10,  # Ac
        1.30,  # Th
        1.50,  # Pa
        1.38,  # U
        1.36,  # Np
        1.28,  # Pu
        1.13,  # Am
        1.28,  # Cm
        1.30,  # Bk
        1.30,  # Cf
        1.30,  # Es
        1.30,  # Fm
        1.30,  # Md
        1.30,  # No
        1.30,  # Lr
        1.50,  # Rf (only dummies from here)
        1.50,  # Db
        1.50,  # Sg
        1.50,  # Bh
        1.50,  # Hs
        1.50,  # Mt
        1.50,  # Ds
        1.50,  # Rg
        1.50,  # Cn
        1.50,  # Nh
        1.50,  # Fl
        1.50,  # Lv
        1.50,  # Mc
        1.50,  # Ts
        1.50,  # Og
    ]
)

# Semiempirical Evaluation of the GlobalHardness of the Atoms of 103 Elements of
# the Periodic Table Using the Most Probable Radii as their Size Descriptors
# DULAL C. GHOSH, NAZMUL ISLAM 2009 in Wiley InterScience. DOI 10.1002/qua.22202
#
# values in the paper multiplied by two because (ii:ii)=(IP-EA)=d^2 E/dN^2
# but the hardness definition they use is 1/2d^2 E/dN^2 (in Eh)
chemical_hardness = np.array(
    [
        0.0,  # None
        0.47259288,  # H
        0.92203391,  # He
        0.17452888,  # Li (2nd)
        0.25700733,  # Be
        0.33949086,  # B
        0.42195412,  # C
        0.50438193,  # N
        0.58691863,  # O
        0.66931351,  # F
        0.75191607,  # Ne
        0.17964105,  # Na (3rd)
        0.22157276,  # Mg
        0.26348578,  # Al
        0.30539645,  # Si
        0.34734014,  # P
        0.38924725,  # S
        0.43115670,  # Cl
        0.47308269,  # Ar
        0.17105469,  # K  (4th)
        0.20276244,  # Ca
        0.21007322,  # Sc
        0.21739647,  # Ti
        0.22471039,  # V
        0.23201501,  # Cr
        0.23933969,  # Mn
        0.24665638,  # Fe
        0.25398255,  # Co
        0.26128863,  # Ni
        0.26859476,  # Cu
        0.27592565,  # Zn
        0.30762999,  # Ga
        0.33931580,  # Ge
        0.37235985,  # As
        0.40273549,  # Se
        0.43445776,  # Br
        0.46611708,  # Kr
        0.15585079,  # Rb (5th)
        0.18649324,  # Sr
        0.19356210,  # Y
        0.20063311,  # Zr
        0.20770522,  # Nb
        0.21477254,  # Mo
        0.22184614,  # Tc
        0.22891872,  # Ru
        0.23598621,  # Rh
        0.24305612,  # Pd
        0.25013018,  # Ag
        0.25719937,  # Cd
        0.28784780,  # In
        0.31848673,  # Sn
        0.34912431,  # Sb
        0.37976593,  # Te
        0.41040808,  # I
        0.44105777,  # Xe
        0.05019332,  # Cs (6th)
        0.06762570,  # Ba
        0.08504445,  # La
        0.10247736,  # Ce
        0.11991105,  # Pr
        0.13732772,  # Nd
        0.15476297,  # Pm
        0.17218265,  # Sm
        0.18961288,  # Eu
        0.20704760,  # Gd
        0.22446752,  # Tb
        0.24189645,  # Dy
        0.25932503,  # Ho
        0.27676094,  # Er
        0.29418231,  # Tm
        0.31159587,  # Yb
        0.32902274,  # Lu
        0.34592298,  # Hf
        0.36388048,  # Ta
        0.38130586,  # W
        0.39877476,  # Re
        0.41614298,  # Os
        0.43364510,  # Ir
        0.45104014,  # Pt
        0.46848986,  # Au
        0.48584550,  # Hg
        0.12526730,  # Tl
        0.14268677,  # Pb
        0.16011615,  # Bi
        0.17755889,  # Po
        0.19497557,  # At
        0.21240778,  # Rn
        0.07263525,  # Fr (7th)
        0.09422158,  # Ra
        0.09920295,  # Ac
        0.10418621,  # Th
        0.14235633,  # Pa
        0.16394294,  # U
        0.18551941,  # Np
        0.22370139,  # Pu
        0.25110000,  # Am
        0.25030000,  # Cm
        0.28840000,  # Bk
        0.31000000,  # Cf
        0.33160000,  # Es
        0.35320000,  # Fm
        0.36820000,  # Md
        0.39630000,  # No
        0.40140000,  # Lr
        0.00000000,  # Rf
        0.00000000,  # Db
        0.00000000,  # Sg
        0.00000000,  # Bh
        0.00000000,  # Hs
        0.00000000,  # Mt
        0.00000000,  # Ds
        0.00000000,  # Rg
        0.00000000,  # Cn
        0.00000000,  # Nh
        0.00000000,  # Fl
        0.00000000,  # Lv
        0.00000000,  # Mc
        0.00000000,  # Ts
        0.00000000,  # Og
    ]
)

def2ecp_nuclear_charges = np.array(
    [
        0,  # None
        1,  # H
        2,  # He
        3,  # Li (2nd)
        4,  # Be
        5,  # B
        6,  # C
        7,  # N
        8,  # O
        9,  # F
        10,  # Ne
        11,  # Na (3rd)
        12,  # Mg
        13,  # Al
        14,  # Si
        15,  # P
        16,  # S
        17,  # Cl
        18,  # Ar
        19,  # K  (4th)
        20,  # Ca
        21,  # Sc
        22,  # Ti
        23,  # V
        24,  # Cr
        25,  # Mn
        26,  # Fe
        27,  # Co
        28,  # Ni
        29,  # Cu
        30,  # Zn
        31,  # Ga
        32,  # Ge
        33,  # As
        34,  # Se
        35,  # Br
        36,  # Kr
        9,  # Rb (5th)
        10,  # Sr
        11,  # Y
        12,  # Zr
        13,  # Nb
        14,  # Mo
        15,  # Tc
        16,  # Ru
        17,  # Rh
        18,  # Pd
        19,  # Ag
        20,  # Cd
        21,  # In
        22,  # Sn
        23,  # Sb
        24,  # Te
        25,  # I
        26,  # Xe
        9,  # Cs (6th)
        10,  # Ba
        11,  # La
        30,  # Ce
        31,  # Pr
        32,  # Nd
        33,  # Pm
        34,  # Sm
        35,  # Eu
        36,  # Gd
        37,  # Tb
        38,  # Dy
        39,  # Ho
        40,  # Er
        41,  # Tm
        42,  # Yb
        43,  # Lu
        12,  # Hf
        13,  # Ta
        14,  # W
        15,  # Re
        16,  # Os
        17,  # Ir
        18,  # Pt
        19,  # Au
        20,  # Hg
        21,  # Tl
        22,  # Pb
        23,  # Bi
        24,  # Po
        25,  # At
        26,  # Rn
        9,  # Fr (7th)
        10,  # Ra
        11,  # Ac
        30,  # Th
        31,  # Pa
        32,  # U
        33,  # Np
        34,  # Pu
        35,  # Am
        36,  # Cm
        37,  # Bk
        38,  # Cf
        39,  # Es
        40,  # Fm
        41,  # Md
        42,  # No
        43,  # Lr
        12,  # Rf
        13,  # Db
        14,  # Sg
        15,  # Bh
        16,  # Hs
        17,  # Mt
        18,  # Ds
        19,  # Rg
        20,  # Cn
        21,  # Nh
        22,  # Fl
        23,  # Lv
        24,  # Mc
        25,  # Ts
        26,  # Og
    ]
)

# PBE0/def2-QZVP atomic values calculated by S. Grimme in Gaussian (2010)
# rare gases recalculated by J. Mewes with PBE0/aug-cc-pVQZ in Dirac (2018)
# He: 3.4698 -> 3.5544, Ne: 3.1036 -> 3.7943, Ar: 5.6004 -> 5.6638,
# Kr: 6.1971 -> 6.2312, Xe: 7.5152 -> 8.8367
# also new super heavies Cn--Og
r4_over_r2 = np.array(
    [
        0.0,  # None
        8.0589,  # H
        3.4698,  # He
        29.0974,  # Li (2nd)
        14.8517,  # Be
        11.8799,  # B
        7.8715,  # C
        5.5588,  # N
        4.7566,  # O
        3.8025,  # F
        3.1036,  # Ne
        26.1552,  # Na (3rd)
        17.2304,  # Mg
        17.7210,  # Al
        12.7442,  # Si
        9.5361,  # P
        8.1652,  # S
        6.7463,  # Cl
        5.6004,  # Ar
        29.2012,  # K  (4th)
        22.3934,  # Ca
        19.0598,  # Sc
        16.8590,  # Ti
        15.4023,  # V
        12.5589,  # Cr
        13.4788,  # Mn
        12.2309,  # Fe
        11.2809,  # Co
        10.5569,  # Ni
        10.1428,  # Cu
        9.4907,  # Zn
        13.4606,  # Ga
        10.8544,  # Ge
        8.9386,  # As
        8.1350,  # Se
        7.1251,  # Br
        6.1971,  # Kr
        30.0162,  # Rb (5th)
        24.4103,  # Sr
        20.3537,  # Y
        17.4780,  # Zr
        13.5528,  # Nb
        11.8451,  # Mo
        11.0355,  # Tc
        10.1997,  # Ru
        9.5414,  # Rh
        9.0061,  # Pd
        8.6417,  # Ag
        8.9975,  # Cd
        14.0834,  # In
        11.8333,  # Sn
        10.0179,  # Sb
        9.3844,  # Te
        8.4110,  # I
        7.5152,  # Xe
        32.7622,  # Cs (6th)
        27.5708,  # Ba
        23.1671,  # La
        21.6003,  # Ce
        20.9615,  # Pr
        20.4562,  # Nd
        20.1010,  # Pm
        19.7475,  # Sm
        19.4828,  # Eu
        15.6013,  # Gd
        19.2362,  # Tb
        17.4717,  # Dy
        17.8321,  # Ho
        17.4237,  # Er
        17.1954,  # Tm
        17.1631,  # Yb
        14.5716,  # Lu
        15.8758,  # Hf
        13.8989,  # Ta
        12.4834,  # W
        11.4421,  # Re
        10.2671,  # Os
        8.3549,  # Ir
        7.8496,  # Pt
        7.3278,  # Au
        7.4820,  # Hg
        13.5124,  # Tl
        11.6554,  # Pb
        10.0959,  # Bi
        9.7340,  # Po
        8.8584,  # At
        8.0125,  # Rn
        29.8135,  # Fr (7th)
        26.3157,  # Ra
        19.1885,  # Ac
        15.8542,  # Th
        16.1305,  # Pa
        15.6161,  # U
        15.1226,  # Np
        16.1576,  # Pu
        14.6510,  # Am
        14.7178,  # Cm
        13.9108,  # Bk
        13.5623,  # Cf
        13.2326,  # Es
        12.9189,  # Fm
        12.6133,  # Md
        12.3142,  # No
        14.8326,  # Lr
        12.3771,  # Rf
        10.6378,  # Db
        9.3638,  # Sg
        8.2297,  # Bh
        7.5667,  # Hs
        6.9456,  # Mt
        6.3946,  # Ds
        5.9159,  # Rg
        5.4929,  # Cn
        6.7286,  # Nh
        6.5144,  # Fl
        10.9169,  # Lv
        10.3600,  # Mc
        9.4723,  # Ts
        8.6641,  # Og
    ]
)

sqrt_z_r4_over_r2 = np.sqrt(
    np.array([0.5 * sqrt(z) for z in range(0, 119)]) * r4_over_r2
)

covalent_radii_d3 = 4.0 / 3.0 * covalent_radii
