# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.
"""
Tests for the ASE Calculator

.. important::

    The ASE calculator is taking input coordinates in Angstrom and returns
    energies in eV (while the library is working internally in atomic units).
"""

import numpy as np
import pytest
from pytest import approx

try:
    import ase
    from ase.atoms import Atoms
    from tblite.ase import TBLite
except ModuleNotFoundError:
    ase = None


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_gfn2_xtb_0d():
    """Test ASE interface to GFN2-xTB"""
    thr = 1.0e-5

    atoms = Atoms(
        symbols="CHOCH2CH4CH2OH",
        positions=np.array(
            [
                [+1.578385, +0.147690, +0.343809],
                [+1.394750, +0.012968, +1.413545],
                [+1.359929, -1.086203, -0.359782],
                [+0.653845, +1.215099, -0.221322],
                [+1.057827, +2.180283, +0.093924],
                [+0.729693, +1.184864, -1.311438],
                [-0.817334, +1.152127, +0.208156],
                [-1.303525, +2.065738, -0.145828],
                [-0.883765, +1.159762, +1.299260],
                [+1.984120, -1.734446, -0.021385],
                [+2.616286, +0.458948, +0.206544],
                [-1.627725, -0.034052, -0.311301],
                [-2.684229, +0.151015, -0.118566],
                [-1.501868, -0.118146, -1.397506],
                [-1.324262, -1.260154, +0.333377],
                [-0.417651, -1.475314, +0.076637],
            ]
        ),
    )
    forces = np.array(
        [
            [-0.28558053, -0.48587202, -0.28390121],
            [-0.05526630, -0.07403567, +0.09666425],
            [+0.18832516, +0.39158805, +0.33883777],
            [-0.01166544, +0.24654770, -0.16778463],
            [-0.05367098, +0.03368315, +0.05115390],
            [-0.09604839, -0.10000234, +0.01095881],
            [+0.09420778, +0.18571230, +0.17795940],
            [+0.02650910, -0.04070420, -0.03657800],
            [+0.07887078, -0.05806412, +0.00652814],
            [-0.00787394, -0.07974253, -0.03180615],
            [+0.13328525, +0.03209095, -0.04638655],
            [+0.08200627, -0.39161140, +0.12110936],
            [-0.11455254, -0.01484292, +0.09974688],
            [+0.09786308, -0.09130266, -0.05742520],
            [-0.26641802, +0.47603925, -0.27857176],
            [+0.19000871, -0.02948355, -0.00050500],
        ]
    )

    TBLite(method="GFN2-xTB", atoms=atoms)

    assert approx(atoms.get_potential_energy(), abs=thr) == -592.6794366990786
    assert approx(atoms.get_forces(), abs=thr) == forces


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_gfn1_xtb_0d():
    """Test ASE interface to GFN1-xTB"""
    thr = 1.0e-5

    atoms = Atoms(
        symbols="CHOCH2CH4CH2OH",
        positions=np.array(
            [
                [+1.578385, +0.147690, +0.343809],
                [+1.394750, +0.012968, +1.413545],
                [+1.359929, -1.086203, -0.359782],
                [+0.653845, +1.215099, -0.221322],
                [+1.057827, +2.180283, +0.093924],
                [+0.729693, +1.184864, -1.311438],
                [-0.817334, +1.152127, +0.208156],
                [-1.303525, +2.065738, -0.145828],
                [-0.883765, +1.159762, +1.299260],
                [+1.984120, -1.734446, -0.021385],
                [+2.616286, +0.458948, +0.206544],
                [-1.627725, -0.034052, -0.311301],
                [-2.684229, +0.151015, -0.118566],
                [-1.501868, -0.118146, -1.397506],
                [-1.324262, -1.260154, +0.333377],
                [-0.417651, -1.475314, +0.076637],
            ]
        ),
    )
    forces = np.array(
        [
            [-0.37070590, -0.51067739, -0.27981764],
            [-0.04339461, -0.09290876, +0.22940156],
            [+0.11141234, +0.46678720, +0.24552625],
            [+0.04255709, +0.19019316, -0.23531997],
            [-0.01897377, +0.10810803, +0.05314982],
            [-0.07150720, -0.05182148, -0.08413638],
            [+0.06631826, +0.10587709, +0.29833479],
            [-0.01062355, +0.02301460, -0.04964730],
            [+0.06610108, -0.02724994, +0.09234280],
            [+0.06519070, -0.19311773, -0.01152205],
            [+0.23879786, +0.09871398, -0.04009526],
            [-0.04381577, -0.49997745, +0.08672818],
            [-0.23259608, +0.13735636, +0.06783414],
            [+0.08297636, -0.09566973, -0.20602954],
            [-0.23686052, +0.57454371, -0.17194215],
            [+0.35512370, -0.23317164, +0.00519275],
        ]
    )

    atoms.calc = TBLite(method="GFN1-xTB")

    assert approx(atoms.get_potential_energy(), abs=thr) == -632.7363734598027
    assert approx(atoms.get_forces(), abs=thr) == forces


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_gfn1_xtb_3d():
    """Test ASE interface to GFN1-xTB"""
    thr = 5.0e-6

    atoms = Atoms(
        symbols="C4O8",
        positions=np.array(
            [
                [0.9441259872, 0.9437851680, 0.9543505632],
                [3.7179966528, 0.9556570368, 3.7316862240],
                [3.7159517376, 3.7149292800, 0.9692330016],
                [0.9529872864, 3.7220864832, 3.7296981120],
                [1.6213905408, 1.6190616096, 1.6313879040],
                [0.2656685664, 0.2694175776, 0.2776540416],
                [4.3914553920, 1.6346256864, 3.0545920000],
                [3.0440834880, 0.2764611744, 4.4080419264],
                [4.3910577696, 3.0416409504, 0.2881058304],
                [3.0399936576, 4.3879335936, 1.6497353376],
                [0.2741322432, 4.4003734944, 3.0573754368],
                [1.6312174944, 3.0434586528, 4.4023048032],
            ]
        ),
        cell=np.array([5.68032, 5.68032, 5.68032]),
        pbc=np.array([True, True, True]),
    )
    forces = np.array(
        [
            [-0.08831857, -0.07001440, -0.07468802],
            [-0.03556650, -0.02242172, +0.03047670],
            [+0.03228741, -0.03948068, -0.02892496],
            [-0.02568948, +0.03372976, -0.03161895],
            [-1.90307589, -1.90237488, -1.90613228],
            [+1.98862186, +1.96959087, +1.97850376],
            [-1.88899728, -1.93509758, +1.91693054],
            [+1.92988790, +1.95062169, -1.94117459],
            [-1.93845751, +1.93069936, +1.96027450],
            [+1.91146979, -1.88621511, -1.93902798],
            [+1.94937085, -1.94761107, +1.92150623],
            [-1.93153260, +1.91857377, -1.88612496],
        ]
    )
    stress = np.array(
        [
            +4.49045792e-02,
            +4.49168887e-02,
            +4.49566951e-02,
            +3.38245641e-05,
            +1.52117499e-05,
            +1.13328271e-04,
        ]
    )

    atoms.calc = TBLite(method="GFN1-xTB")
    assert atoms.pbc.all()

    assert approx(atoms.get_potential_energy(), abs=thr) == -1257.0801067985549
    assert approx(atoms.get_forces(), abs=thr) == forces
    assert approx(atoms.get_stress(), abs=thr) == stress


def get_crcp2():
    """Get structure for CrCP2"""

    atoms = Atoms(
        symbols="CrC5H5C3HCHCH3",
        positions=ase.units.Bohr
        * np.array(
            [
                [+0.00000000000000, +0.00000000000000, -0.06044684528305],
                [+0.00000000000000, +3.19613712523833, +2.30877824528580],
                [+2.18828801115897, +3.32943780995850, +0.70249948585735],
                [+1.33235791539260, +3.55640652898451, -1.83908673090077],
                [-1.33235791539260, +3.55640652898451, -1.83908673090077],
                [-2.18828801115897, +3.32943780995850, +0.70249948585735],
                [+0.00000000000000, +3.10509505378016, +4.34935395653655],
                [+4.13810718850644, +3.28428734944129, +1.31235006648465],
                [+2.52190264478215, +3.60569548880831, -3.50208900904436],
                [-2.52190264478215, +3.60569548880831, -3.50208900904436],
                [-4.13810718850644, +3.28428734944129, +1.31235006648465],
                [+2.18828801115897, -3.32943780995850, +0.70249948585735],
                [+0.00000000000000, -3.19613712523833, +2.30877824528580],
                [+1.33235791539260, -3.55640652898451, -1.83908673090077],
                [+4.13810718850644, -3.28428734944129, +1.31235006648465],
                [-2.18828801115897, -3.32943780995850, +0.70249948585735],
                [+0.00000000000000, -3.10509505378016, +4.34935395653655],
                [-1.33235791539260, -3.55640652898451, -1.83908673090077],
                [+2.52190264478215, -3.60569548880831, -3.50208900904436],
                [-4.13810718850644, -3.28428734944129, +1.31235006648465],
                [-2.52190264478215, -3.60569548880831, -3.50208900904436],
            ]
        ),
    )

    return atoms


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_spgfn1_xtb():
    """Test ASE interface to spGFN1-xTB"""
    thr = 5.0e-6

    atoms = get_crcp2()

    atoms.calc = TBLite(method="GFN1-xTB")
    assert approx(atoms.get_potential_energy(), abs=thr) == -771.4322856679416

    atoms.calc.set(spin_polarization=1.0)
    assert approx(atoms.get_potential_energy(), abs=thr) == -771.4322856679416

    atoms.calc.set(multiplicity=3)
    assert approx(atoms.get_potential_energy(), abs=thr) == -772.1635105495686


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_solvation_gfn2_xtb_cpcm():
    """Test CPCM solvation with GFN2-xTB"""
    thr = 5.0e-5 # currently loose testing due to non-variational CPCM

    atoms = get_crcp2()

    atoms.calc = TBLite(method="GFN2-xTB")
    atoms.calc.set(accuracy=0.1)

    atoms.calc.set(cpcm_solvation=7.0)
    assert approx(atoms.get_potential_energy(), abs=thr) == -773.6978494954839
                                                            

@pytest.mark.skipif(ase is None, reason="requires ase")
def test_solvation_gfn2_xtb_alpb():
    """Test ALPB solvation with GFN2-xTB"""
    thr = 5.0e-6

    atoms = get_crcp2()

    atoms.calc = TBLite(method="GFN2-xTB")
    atoms.calc.set(accuracy=1.0)

    atoms.calc.set(alpb_solvation="ethanol")
    assert approx(atoms.get_potential_energy(), abs=thr) == -774.1242966319087

    atoms.calc.set(alpb_solvation=("ethanol", "bar1mol"))
    assert approx(atoms.get_potential_energy(), abs=thr) == -774.0418125853236

    atoms.calc.set(alpb_solvation=("ethanol", "reference"))
    assert approx(atoms.get_potential_energy(), abs=thr) == -773.9688203669275


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_solvation_gfn1_xtb_alpb():
    """Test ALPB solvation with GFN1-xTB"""
    thr = 5.0e-6

    atoms = get_crcp2()

    atoms.calc = TBLite(method="GFN1-xTB")
    atoms.calc.set(accuracy=1.0)

    atoms.calc.set(alpb_solvation="dmf")
    print(atoms.get_potential_energy())
    assert approx(atoms.get_potential_energy(), abs=thr) == -771.7287431921513

    atoms.calc.set(alpb_solvation=("dmf", "bar1mol"))
    print(atoms.get_potential_energy())
    assert approx(atoms.get_potential_energy(), abs=thr) == -771.6462591455721

    atoms.calc.set(alpb_solvation=("dmf", "reference"))
    print(atoms.get_potential_energy())
    assert approx(atoms.get_potential_energy(), abs=thr) == -771.5803670939658


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_solvation_gfn1_xtb_gbe():
    """Test GBE solvation with GFN1-xTB"""
    thr = 5.0e-6

    atoms = get_crcp2()

    atoms.calc = TBLite(method="GFN1-xTB")
    atoms.calc.set(accuracy=1.0)

    atoms.calc.set(gbe_solvation=(7.0, "p16"))
    assert approx(atoms.get_potential_energy(), abs=thr) == -771.4434811395378


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_solvation_gfn2_xtb_gbsa():
    """Test GBSA solvation with GFN2-xTB"""
    thr = 5.0e-6

    atoms = get_crcp2()

    atoms.calc = TBLite(method="GFN2-xTB")
    atoms.calc.set(accuracy=1.0)

    atoms.calc.set(gbsa_solvation="water")
    assert approx(atoms.get_potential_energy(), abs=thr) == -773.8895533357601

    atoms.calc.set(gbsa_solvation=("water", "gsolv"))
    assert approx(atoms.get_potential_energy(), abs=thr) == -773.8895533357601

    atoms.calc.set(gbsa_solvation=("water", "bar1mol"))
    assert approx(atoms.get_potential_energy(), abs=thr) == -773.8070696428758


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_solvation_gfn2_xtb_gb():
    """Test GB solvation with GFN2-xTB"""
    thr = 5.0e-6

    atoms = get_crcp2()

    atoms.calc = TBLite(method="GFN2-xTB")
    atoms.calc.set(accuracy=1.0)

    atoms.calc.set(gb_solvation=(7.0, "still"))
    assert approx(atoms.get_potential_energy(), abs=thr) == -773.8038793721613
