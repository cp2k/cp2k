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

# pyright: reportPossiblyUnboundVariable=false
# pylint: disable=possibly-used-before-assignment

import numpy as np
from pytest import approx, mark

try:
    from ase.build import molecule
    from ase.calculators.emt import EMT
    from dftd4.ase import DFTD4

    has_ase = True
except ModuleNotFoundError:
    has_ase = False

pytestmark = mark.skipif(not has_ase, reason="requires ase")


def test_ase_scand4() -> None:
    thr = 1.0e-6

    forces = np.array(
        [
            [-7.90552684e-21, -1.15811595e-19, -2.80061133e-05],
            [+4.61502216e-20, +4.26735028e-04, +4.94269127e-04],
            [+2.28682834e-19, -4.26735028e-04, +4.94269127e-04],
            [+1.39725405e-20, -9.71924142e-20, -7.95205193e-04],
            [-1.08042249e-04, -8.23929519e-05, +2.31098749e-05],
            [+1.08042249e-04, -8.23929519e-05, +2.31098749e-05],
            [+1.08042249e-04, +8.23929519e-05, +2.31098749e-05],
            [-1.08042249e-04, +8.23929519e-05, +2.31098749e-05],
            [+1.07391220e-20, -4.98420762e-05, -1.28883224e-04],
            [+5.97028977e-21, +4.98420762e-05, -1.28883224e-04],
        ]
    )

    atoms = molecule("methylenecyclopropane")
    atoms.calc = DFTD4(method="SCAN")

    assert (
        approx(atoms.get_potential_energy(), abs=thr) == -0.021665446836610567
    )
    assert approx(atoms.get_forces(), abs=thr) == forces

    atoms.calc = DFTD4(method="SCAN").add_calculator(EMT())
    assert approx(atoms.get_potential_energy(), abs=thr) == 3.6624398683434225

    if hasattr(atoms.calc, "calcs"):
        calcs = atoms.calc.calcs  # type: ignore[attr-defined]
    else:
        calcs = atoms.calc.mixer.calcs  # type: ignore[attr-defined]

    energies = [calc.get_potential_energy() for calc in calcs]
    assert approx(energies, abs=thr) == [
        -0.021665446836610563,
        3.684105315180033,
    ]


def test_ase_pbed4s() -> None:
    thr = 1.0e-6

    forces = np.array(
        [
            [-1.05900228e-18, -1.47186433e-03, -2.08505399e-03],
            [-8.94454333e-19, 1.47186433e-03, -2.08505399e-03],
            [-2.48990627e-03, -1.06729670e-18, 2.15375860e-03],
            [2.48990627e-03, -1.91047493e-19, 2.15375860e-03],
            [-4.33024525e-22, -2.06112648e-03, -6.95481292e-05],
            [2.47579260e-21, -2.26182709e-03, -5.56758674e-04],
            [3.90654502e-20, 2.06112648e-03, -6.95481292e-05],
            [4.57489295e-20, 2.26182709e-03, -5.56758674e-04],
            [-1.57814347e-03, 6.44754156e-19, 5.57602199e-04],
            [1.57814347e-03, 5.62735538e-19, 5.57602199e-04],
        ]
    )

    atoms = molecule("bicyclobutane")
    atoms.calc = DFTD4(method="PBE", model="d4s")

    assert approx(atoms.get_potential_energy(), abs=thr) == -0.16377494406788423
    assert approx(atoms.get_forces(), abs=thr) == forces

    atoms.calc = DFTD4(method="PBE", model="d4s").add_calculator(EMT())
    assert approx(atoms.get_potential_energy(), abs=thr) == 3.364290912363593

    if hasattr(atoms.calc, "calcs"):
        calcs = atoms.calc.calcs  # type: ignore[attr-defined]
    else:
        calcs = atoms.calc.mixer.calcs  # type: ignore[attr-defined]

    energies = [calc.get_potential_energy() for calc in calcs]
    assert approx(energies, abs=thr) == [
        -0.16377494406788423,
        3.5280658564314775,
    ]

def test_ase_pbed4s_custom_damping() -> None:
    thr = 1.0e-12
    atoms = molecule("bicyclobutane")

    atoms.calc = DFTD4(method="PBE", model="d4s", damping_hint={"2b": "rational", "3b": "zero-avg"})

    print(atoms.get_potential_energy())

    assert approx(atoms.get_potential_energy(), abs=thr) == -0.16377494406788423


def test_ase_pbed4_noatm() -> None:
    thr = 1.0e-12
    atoms = molecule("bicyclobutane")

    atoms.calc = DFTD4(method="PBE", model="d4", damping_hint={"3b": "none"})
    assert approx(atoms.get_potential_energy(), abs=thr) == -0.15427642282371362


def test_ase_tpssd4() -> None:
    thr = 1.0e-6

    forces = np.array(
        [
            [-8.27697238e-04, -1.74116189e-02, +9.50315402e-05],
            [-3.00615825e-04, +3.51410800e-04, -3.12339518e-03],
            [-9.55691674e-04, -4.05325537e-03, +7.33934018e-05],
            [+2.70525425e-04, -1.91784287e-03, -6.13796425e-04],
            [+1.56444473e-02, +5.71643192e-03, +6.29706049e-04],
            [-1.43399413e-02, +7.36397630e-03, +7.87584027e-04],
            [+4.41551907e-03, +4.04705396e-04, +1.42098826e-03],
            [-4.17670039e-03, +1.57923335e-03, +1.60488604e-03],
            [+4.66065256e-03, +7.97764912e-05, -3.60249901e-05],
            [-4.95386767e-03, +1.38115911e-03, -1.63079386e-05],
            [+5.36717422e-03, +2.78165913e-03, -4.68647341e-04],
            [-4.80380448e-03, +3.72436466e-03, -3.53417437e-04],
        ]
    )

    atoms = molecule("C2H6CHOH")
    atoms.calc = DFTD4(method="TPSS")

    assert approx(atoms.get_potential_energy(), abs=thr) == -0.24206732765720423
    assert approx(atoms.get_forces(), abs=thr) == forces

    atoms.calc = DFTD4(method="TPSS").add_calculator(EMT())
    assert approx(atoms.get_potential_energy(), abs=thr) == 4.864016486351274

    if hasattr(atoms.calc, "calcs"):
        calcs = atoms.calc.calcs  # type: ignore[attr-defined]
    else:
        calcs = atoms.calc.mixer.calcs  # type: ignore[attr-defined]

    energies = [calc.get_potential_energy() for calc in calcs]
    assert approx(energies, abs=thr) == [
        -0.24206732765720396,
        5.106083814008478,
    ]
