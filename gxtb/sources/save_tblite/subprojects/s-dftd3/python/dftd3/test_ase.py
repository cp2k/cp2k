# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

from typing import Iterator

import numpy as np
import pytest
from pytest import approx

try:
    import ase
    from dftd3.ase import DFTD3
    from ase.build import molecule
    from ase.calculators.emt import EMT
except ModuleNotFoundError:
    ase = None


def get_calcs(calc) -> Iterator[ase.calculators.calculator.Calculator]:
    if hasattr(calc, "mixer"):
        calc = calc.mixer
    yield from calc.calcs


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_ase_scand3_atm():
    thr = 1.0e-6

    forces = np.array(
        [
            [-0.00000000e-00, -0.00000000e-00, -6.83426991e-05],
            [-0.00000000e-00, +3.44839555e-04, +7.21176947e-04],
            [+6.80565391e-22, -3.44839555e-04, +7.21176947e-04],
            [+2.12676685e-23, +3.26671388e-20, -1.60555514e-03],
            [+3.21599856e-04, +5.54947267e-04, +6.86106874e-04],
            [-3.21599856e-04, +5.54947267e-04, +6.86106874e-04],
            [-3.21599856e-04, -5.54947267e-04, +6.86106874e-04],
            [+3.21599856e-04, -5.54947267e-04, +6.86106874e-04],
            [+1.87155483e-21, +2.87678390e-04, -1.25644177e-03],
            [-3.40282696e-22, -2.87678390e-04, -1.25644177e-03],
        ]
    )

    atoms = molecule("methylenecyclopropane")
    atoms.calc = DFTD3(
        method="SCAN", damping="d3bj", params_tweaks={"method": "SCAN", "atm": True}
    )

    assert atoms.get_potential_energy() == approx(-0.03880921894019244, abs=thr)
    assert atoms.get_forces() == approx(forces, abs=thr)

    atoms.calc = DFTD3(
        method="SCAN", damping="d3bj", params_tweaks={"method": "SCAN", "atm": True}
    ).add_calculator(EMT())
    assert atoms.get_potential_energy() == approx(3.6452960962398406, abs=thr)
    energies = [calc.get_potential_energy() for calc in get_calcs(atoms.calc)]
    assert energies == approx([-0.03880921894019244, 3.684105315180033], abs=thr)


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_ase_scand3():
    thr = 1.0e-6

    forces = np.array(
        [
            [-0.00000000e00, -1.36113078e-21, -6.16745949e-05],
            [-1.08890463e-20, +3.10643080e-04, +6.89816195e-04],
            [+1.08890463e-20, -3.10643080e-04, +6.89816195e-04],
            [-0.00000000e00, -0.00000000e00, -1.54664530e-03],
            [+2.90862635e-04, +5.36974154e-04, +6.87685696e-04],
            [-2.90862635e-04, +5.36974154e-04, +6.87685696e-04],
            [-2.90862635e-04, -5.36974154e-04, +6.87685696e-04],
            [+2.90862635e-04, -5.36974154e-04, +6.87685696e-04],
            [-0.00000000e00, +2.75841655e-04, -1.26102764e-03],
            [-0.00000000e00, -2.75841655e-04, -1.26102764e-03],
        ]
    )

    atoms = molecule("methylenecyclopropane")
    atoms.calc = DFTD3(method="SCAN", damping="d3bj")

    assert atoms.get_potential_energy() == approx(-0.03898577243903914, abs=thr)
    assert atoms.get_forces() == approx(forces, abs=thr)

    atoms.calc = DFTD3(method="SCAN", damping="d3bj").add_calculator(EMT())
    assert atoms.get_potential_energy() == approx(3.645119542740994, abs=thr)
    energies = [calc.get_potential_energy() for calc in get_calcs(atoms.calc)]
    assert energies == approx([-0.03898577243903914, 3.684105315180033], abs=thr)


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_ase_tpssd3_atm():
    thr = 1.0e-6

    forces = np.array(
        [
            [+1.21727790e-03, +1.98579200e-03, -1.16371697e-02],
            [-5.82484114e-04, +9.01770290e-03, +7.78537640e-03],
            [-4.30031958e-03, +4.63213536e-03, -4.56657109e-03],
            [-1.16941383e-03, -8.39071556e-03, +1.60593512e-02],
            [-6.90354443e-03, -5.07801933e-03, -1.75396161e-03],
            [+1.03561818e-02, -1.68908740e-02, -2.74225314e-03],
            [+5.59001294e-03, +3.35129491e-03, -9.24429928e-04],
            [-5.13316989e-03, +6.07626858e-03, +3.89454026e-05],
            [+3.35952011e-03, +3.95424504e-03, -5.65438002e-04],
            [-2.13140242e-03, +2.77295425e-03, -4.76829804e-04],
            [+4.33961724e-03, -1.51731003e-03, -7.01598391e-04],
            [-4.64227572e-03, +8.65258554e-05, -5.15421318e-04],
        ]
    )

    atoms = molecule("C2H6CHOH")
    atoms.calc = DFTD3(
        method="TPSS", damping="d3zero", params_tweaks={"method": "TPSS", "atm": True}
    )

    assert atoms.get_potential_energy() == approx(-0.14230914516094673, abs=thr)
    assert atoms.get_forces() == approx(forces, abs=thr)

    atoms.calc = DFTD3(
        method="TPSS", damping="d3zero", params_tweaks={"method": "TPSS", "atm": True}
    ).add_calculator(EMT())
    assert atoms.get_potential_energy() == approx(4.963774668847532, abs=thr)
    energies = [calc.get_potential_energy() for calc in get_calcs(atoms.calc)]
    assert energies == approx([-0.14230914516094673, 5.106083814008478], abs=thr)


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_ase_realspace_cutoff():
    """Test that realspace_cutoff parameter works and can be updated with cache_api=True"""
    from ase.units import Bohr

    thr = 1.0e-6
    atoms = molecule("H2O")

    # Test with default cutoffs
    calc_default = DFTD3(method="PBE", damping="d3bj", cache_api=True)
    atoms.calc = calc_default
    energy_default = atoms.get_potential_energy()
    forces_default = atoms.get_forces()

    # Test with very small cutoffs (smaller than H2O bond length ~1 Angstrom) - should give zero or very small interactions
    calc_custom = DFTD3(
        method="PBE",
        damping="d3bj",
        realspace_cutoff={"disp2": 0.5, "disp3": 0.5, "cn": 0.5},
        cache_api=True,
    )
    atoms.calc = calc_custom
    energy_custom = atoms.get_potential_energy()
    forces_custom = atoms.get_forces()

    # With very small cutoffs, energy should be much smaller (closer to zero)
    assert abs(energy_custom) < abs(energy_default)
    assert not np.allclose(forces_custom, forces_default, atol=thr)

    # Test updating cutoff via set() with cache_api=True
    # Reset to default-like cutoffs using set()
    calc_custom.set(
        realspace_cutoff={"disp2": 60.0 * Bohr, "disp3": 40.0 * Bohr, "cn": 40.0 * Bohr}
    )
    energy_updated = atoms.get_potential_energy()
    forces_updated = atoms.get_forces()

    # Should match the default energy/forces
    assert energy_updated == approx(energy_default, abs=thr)
    assert forces_updated == approx(forces_default, abs=thr)

    # Test with empty dict (should use library defaults)
    calc_empty = DFTD3(
        method="PBE", damping="d3bj", realspace_cutoff={}, cache_api=True
    )
    atoms.calc = calc_empty
    energy_empty = atoms.get_potential_energy()

    # Empty dict should behave like no cutoff override (same as default)
    assert energy_empty == approx(energy_default, abs=thr)


@pytest.mark.skipif(ase is None, reason="requires ase")
def test_ase_tpssd3():
    thr = 1.0e-6

    forces = np.array(
        [
            [+1.21877896e-03, +2.15213419e-03, -1.16587717e-02],
            [-5.83752484e-04, +9.03161776e-03, +7.80662422e-03],
            [-4.29910142e-03, +4.66684847e-03, -4.56953656e-03],
            [-1.17172870e-03, -8.37195662e-03, +1.62139446e-02],
            [-6.96988440e-03, -5.09878876e-03, -1.75893974e-03],
            [+1.04161844e-02, -1.69289311e-02, -2.74875712e-03],
            [+5.52837716e-03, +3.32800644e-03, -1.02268960e-03],
            [-5.07390896e-03, +6.03647311e-03, -5.99714947e-05],
            [+3.34598807e-03, +3.95742085e-03, -5.56259794e-04],
            [-2.11534607e-03, +2.76307015e-03, -4.68480172e-04],
            [+4.31430409e-03, -1.56090240e-03, -6.80734794e-04],
            [-4.60991066e-03, +2.50079580e-05, -4.96427929e-04],
        ]
    )

    atoms = molecule("C2H6CHOH")
    atoms.calc = DFTD3(method="TPSS", damping="d3zero")

    assert atoms.get_potential_energy() == approx(-0.14270927858918006, abs=thr)
    assert atoms.get_forces() == approx(forces, abs=thr)

    atoms.calc = DFTD3(method="TPSS", damping="d3zero").add_calculator(EMT())
    assert atoms.get_potential_energy() == approx(4.963374535419293, abs=thr)
    energies = [calc.get_potential_energy() for calc in get_calcs(atoms.calc)]
    assert energies == approx([-0.14270927858918006, 5.106083814008478], abs=thr)
