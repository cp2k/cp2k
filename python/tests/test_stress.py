# SPDX-License-Identifier: GPL-2.0-or-later

from copy import deepcopy

import numpy as np
from ase import Atoms
from ase.optimize import BFGS

from cp2k._units import BOHR_TO_ANGSTROM as Bohr, HARTREE_TO_EV as Hartree
from cp2k.ase import CP2KCalculator
import pytest


def test_stress_lifecycle(runtime, fake_library, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with runtime.create_force_env({"GLOBAL": {}}) as env:
        with pytest.raises(RuntimeError, match="stress=True"):
            _ = env.stress
        result = env.calculate(stress=True)
        np.testing.assert_array_equal(result.stress, fake_library.stress)
        np.testing.assert_allclose(
            result.virial, result.stress * np.linalg.det(env.cell)
        )
        env.positions = env.positions
        with pytest.raises(RuntimeError, match="stress=True"):
            _ = env.virial
        assert env.calculate(forces=False, stress=True).forces is None
        env.calculate(forces=False)
        with pytest.raises(RuntimeError, match="stress=True"):
            _ = env.stress
        fake_library.stress_available = 0
        with pytest.raises(RuntimeError, match="STRESS_TENSOR"):
            env.calculate(stress=True)
        with pytest.raises(RuntimeError):
            _ = env.potential_energy


def test_older_library(runtime, fake_library, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    del fake_library.cp2k_get_stress_tensor
    with runtime.create_force_env({"GLOBAL": {}}) as env:
        assert env.calculate().energy == -1
        with pytest.raises(RuntimeError, match="lacks"):
            env.calculate(stress=True)


@pytest.mark.integration
def test_native_stress_finite_strain(real_runtime, lj_input, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with real_runtime.create_force_env(lj_input, output_file="stress.out") as env:
        positions, cell = env.positions, env.cell
        result = env.calculate(stress=True)
        np.testing.assert_allclose(result.stress, result.stress.T, atol=1e-14)
        delta = 1e-5
        for i, j in ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)):
            energies = []
            for sign in (1, -1):
                deformation = np.eye(3)
                deformation[i, j] += sign * delta
                env.cell = cell @ deformation.T
                env.positions = positions @ deformation.T
                energies.append(env.calculate(forces=False).energy)
            derivative = (energies[0] - energies[1]) / (2 * delta)
            np.testing.assert_allclose(
                -derivative, result.virial[i, j], rtol=2e-5, atol=1e-10
            )


@pytest.mark.integration
def test_ase_cell_optimization(real_runtime, lj_input, tmp_path, monkeypatch):
    from ase.filters import FrechetCellFilter
    from cp2k._units import BOHR_TO_ANGSTROM

    monkeypatch.chdir(tmp_path)
    with real_runtime.create_force_env(lj_input, output_file="reference.out") as env:
        cell, positions = env.cell * BOHR_TO_ANGSTROM, env.positions * BOHR_TO_ANGSTROM
        reference = env.calculate(stress=True)
    inp = deepcopy(lj_input)
    for key in ("CELL", "COORD", "TOPOLOGY"):
        del inp["FORCE_EVAL"]["SUBSYS"][key]
    atoms = Atoms("Ar2", positions=positions, cell=cell, pbc=True)
    with CP2KCalculator(real_runtime, inp, output_file="ase-cell.out") as calc:
        atoms.calc = calc
        np.testing.assert_allclose(
            atoms.get_stress(voigt=False),
            -reference.stress * Hartree / Bohr**3,
            rtol=1e-4,
            atol=1e-9,
        )
        initial = atoms.get_potential_energy()
        BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-4, steps=5)
        assert atoms.get_potential_energy() < initial
        assert not np.allclose(atoms.cell.array, cell)
