# SPDX-License-Identifier: GPL-2.0-or-later

from copy import deepcopy

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import CalculationFailed

from cp2k._units import BOHR_TO_ANGSTROM as Bohr, HARTREE_TO_EV as Hartree
from cp2k.ase import CP2KCalculator


def test_ase_rejects_unconverged_result(runtime, fake_library, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0.8, 0, 0]], cell=[8] * 3)
    with CP2KCalculator(runtime, {}) as calc:
        atoms.calc = calc
        atoms.get_forces()
        atoms.positions[1, 0] += 0.1
        fake_library.scf_status = 0
        with pytest.raises(CalculationFailed, match="SCF did not converge"):
            atoms.get_forces()
        assert calc.results == {}


def test_ase_conversion_and_caching(runtime, fake_library, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = Atoms(
        "H2",
        positions=[[0, 0, 0], [0.8, 0, 0]],
        cell=[[5, 0, 0], [0.4, 6, 0], [0.2, 0.3, 7]],
        pbc=True,
    )
    inp = {"force_eval": {"subsys": {"kind": {"_": "H"}}}}
    saved = deepcopy(inp)
    with CP2KCalculator(runtime, inp) as calc:
        atoms.calc = calc
        assert atoms.get_potential_energy() == -Hartree
        np.testing.assert_allclose(atoms.get_forces(), Hartree / Bohr)
        np.testing.assert_allclose(fake_library.positions, atoms.positions / Bohr)
        np.testing.assert_allclose(fake_library.cell, atoms.cell.array / Bohr)
        assert (
            sum(name == "cp2k_calc_energy_force" for name, _ in fake_library.calls) == 1
        )
        atoms.positions[1, 0] += 0.1
        atoms.get_forces()
        assert (
            sum(name == "cp2k_calc_energy_force" for name, _ in fake_library.calls) == 2
        )
        assert fake_library.counter == 1
        atoms.pbc = False
        atoms.get_potential_energy()
        assert fake_library.counter == 2
        calc.set(inp={"FORCE_EVAL": {"METHOD": "FIST"}})
        atoms.get_potential_energy()
        assert fake_library.counter == 3
        np.testing.assert_allclose(
            atoms.get_stress(),
            -fake_library.stress.flat[[0, 4, 8, 5, 2, 1]] * Hartree / Bohr**3,
        )
    assert inp == saved
    with pytest.raises(RuntimeError, match="closed"):
        atoms.get_potential_energy()
    assert not runtime._closed


def test_project_path_is_quoted(runtime, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = Atoms("H", cell=[8] * 3)
    with CP2KCalculator(
        runtime, {}, directory=tmp_path / "directory with spaces"
    ) as calc:
        project = calc._make_input(atoms)["GLOBAL"]["PROJECT"]
        assert project == '"directory with spaces/cp2k"'
    with CP2KCalculator(runtime, {}, label="x" * 81) as calc:
        with pytest.raises(ValueError, match="80 bytes"):
            calc._make_input(atoms)


@pytest.mark.parametrize(
    "bad", ["CELL", "COORD", "TOPOLOGY", "VELOCITY", "SHELL", "CORE"]
)
def test_ase_rejects_geometry_input(runtime, bad, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0.8, 0, 0]], cell=[8] * 3)
    with CP2KCalculator(runtime, {"FORCE_EVAL": {"SUBSYS": {bad: {}}}}) as calc:
        atoms.calc = calc
        with pytest.raises(ValueError, match="ASE owns"):
            atoms.get_potential_energy()


def test_invalid_atoms_and_parameters(runtime, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(TypeError, match="Unknown"):
        CP2KCalculator(runtime, {}, xc="PBE")
    with CP2KCalculator(runtime, {}) as calc:
        atoms = Atoms("H", calculator=calc)
        with pytest.raises(ValueError, match="cell"):
            atoms.get_potential_energy()
        atoms.cell = [8] * 3
        atoms.set_initial_magnetic_moments([1])
        with pytest.raises(ValueError, match="spin"):
            atoms.get_potential_energy()
