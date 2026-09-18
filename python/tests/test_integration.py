"""Real library tests; enable with CP2K_TEST_LIBRARY=/absolute/libcp2k.so."""

# SPDX-License-Identifier: GPL-2.0-or-later

from copy import deepcopy
import os

import numpy as np
import pytest
from ase import Atoms
from ase.optimize import BFGS
from ase.units import Bohr, Hartree

from cp2k import CP2K, input_to_string
from cp2k.ase import CP2KCalculator

pytestmark = pytest.mark.integration


def test_native_energy_force_and_cell(real_runtime, h2_input, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert "CP2K version" in real_runtime.version
    with real_runtime.create_force_env(h2_input, output_file="native.out") as env:
        assert env.natom == env.nparticle == 2
        positions = env.positions
        np.testing.assert_allclose(
            positions, np.array([[3.6, 4, 4], [4.4, 4, 4]]) / Bohr, rtol=1e-7
        )
        result = env.calculate()
        assert -1.3 < result.energy < -0.8
        np.testing.assert_allclose(result.forces.sum(axis=0), 0, atol=1e-7)
        delta = 1e-3
        plus = positions.copy()
        plus[1, 0] += delta
        env.positions = plus
        eplus = env.calculate(forces=False).energy
        minus = positions.copy()
        minus[1, 0] -= delta
        env.positions = minus
        eminus = env.calculate(forces=False).energy
        np.testing.assert_allclose(
            result.forces[1, 0], -(eplus - eminus) / (2 * delta), atol=2e-5
        )
        # A non-symmetric cell catches C/Fortran transposition errors.
        cell = np.array([[15.0, 0.0, 0.0], [0.4, 16.0, 0.0], [0.3, 0.2, 17.0]])
        env.cell = cell
        np.testing.assert_allclose(env.cell, cell)
        env.positions = positions
        assert np.isfinite(env.calculate().energy)


@pytest.mark.parametrize("units", ["metal", "real"])
def test_lammps_native(real_runtime, lj_input, tmp_path, monkeypatch, units):
    if not os.environ.get("CP2K_TEST_LAMMPS"):
        pytest.skip("Set CP2K_TEST_LAMMPS for a compatible LAMMPS library")
    from lammps import lammps
    from cp2k.lammps import ExternalForce
    from cp2k._units import BOHR_TO_ANGSTROM, HARTREE_TO_EV, HARTREE_TO_KCALMOL

    monkeypatch.chdir(tmp_path)
    energy_factor = HARTREE_TO_EV if units == "metal" else HARTREE_TO_KCALMOL
    with real_runtime.create_force_env(lj_input, output_file="lammps.out") as env:
        reference = env.calculate(stress=True)
        lmp = lammps(
            comm=real_runtime._comm, cmdargs=["-log", "none", "-screen", "none"]
        )
        try:
            lmp.commands_string(f"""
units {units}
atom_style atomic
boundary p p p
region box prism 0 20 0 21 0 22 1 2 3
create_box 1 box
create_atoms 1 single 7.2 4.8 4.5
create_atoms 1 single 4 4 4
mass 1 39.948
pair_style zero 8.0
pair_coeff * *
compute cp_pressure all pressure NULL virial
thermo_style custom step pe c_cp_pressure[1] c_cp_pressure[2] c_cp_pressure[3] c_cp_pressure[4] c_cp_pressure[5] c_cp_pressure[6]
thermo_modify norm no
""")
            with ExternalForce(lmp, env, atom_ids=[2, 1]) as callback:
                callback.run(0)
                np.testing.assert_allclose(
                    lmp.get_thermo("pe"), reference.energy * energy_factor, atol=1e-10
                )
                tags = lmp.numpy.extract_atom("id")[:2]
                indices = np.array([1 if tag == 1 else 0 for tag in tags])
                np.testing.assert_allclose(
                    lmp.numpy.extract_atom("f")[:2],
                    reference.forces[indices] * energy_factor / BOHR_TO_ANGSTROM,
                    rtol=1e-8,
                )
                pressure = lmp.numpy.extract_compute("cp_pressure", 0, 1).copy()
                expected = (
                    reference.virial.flat[[0, 4, 8, 1, 2, 5]]
                    * energy_factor
                    / np.linalg.det(env.cell * BOHR_TO_ANGSTROM)
                    * lmp.extract_global("nktv2p")
                )
                np.testing.assert_allclose(pressure, expected, rtol=1e-8, atol=1e-10)
                callback.command("change_box all x scale 1.01 remap")
                callback.run(0)
                assert env.cell[0, 0] * BOHR_TO_ANGSTROM == pytest.approx(20.2)
                callback.command(
                    "velocity all create 10.0 731 mom yes rot no dist gaussian"
                )
                callback.command("fix thermostat all npt temp 10 10 100 iso 0 0 1000")
                callback.command(
                    "timestep " + ("0.0001" if units == "metal" else "0.1")
                )
                callback.run(3)
                assert np.isfinite(lmp.get_thermo("pe"))
                callback.command("unfix thermostat")
            assert not lmp.has_id("fix", "cp2k")
            with ExternalForce(lmp, env, stress=False) as callback:
                with pytest.raises(RuntimeError, match="callback failed"):
                    callback.command("change_box all x scale 1.01 remap")
                    callback.run(0)
            assert np.isfinite(env.calculate().energy)
        finally:
            lmp.close()
    assert not list(tmp_path.glob("cp2k-python-*.inp"))


def test_ase_optimization(real_runtime, h2_input, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with real_runtime.create_force_env(h2_input, output_file="reference.out") as env:
        reference = env.calculate()
    inp = deepcopy(h2_input)
    del inp["FORCE_EVAL"]["SUBSYS"]["CELL"]
    del inp["FORCE_EVAL"]["SUBSYS"]["COORD"]
    del inp["GLOBAL"]  # Exercise generated PROJECT, including paths with spaces.
    atoms = Atoms("H2", positions=[[3.6, 4, 4], [4.4, 4, 4]], cell=[8] * 3, pbc=True)
    with CP2KCalculator(real_runtime, inp, directory="ase output", label="h2") as calc:
        atoms.calc = calc
        e0 = atoms.get_potential_energy()
        np.testing.assert_allclose(e0, reference.energy * Hartree, atol=1e-6)
        np.testing.assert_allclose(
            atoms.get_forces(), reference.forces * Hartree / Bohr, atol=1e-5
        )
        optimizer = BFGS(atoms, logfile=None)
        optimizer.run(fmax=0.1, steps=3)
        assert atoms.get_potential_energy() < e0


def test_run_input_and_file_input(real_runtime, h2_input, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    h2_input["GLOBAL"].update(RUN_TYPE="ENERGY_FORCE", PRINT_LEVEL="LOW")
    source = tmp_path / "h2.inp"
    source.write_text(input_to_string(h2_input))
    real_runtime.run_input(input_file=source, output_file="run.out")
    output = (tmp_path / "run.out").read_text()
    assert "ENERGY|" in output
    assert "PROGRAM ENDED AT" in output
    real_runtime.run_input(input_file=source, output_file="run.out")
    assert (tmp_path / "run.out").read_text().startswith(output)
    assert (tmp_path / "run.out").stat().st_size > len(output)
    with real_runtime.create_force_env(
        input_file=source, output_file="file.out"
    ) as env:
        assert np.isfinite(env.calculate().energy)


def test_sequential_environments(real_runtime, h2_input, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    second_input = deepcopy(h2_input)
    second_input["GLOBAL"]["PROJECT"] = "second-h2"
    second_input["FORCE_EVAL"]["SUBSYS"]["COORD"]["_lines"][1] = "H 4.5 4 4"
    with real_runtime.create_force_env(h2_input, output_file="first.out") as first:
        reference = first.calculate()
        # Native teardown frees shared integral tables; reject overlapping
        # environments before entering CP2K, rather than risking a segfault.
        with pytest.raises(RuntimeError, match="current force environment"):
            real_runtime.create_force_env(second_input, output_file="second.out")
        assert np.isfinite(first.calculate().energy)
    with real_runtime.create_force_env(
        second_input, output_file="second.out"
    ) as second:
        assert abs(second.calculate().energy - reference.energy) > 1e-4
    with real_runtime.create_force_env(
        h2_input, output_file="repeated.out"
    ) as repeated:
        result = repeated.calculate()
        np.testing.assert_allclose(result.energy, reference.energy, atol=1e-8)
        np.testing.assert_allclose(result.forces, reference.forces, atol=1e-7)


def test_native_md(real_runtime, h2_input, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    h2_input["GLOBAL"].update(RUN_TYPE="MD", PRINT_LEVEL="LOW", PROJECT="python-md")
    h2_input["MOTION"] = {
        "MD": {"ENSEMBLE": "NVE", "STEPS": 2, "TIMESTEP": 0.1, "TEMPERATURE": 100}
    }
    real_runtime.run_input(h2_input, output_file="md.out")
    output = (tmp_path / "md.out").read_text()
    assert "PROGRAM ENDED AT" in output
    trajectory = np.loadtxt(tmp_path / "python-md-1.ener")
    np.testing.assert_array_equal(trajectory[:, 0], [0, 1, 2])
    assert np.isfinite(trajectory).all()


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


@pytest.mark.parametrize("platform_name", ["Reference", "CPU"])
def test_openmm_native(real_runtime, lj_input, tmp_path, monkeypatch, platform_name):
    openmm = pytest.importorskip("openmm", minversion="8.6.1")
    from openmm import unit
    from cp2k.openmm import create_force
    from cp2k._units import BOHR_TO_NM, HARTREE_TO_KJMOL

    monkeypatch.chdir(tmp_path)
    with real_runtime.create_force_env(lj_input, output_file="openmm.out") as env:
        reference = env.calculate()
        system = openmm.System()
        for _ in range(env.nparticle):
            system.addParticle(39.948)
        system.setDefaultPeriodicBoxVectors(*env.cell * BOHR_TO_NM)
        force = create_force(env, periodic=True)
        system.addForce(force)
        integrator = openmm.VerletIntegrator(0.0001)
        context = openmm.Context(
            system, integrator, openmm.Platform.getPlatformByName(platform_name)
        )
        context.setPositions(env.positions * BOHR_TO_NM)
        state = context.getState(getEnergy=True, getForces=True)
        np.testing.assert_allclose(
            state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole),
            reference.energy * HARTREE_TO_KJMOL,
            atol=2e-6,
        )
        np.testing.assert_allclose(
            state.getForces(asNumpy=True).value_in_unit(
                unit.kilojoules_per_mole / unit.nanometer
            ),
            reference.forces * HARTREE_TO_KJMOL / BOHR_TO_NM,
            rtol=2e-5,
            atol=2e-5,
        )
        context.setVelocitiesToTemperature(10, 17)
        integrator.step(3)
        assert np.isfinite(context.getState(getEnergy=True).getPotentialEnergy()._value)
        new_cell = env.cell * 1.01
        context.setPeriodicBoxVectors(*new_cell * BOHR_TO_NM)
        context.getState(getEnergy=True)
        np.testing.assert_allclose(env.cell, new_cell)
        del context, integrator
        # The caller retains ownership of both the environment and runtime.
        assert np.isfinite(env.calculate().energy)
