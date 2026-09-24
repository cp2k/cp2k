"""Real library tests; enable with CP2K_TEST_LIBRARY=/absolute/libcp2k.so."""

# SPDX-License-Identifier: GPL-2.0-or-later

from copy import deepcopy
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import numpy as np
import pytest
from ase import Atoms
from ase.optimize import BFGS

from cp2k import CP2K, SCFConvergenceError, input_to_string
from cp2k._units import BOHR_TO_ANGSTROM as Bohr, HARTREE_TO_EV as Hartree
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
        assert result.scf_converged is True
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
    assert not list(tmp_path.glob("cp2k-python-*.inp"))


@pytest.mark.parametrize("outer", [False, True])
def test_native_scf_failure(real_runtime, h2_input, tmp_path, monkeypatch, outer):
    monkeypatch.chdir(tmp_path)
    scf = h2_input["FORCE_EVAL"]["DFT"]["SCF"]
    scf.update(MAX_SCF=1, EPS_SCF=1e-30, IGNORE_CONVERGENCE_FAILURE=True)
    if outer:
        scf["OUTER_SCF"] = {"MAX_SCF": 1, "EPS_SCF": 1e-30}
    with real_runtime.create_force_env(h2_input, output_file="failed.out") as env:
        with pytest.raises(SCFConvergenceError):
            env.calculate()
        assert env.scf_converged is False
        with pytest.raises(RuntimeError, match="calculate"):
            _ = env.forces
        result = env.calculate(forces=False, check_convergence=False)
        assert np.isfinite(result.energy)
        assert result.scf_converged is False


def test_native_no_scf(real_runtime, h2_input, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    h2_input["FORCE_EVAL"]["DFT"]["SCF"]["MAX_SCF"] = 0
    with real_runtime.create_force_env(h2_input, output_file="no-scf.out") as env:
        assert env.calculate(forces=False).scf_converged is None


@pytest.mark.parametrize("target, expected", [(0.0, True), (0.5, False)])
def test_native_cdft_status(
    real_runtime, h2_input, tmp_path, monkeypatch, target, expected
):
    monkeypatch.chdir(tmp_path)
    dft = h2_input["FORCE_EVAL"]["DFT"]
    dft["SCF"]["IGNORE_CONVERGENCE_FAILURE"] = True
    dft["SCF"]["OUTER_SCF"] = {"MAX_SCF": 5, "EPS_SCF": 1e-9}
    dft["QS"]["CDFT"] = {
        "TYPE_OF_CONSTRAINT": "BECKE",
        "STRENGTH": 0.0,
        "TARGET": target,
        "ATOM_GROUP": {"ATOMS": [1, 2], "COEFF": [1, -1]},
        "BECKE_CONSTRAINT": {"IN_MEMORY": True},
        "OUTER_SCF": {
            "EPS_SCF": 1e-4,
            "MAX_SCF": 1,
            "TYPE": "CDFT_CONSTRAINT",
            "OPTIMIZER": "BISECT",
        },
    }
    with real_runtime.create_force_env(h2_input, output_file="cdft.out") as env:
        result = env.calculate(forces=False, check_convergence=False)
        assert result.scf_converged is expected


@pytest.mark.parametrize("solver", ["diagonalization", "outer_scf"])
def test_native_scf_solvers(real_runtime, h2_input, tmp_path, monkeypatch, solver):
    monkeypatch.chdir(tmp_path)
    scf = h2_input["FORCE_EVAL"]["DFT"]["SCF"]
    if solver == "diagonalization":
        del scf["OT"]
        scf["DIAGONALIZATION"] = {}
    else:
        scf["OUTER_SCF"] = {"MAX_SCF": 5, "EPS_SCF": 1e-9}
    with real_runtime.create_force_env(h2_input, output_file="solver.out") as env:
        assert env.calculate().scf_converged is True


def test_ase_shell_comparison(tmp_path):
    library = os.environ.get("CP2K_TEST_LIBRARY")
    executable = os.environ.get("CP2K_TEST_EXECUTABLE")
    if not library or not executable:
        pytest.skip("Set CP2K_TEST_LIBRARY and CP2K_TEST_EXECUTABLE")
    script = Path(__file__).resolve().parents[1] / "examples" / "compare_ase.py"
    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            "--library",
            library,
            "--command",
            f"{shlex.quote(executable)} -s",
            "--steps",
            "2",
        ],
        cwd=tmp_path,
        check=False,
        timeout=120,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    report = json.loads((tmp_path / "ase-comparison.json").read_text())
    assert report["max_energy_difference_eV"] < 1e-6
    assert report["max_force_difference_eV_per_angstrom"] < 1e-5
    for backend in ("shell", "direct"):
        assert len(report[backend]["warm_evaluation_seconds"]) == 2


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
