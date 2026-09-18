"""Native PLUMED regressions; enable with CP2K_TEST_PLUMED=1."""

# SPDX-License-Identifier: GPL-2.0-or-later

from copy import deepcopy
import os

import numpy as np
import pytest

from cp2k._units import BOHR_TO_ANGSTROM, HARTREE_TO_JOULE

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not os.environ.get("CP2K_TEST_PLUMED"),
        reason="Set CP2K_TEST_PLUMED for a PLUMED-enabled library",
    ),
]


def run_md(runtime, inp, path, plumed, steps=0, restart=False):
    path.mkdir(exist_ok=restart)
    (path / "plumed.dat").write_text(plumed)
    inp = deepcopy(inp)
    inp["GLOBAL"].update(RUN_TYPE="MD", PROJECT="argon")
    inp["FORCE_EVAL"]["SUBSYS"]["VELOCITY"] = {"_lines": ["0 0 0", "0 0 0"]}
    inp["MOTION"] = {
        "FREE_ENERGY": {
            "METADYN": {"USE_PLUMED": True, "PLUMED_INPUT_FILE": "plumed.dat"}
        },
        "MD": {"ENSEMBLE": "NVE", "STEPS": steps, "TIMESTEP": 0.1, "TEMPERATURE": 300},
        "PRINT": {"FORCES": {"_": "ON"}, "STRESS": {"_": "ON"}},
    }
    if restart:
        inp["EXT_RESTART"] = {"RESTART_FILE_NAME": "argon-1.restart"}
    old = os.getcwd()
    try:
        os.chdir(path)
        runtime.run_input(inp, output_file="md.out")
    finally:
        os.chdir(old)
    energy = np.atleast_2d(np.loadtxt(path / "argon-1.ener"))
    forces = np.loadtxt(
        path / "argon-frc-1.xyz", skiprows=2, max_rows=2, usecols=(1, 2, 3)
    )
    pressure = np.atleast_2d(np.loadtxt(path / "argon-1.stress"))[0, 2:].reshape(3, 3)
    return energy, forces, pressure


def test_plumed_energy_bias(real_runtime, lj_input, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with real_runtime.create_force_env(lj_input, output_file="reference.out") as env:
        reference = env.calculate(stress=True)
    # BIASVALUE(ENERGY) must double the potential, physical forces AND virial.
    energy, forces, pressure = run_md(
        real_runtime,
        lj_input,
        tmp_path / "energy",
        "e: ENERGY\nb: BIASVALUE ARG=e\n",
        steps=2,
    )
    np.testing.assert_allclose(energy[0, 4], 2 * reference.energy, atol=1e-9)
    np.testing.assert_allclose(forces, 2 * reference.forces, rtol=2e-5, atol=1e-9)
    bar_factor = HARTREE_TO_JOULE / (BOHR_TO_ANGSTROM * 1e-10) ** 3 / 1e5
    np.testing.assert_allclose(pressure, 2 * reference.stress * bar_factor, rtol=2e-6)
    # A static bias is conservative; no previous bias may feed back into ENERGY.
    assert np.ptp(energy[:, 5]) < 1e-8


def test_plumed_distance_bias(real_runtime, lj_input, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with real_runtime.create_force_env(lj_input, output_file="reference.out") as env:
        reference = env.calculate(stress=True)
        delta = env.positions[1] - env.positions[0]
        volume = np.linalg.det(env.cell)
    distance = np.linalg.norm(delta)
    target, spring = 5.0, 0.001
    bias_energy = 0.5 * spring * (distance - target) ** 2
    force1 = -spring * (distance - target) * delta / distance
    bias_force = np.array([-force1, force1])
    bias_virial = np.outer(delta, force1)
    data = (
        "UNITS LENGTH=A ENERGY=Ha TIME=fs\nd: DISTANCE ATOMS=1,2\nb: RESTRAINT ARG=d AT=%.17g KAPPA=%.17g\n"
        % (target * BOHR_TO_ANGSTROM, spring / BOHR_TO_ANGSTROM**2)
    )
    energy, forces, pressure = run_md(
        real_runtime, lj_input, tmp_path / "distance", data
    )
    np.testing.assert_allclose(energy[0, 4], reference.energy + bias_energy, atol=1e-9)
    np.testing.assert_allclose(
        forces, reference.forces + bias_force, rtol=3e-5, atol=1e-9
    )
    bar_factor = HARTREE_TO_JOULE / (BOHR_TO_ANGSTROM * 1e-10) ** 3 / 1e5
    np.testing.assert_allclose(
        pressure, (reference.stress + bias_virial / volume) * bar_factor, rtol=3e-5
    )


def test_plumed_restart(real_runtime, lj_input, tmp_path):
    # Omit TEMP intentionally: well-tempered metadynamics must receive CP2K's kBT.
    data = "UNITS LENGTH=A ENERGY=Ha TIME=fs\nd: DISTANCE ATOMS=1,2\nm: METAD ARG=d SIGMA=0.2 HEIGHT=0.0001 PACE=1 BIASFACTOR=10 FILE=HILLS\n"
    run_md(real_runtime, lj_input, tmp_path / "continuous", data, steps=4)
    run_md(real_runtime, lj_input, tmp_path / "split", data, steps=2)
    before = np.atleast_2d(np.loadtxt(tmp_path / "split/HILLS"))
    run_md(real_runtime, lj_input, tmp_path / "split", data, steps=2, restart=True)
    after = np.atleast_2d(np.loadtxt(tmp_path / "split/HILLS"))
    continuous = np.atleast_2d(np.loadtxt(tmp_path / "continuous/HILLS"))
    assert len(before) == 2 and len(after) == 4
    np.testing.assert_array_equal(after[:2], before)
    np.testing.assert_allclose(after, continuous, rtol=1e-10, atol=1e-12)
    assert np.all(np.diff(after[:, 0]) > 0)  # no repeated hill at the restart step
    split_energy = np.loadtxt(tmp_path / "split/argon-1.ener")[-1, 2:6]
    continuous_energy = np.loadtxt(tmp_path / "continuous/argon-1.ener")[-1, 2:6]
    np.testing.assert_allclose(split_energy, continuous_energy, atol=1e-8, rtol=1e-8)
