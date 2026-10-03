"""Run on one and two ranks in separate directories; compare argon-1.cell/ener."""

# SPDX-License-Identifier: GPL-2.0-or-later

from copy import deepcopy
from pathlib import Path

from mpi4py import MPI
import numpy as np

from cp2k import CP2K
from plumed_input import plumed_input

comm = MPI.COMM_WORLD
inp = plumed_input()
inp["GLOBAL"].update(RUN_TYPE="MD", PROJECT="reference")
inp["FORCE_EVAL"]["SUBSYS"]["VELOCITY"] = {"_lines": ["0 0 0", "0 0 0"]}
inp["MOTION"] = {
    "MD": {"ENSEMBLE": "NVE", "STEPS": 0, "TIMESTEP": 0.1, "TEMPERATURE": 300},
    "PRINT": {"STRESS": {"_": "ON"}, "CELL": {"_": "ON"}},
}
with CP2K(comm=comm) as runtime:
    # Compare against native unbiased output, independent of the new stress API.
    runtime.run_input(inp, output_file="reference.out")
    comm.Barrier()
    reference_energy = np.atleast_2d(np.loadtxt("reference-1.ener"))
    reference_pressure = np.atleast_2d(np.loadtxt("reference-1.stress"))[0, 2:]
    biased = deepcopy(inp)
    biased["GLOBAL"]["PROJECT"] = "argon"
    biased["MOTION"]["MD"].update(
        ENSEMBLE="NPE_I", STEPS=5, BAROSTAT={"PRESSURE": 1, "TIMECON": 100}
    )
    biased["MOTION"]["FREE_ENERGY"] = {
        "METADYN": {"USE_PLUMED": True, "PLUMED_INPUT_FILE": "plumed.dat"}
    }
    if comm.rank == 0:
        Path("plumed.dat").write_text("e: ENERGY\nb: BIASVALUE ARG=e\n")
    comm.Barrier()
    runtime.run_input(biased, output_file="md.out")
    comm.Barrier()
    energy = np.loadtxt("argon-1.ener")
    pressure = np.loadtxt("argon-1.stress")[0, 2:]
    np.testing.assert_allclose(energy[0, 4], 2 * reference_energy[0, 4], atol=1e-9)
    np.testing.assert_allclose(pressure, 2 * reference_pressure, rtol=2e-6)
    assert np.isfinite(energy).all()
if comm.rank == 0:
    print("PLUMED MPI: ENERGY bias, all virial components and variable-cell MD passed")
