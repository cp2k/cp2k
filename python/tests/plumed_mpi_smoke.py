"""Run on one and two ranks in separate directories; compare argon-1.cell/ener."""

# SPDX-License-Identifier: GPL-2.0-or-later

from pathlib import Path

from mpi4py import MPI
import numpy as np

from cp2k import CP2K
from cp2k._units import BOHR_TO_ANGSTROM, HARTREE_TO_JOULE

comm = MPI.COMM_WORLD
source = Path(__file__).resolve().parents[1] / "examples" / "argon.inp"
with CP2K(comm=comm) as runtime:
    with runtime.create_force_env(
        input_file=source, output_file="reference.out"
    ) as env:
        reference = env.calculate(stress=True)
    inp = source.read_text().replace("PRINT_LEVEL LOW", "PRINT_LEVEL LOW\n RUN_TYPE MD")
    inp = inp.replace("&COORD", "&VELOCITY\n 0 0 0\n 0 0 0\n &END VELOCITY\n &COORD")
    inp += """
&MOTION
 &FREE_ENERGY
  &METADYN
   USE_PLUMED T
   PLUMED_INPUT_FILE plumed.dat
  &END
 &END
 &MD
  ENSEMBLE NPE_I
  STEPS 5
  TIMESTEP 0.1
  TEMPERATURE 300
  &BAROSTAT
   PRESSURE 1
   TIMECON 100
  &END
 &END
 &PRINT
  &STRESS ON
  &END
  &CELL ON
  &END
 &END
&END
"""
    if comm.rank == 0:
        Path("plumed.dat").write_text("e: ENERGY\nb: BIASVALUE ARG=e\n")
    comm.Barrier()
    runtime.run_input(inp, output_file="md.out")
    comm.Barrier()
    energy = np.loadtxt("argon-1.ener")
    pressure = np.loadtxt("argon-1.stress")[0, 2:].reshape(3, 3)
    np.testing.assert_allclose(energy[0, 4], 2 * reference.energy, atol=1e-9)
    bar_factor = HARTREE_TO_JOULE / (BOHR_TO_ANGSTROM * 1e-10) ** 3 / 1e5
    np.testing.assert_allclose(pressure, 2 * reference.stress * bar_factor, rtol=2e-6)
    assert np.isfinite(energy).all()
if comm.rank == 0:
    print("PLUMED MPI: ENERGY bias, all virial components and variable-cell MD passed")
