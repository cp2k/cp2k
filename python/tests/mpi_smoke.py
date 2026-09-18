"""Run via mpiexec -n 2 python mpi_smoke.py [world|split|implicit].

This deliberately runs outside pytest: MPI must own initialization before CP2K.
Use a dedicated working directory; output files are retained for inspection.
"""

# SPDX-License-Identifier: GPL-2.0-or-later

from pathlib import Path
import sys

from mpi4py import MPI
import numpy as np

from cp2k import CP2K

mode = sys.argv[1] if len(sys.argv) > 1 else "world"
if mode not in ("world", "split", "implicit"):
    raise ValueError("Expected world, split or implicit")
world = MPI.COMM_WORLD
# Only rank 0 participates in CP2K in split mode. An accidental use of WORLD
# inside the runtime would deadlock; run the smoke test with a timeout.
comm = (
    world.Split(0 if world.rank == 0 else MPI.UNDEFINED) if mode == "split" else world
)
if comm != MPI.COMM_NULL:
    with CP2K(comm=None if mode == "implicit" else comm) as cp:
        inp = Path(__file__).resolve().parents[1] / "examples" / "h2.inp"
        with cp.create_force_env(input_file=inp, output_file=f"mpi-{mode}.out") as env:
            result = env.calculate()
            energies = comm.allgather(result.energy)
            np.testing.assert_allclose(energies, result.energy, atol=1e-12)
            assert -1.3 < result.energy < -0.8
            np.testing.assert_allclose(result.forces.sum(axis=0), 0, atol=1e-7)
        cp.run_input(input_file=inp, output_file=f"mpi-{mode}-run.out")
        with cp.create_force_env(
            input_file=inp, output_file=f"mpi-{mode}-after.out"
        ) as env:
            np.testing.assert_allclose(env.calculate().energy, result.energy, atol=1e-8)
    assert not MPI.Is_finalized()
    if comm.rank == 0:
        assert "PROGRAM ENDED AT" in Path(f"mpi-{mode}-run.out").read_text()
    if mode == "split":
        comm.Free()
world.Barrier()
if world.rank == 0:
    print(f"MPI {mode}: energy, forces, run_input and MPI ownership passed")
