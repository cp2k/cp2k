"""Run with mpiexec -n 2 python lammps_mpi_smoke.py in a dedicated directory."""

# SPDX-License-Identifier: GPL-2.0-or-later

from pathlib import Path
from contextlib import contextmanager
import json
import sys

from mpi4py import MPI
from lammps import lammps
import numpy as np

from cp2k import CP2K, SocketEnvironment
from cp2k.lammps import ExternalForce
from cp2k._units import BOHR_TO_ANGSTROM, HARTREE_TO_EV

comm = MPI.COMM_WORLD
source = Path(__file__).resolve().parents[1] / "examples" / "argon.inp"


@contextmanager
def environment():
    if len(sys.argv) == 3 and sys.argv[1] == "--socket-config":
        config = json.loads(Path(sys.argv[2]).read_text())
        with SocketEnvironment(**config, comm=comm) as env:
            yield env
    else:
        with CP2K(comm=comm) as runtime:
            with runtime.create_force_env(
                input_file=source, output_file="mpi-lammps.out"
            ) as env:
                yield env


def exercise(env):
    reference = env.calculate(stress=True)
    lmp = lammps(comm=comm, cmdargs=["-log", "none", "-screen", "none"])
    try:
        lmp.commands_string("""
units metal
atom_style atomic
boundary p p p
region box prism 0 20 0 21 0 22 1 2 3
create_box 1 box
create_atoms 1 single 7.2 4.8 4.5
create_atoms 1 single 4 4 4
mass 1 39.948
pair_style zero 8
pair_coeff * *
compute cp_pressure all pressure NULL virial
thermo_style custom step pe c_cp_pressure[1] c_cp_pressure[2] c_cp_pressure[3] c_cp_pressure[4] c_cp_pressure[5] c_cp_pressure[6]
thermo_modify norm no
""")
        with ExternalForce(lmp, env, atom_ids=[2, 1]) as callback:
            callback.run(0)
            np.testing.assert_allclose(
                lmp.get_thermo("pe"), reference.energy * HARTREE_TO_EV, rtol=1e-10
            )
            expected = (
                reference.virial.flat[[0, 4, 8, 1, 2, 5]]
                * HARTREE_TO_EV
                / np.linalg.det(env.cell * BOHR_TO_ANGSTROM)
                * lmp.extract_global("nktv2p")
            )
            np.testing.assert_allclose(
                lmp.numpy.extract_compute("cp_pressure", 0, 1), expected, rtol=1e-10
            )
            callback.command("velocity all create 10 731 mom yes rot no dist gaussian")
            callback.command("fix thermostat all npt temp 10 10 0.1 iso 0 0 1")
            callback.command("timestep 0.0001")
            callback.run(5)
            assert np.isfinite(lmp.get_thermo("pe"))
            callback.command("unfix thermostat")
            # A rank-local lifecycle failure must terminate collectively, not hang.
            if comm.rank == 0:
                callback._error = ValueError("deliberate rank-local validation failure")
            # Bypass command()'s caller-side check to exercise callback-side validation.
            lmp.command("run 0")
            try:
                callback.check()
            except RuntimeError:
                pass
            else:
                raise AssertionError("Expected a collective callback failure")
    finally:
        lmp.close()


with environment() as env:
    exercise(env)
assert not MPI.Is_finalized()
comm.Barrier()
if comm.rank == 0:
    print(
        "LAMMPS MPI: energy, virial normalization, NPT, empty ranks and failure propagation passed"
    )
