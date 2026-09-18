"""Optional OpenMM PythonForce adapter for a caller-owned CP2K environment."""

# SPDX-License-Identifier: GPL-2.0-or-later

from ._units import BOHR_TO_NM, HARTREE_TO_KJMOL
from .library import _array


class _Computation:
    def __init__(self, environment, periodic):
        environment._check()
        if environment.natom != environment.nparticle:
            raise ValueError("OpenMM requires one CP2K particle per atom")
        comm = environment.communicator
        if comm is not None and comm.size != 1:
            raise ValueError(
                "OpenMM needs a single-rank caller; use SocketEnvironment "
                "for a multi-rank CP2K server"
            )
        self.environment = environment
        self.periodic = bool(periodic)
        self.nparticle = environment.nparticle

    def __call__(self, state):
        env = self.environment
        env._check()
        from openmm import unit

        positions = _array(
            state.getPositions(asNumpy=True).value_in_unit(unit.nanometer),
            (self.nparticle, 3),
            "OpenMM positions",
        )
        if self.periodic:
            env.cell = (
                state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
                / BOHR_TO_NM
            )
        env.positions = positions / BOHR_TO_NM
        result = env.calculate()
        return result.energy * HARTREE_TO_KJMOL, result.forces * (
            HARTREE_TO_KJMOL / BOHR_TO_NM
        )

    def __reduce__(self):
        raise TypeError(
            "A live CP2K environment cannot be serialized. Recreate the runtime, "
            "environment, and PythonForce when restoring an OpenMM simulation."
        )


def create_force(environment, *, periodic):
    """Create an OpenMM force for the entire CP2K system, in identical atom order.

    Requires OpenMM >= 8.6.1 (PythonForce callbacks on the calling thread).
    Call OpenMM from Python's main thread, with a single-rank caller. To use
    multi-rank CP2K, pass a SocketEnvironment connected to an MPI server.
    Keep the caller-owned environment/runtime alive until all Contexts using
    this force are destroyed. Do not use simultaneous Contexts with it.

    Match ``periodic`` to CP2K SUBSYS/CELL/PERIODIC (XYZ or NONE). For periodic
    calculations the current OpenMM box is passed at every evaluation, including
    barostat trial moves. Otherwise CP2K retains its input cell. This supplies
    the potential supplied by the environment. For native CP2K QM/MM, configure
    METHOD QMMM and let CP2K supply the whole potential. Alternatively pass a
    SubtractiveQMMM correction while retaining the host's full MM potential.
    Do not add a whole-system CP2K potential to the same host MM interactions.
    CP2K state is not part of OpenMM's XML/checkpoints: rebuild the adapter and
    use CP2K's own wavefunction restart mechanism if needed.
    """
    import openmm

    version = tuple(int(part) for part in openmm.__version__.split(".")[:3])
    if version < (8, 6, 1) or not hasattr(openmm, "PythonForce"):
        raise RuntimeError("The CP2K adapter requires OpenMM >= 8.6.1")
    force = openmm.PythonForce(_Computation(environment, periodic))
    force.setUsesPeriodicBoundaryConditions(bool(periodic))
    force.setName("CP2K")
    return force
