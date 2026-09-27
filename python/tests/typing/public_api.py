"""Static API contract checks; mypy checks this file without executing it."""

# SPDX-License-Identifier: GPL-2.0-or-later

from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from ase import Atoms

from cp2k import CP2K, CalculationResult, ForceEnvironment, SocketEnvironment
from cp2k.ase import CP2KCalculator


def public_api(runtime: CP2K, atoms: Atoms, remote: SocketEnvironment) -> None:
    env: ForceEnvironment = runtime.create_force_env(
        {"GLOBAL": {"PROJECT": "typing-check"}}, output_file=Path("typing.out")
    )
    version: str = runtime.version
    env.positions = [[0.0, 0.0, 0.0]]
    env.cell = np.eye(3)
    result: CalculationResult = env.calculate(forces=True)
    energy: float = result.energy
    forces: NDArray[np.float64] | None = result.forces
    scf: bool | None = result.scf_converged
    remote_result: CalculationResult = remote.calculate(stress=True)
    remote_forces: NDArray[np.float64] = remote.forces
    with CP2KCalculator(runtime, {}, output_file=Path("ase.out")) as calc:
        calc.calculate(atoms, properties=("energy", "forces"))

    # Each ignore must suppress a real error: strict mode rejects unused ignores.
    # This catches accidentally weakening a public signature or result to Any.
    CP2K(library=42)  # type: ignore[arg-type]
    runtime.create_force_env(output_file=42)  # type: ignore[arg-type]
    env.calculate(forces="yes")  # type: ignore[arg-type]
    remote.calculate(stress="yes")  # type: ignore[arg-type]
    CP2KCalculator(runtime, "not a mapping")  # type: ignore[arg-type]
    energy = version  # type: ignore[assignment]
    forces = energy  # type: ignore[assignment]
    scf = energy  # type: ignore[assignment]
    remote_forces = remote_result.energy  # type: ignore[assignment]
