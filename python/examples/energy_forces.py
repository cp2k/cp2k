"""Run with CP2K_LIBRARY and CP2K_DATA_DIR pointing to your CP2K installation."""

# SPDX-License-Identifier: GPL-2.0-or-later

from pathlib import Path

from cp2k import CP2K

with CP2K() as cp:
    print(cp.version)
    with cp.create_force_env(
        input_file=Path(__file__).with_name("h2.inp"), output_file="h2.out"
    ) as system:
        first = system.calculate()
        print("Energy [hartree]:", first.energy)
        print("Forces [hartree/bohr]:\n", first.forces)
        positions = system.positions
        positions[1, 0] += 0.05
        system.positions = positions
        print("Displaced energy [hartree]:", system.calculate().energy)
