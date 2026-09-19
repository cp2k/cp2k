"""Compare persistent ASE shell and libcp2k calculations on identical H2 points."""

# SPDX-License-Identifier: GPL-2.0-or-later

import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import platform
import time

import ase
from ase import Atoms
from ase.calculators.cp2k import CP2K as ShellCalculator
import numpy as np

from cp2k import CP2K, input_to_string
from cp2k.ase import CP2KCalculator


def calculation_input():
    return {
        "GLOBAL": {"PRINT_LEVEL": "SILENT"},
        "FORCE_EVAL": {
            "METHOD": "Quickstep",
            "DFT": {
                "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                "MGRID": {"CUTOFF": 200},
                "QS": {"EPS_DEFAULT": 1e-10},
                "SCF": {
                    "EPS_SCF": 1e-9,
                    "MAX_SCF": 100,
                    "OT": {"MINIMIZER": "DIIS"},
                    "PRINT": {"RESTART": {"_": "OFF"}},
                },
                "XC": {"XC_FUNCTIONAL": {"_": "PADE"}},
            },
            "SUBSYS": {
                "KIND": {
                    "_": "H",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    "POTENTIAL": "GTH-PADE-q1",
                }
            },
        },
    }


@contextmanager
def direct_calculator(library, inp):
    with CP2K(library=library) as runtime:
        with CP2KCalculator(runtime, inp, label="direct") as calculator:
            yield calculator, runtime.version


@contextmanager
def shell_calculator(command, inp):
    # Disable ASE-generated physical defaults; both backends get this template.
    with ShellCalculator(
        command=command,
        inp=input_to_string(inp),
        label="shell",
        basis_set=None,
        basis_set_file=None,
        potential_file=None,
        pseudo_potential=None,
        cutoff=None,
        max_scf=None,
        xc=None,
        force_eval_method=None,
        print_level=None,
        poisson_solver=None,
        stress_tensor=False,
    ) as calculator:
        yield calculator, None


def measure(factory, steps):
    atoms = Atoms("H2", positions=[[3.6, 4, 4], [4.4, 4, 4]], cell=[8] * 3, pbc=True)
    energies, forces, elapsed = [], [], []
    start = time.perf_counter()
    with factory as (calculator, version):
        setup_seconds = time.perf_counter() - start
        atoms.calc = calculator
        for point in range(steps + 1):
            # Changing every point prevents ASE's result cache from hiding work.
            atoms.positions[1, 0] = 4.4 + point * 0.005
            start = time.perf_counter()
            forces.append(atoms.get_forces().tolist())
            energies.append(atoms.get_potential_energy())
            elapsed.append(time.perf_counter() - start)
    return {
        "setup_seconds": setup_seconds,
        "first_evaluation_seconds": elapsed[0],
        "warm_evaluation_seconds": elapsed[1:],
        "warm_median_seconds": float(np.median(elapsed[1:])),
        "energies_eV": energies,
        "forces_eV_per_angstrom": forces,
        "cp2k_version": version,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--command", required=True, help="Matching CP2K executable with -s"
    )
    parser.add_argument("--library", default=os.environ.get("CP2K_LIBRARY"))
    parser.add_argument("--steps", type=int, default=12)
    parser.add_argument("--output", type=Path, default=Path("ase-comparison.json"))
    args = parser.parse_args()
    if args.steps < 1:
        parser.error("--steps must be positive")
    inp = calculation_input()
    # Finish the subprocess backend before initializing MPI in this process.
    shell = measure(shell_calculator(args.command, inp), args.steps)
    direct = measure(direct_calculator(args.library, inp), args.steps)
    np.testing.assert_allclose(
        direct["energies_eV"], shell["energies_eV"], rtol=0, atol=1e-6
    )
    np.testing.assert_allclose(
        direct["forces_eV_per_angstrom"],
        shell["forces_eV_per_angstrom"],
        rtol=0,
        atol=1e-5,
    )
    report = {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "ase": ase.__version__,
        "numpy": np.__version__,
        "command": args.command,
        "library": args.library,
        "threads": {
            key: os.environ.get(key)
            for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
        },
        "input": inp,
        "shell": shell,
        "direct": direct,
        "max_energy_difference_eV": float(
            np.max(np.abs(np.array(direct["energies_eV"]) - shell["energies_eV"]))
        ),
        "max_force_difference_eV_per_angstrom": float(
            np.max(
                np.abs(
                    np.array(direct["forces_eV_per_angstrom"])
                    - shell["forces_eV_per_angstrom"]
                )
            )
        ),
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Matched energies/forces; wrote {args.output}")


if __name__ == "__main__":
    main()
