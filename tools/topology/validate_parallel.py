"""Compare native topology across MPI/OpenMP layouts and test the MPI adapter."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np
import z2pack

from cp2k_z2pack import CP2KSystem
from validate_cp2k import ROOT
from validate_stanene import compare as compare_stanene


def run(binary, directory, text, ranks, threads):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "input.inp").write_text(text)
    command = ["mpiexec", "-n", str(ranks), str(binary), "-i", "input.inp"]
    with (directory / "run.log").open("w") as output:
        subprocess.run(
            command,
            cwd=directory,
            stdout=output,
            stderr=subprocess.STDOUT,
            check=True,
            env={
                **os.environ,
                "OMP_NUM_THREADS": str(threads),
                "OPENBLAS_NUM_THREADS": "1",
                "CP2K_DATA_DIR": str(ROOT / "data"),
            },
        )
    log = (directory / "run.log").read_text()
    actual_ranks = int(
        re.search(r"Total number of message passing processes\s+(\d+)", log)[1]
    )
    assert actual_ranks == ranks
    return log


def phase_error(a, b):
    assert a.shape == b.shape
    return float(
        max(
            min(
                np.max(
                    np.abs(np.exp(2j * np.pi * x) - np.exp(2j * np.pi * np.roll(y, s)))
                )
                for s in range(len(x))
            )
            for x, y in zip(a, b)
        )
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("workdir", type=Path)
    parser.add_argument(
        "--stanene",
        action="store_true",
        help="Also run the full nontrivial native surface on two ranks",
    )
    args = parser.parse_args()
    binary, root = args.binary.resolve(), args.workdir.resolve()
    examples = Path(__file__).with_name("examples")
    base = (examples / "neon-soc.inp").read_text()
    # A complete scalar basis avoids rank-dependent choices inside a truncated,
    # degenerate virtual manifold when comparing independently diagonalized runs.
    base = base.replace("ADDED_MOS 1", "ADDED_MOS -1")
    base = base.replace(
        "EXCLUDE_BANDS 9 10", "EXCLUDE_BANDS " + " ".join(map(str, range(9, 27)))
    )
    base = base.replace("EPS_SCF 1.0E-9", "EPS_SCF 1.0E-12")
    report = {"layouts": []}
    reference = None
    for ranks, threads in [(1, 1), (2, 1), (2, 2), (4, 1)]:
        directory = root / f"neon-{ranks}r-{threads}t"
        log = run(binary, directory, base, ranks, threads)
        assert "Converged Z2 invariant: 0" in log
        data = np.loadtxt(directory / "neon.wilson", ndmin=2)[:, 2:]
        energy = float(re.findall(r"Total energy:\s+([-\d.]+)", log)[-1])
        if reference is None:
            reference, reference_energy = data, energy
        error = phase_error(reference, data)
        assert error < 1e-9 and abs(energy - reference_energy) < 1e-9
        report["layouts"].append(
            {
                "ranks": ranks,
                "threads": threads,
                "z2": 0,
                "eigenphase_error": error,
                "energy_error_Ha": abs(energy - reference_energy),
            }
        )
        print(report["layouts"][-1], flush=True)

    template = root / "neon-reference.inp"
    template.write_text(
        base.replace(
            "KPOINTS_SOURCE WILSON",
            "KPOINTS_SOURCE NNKP\n        NNKP_FILE loop.nnkp\n        WILSON_LOOP T",
        )
        .replace("Z2 T", "Z2 F")
        .replace("SEED_NAME neon", "SEED_NAME loop")
    )
    system = CP2KSystem(
        input_file=template,
        lattice=5 * np.eye(3),
        command=["mpiexec", "-n", "2", str(binary)],
        workdir=root / "z2pack-mpi",
        num_bands=8,
        env={
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "CP2K_DATA_DIR": str(ROOT / "data"),
        },
    )
    result = z2pack.surface.run(
        system=system,
        surface=lambda s, t: [t, s / 2, 0],
        num_lines=3,
        iterator=[8, 16, 32],
        pos_tol=1e-3,
    )
    assert z2pack.invariant.z2(result) == 0
    assert not result.convergence_report["line"]["PosCheck"]["FAILED"]
    for check in ("MoveCheck", "GapCheck"):
        assert not result.convergence_report["surface"][check]["FAILED"]
    report["mpi_z2pack_z2"] = 0
    report["mpi_z2pack_converged"] = True
    if args.stanene:
        directory = root / "stanene-2r-1t"
        run(binary, directory, (examples / "stanene-soc.inp").read_text(), 2, 1)
        compare_stanene(directory)
        report["stanene"] = json.loads((directory / "comparison.json").read_text())
    (root / "parallel.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
