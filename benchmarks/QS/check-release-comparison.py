#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run the benchmark in this directory with multiple combinations of OpenMP threads and MPI ranks,
calculate the average of the final total energies, the variance and checks that each energy is
within 1e-10 of the average.
"""

import subprocess
import re
import pathlib
import argparse
import tempfile
import os
import sys

SCRIPT_DIR = pathlib.Path(__file__).parent

TEST_COMBINATIONS = [
    # OMP_NUM_THREADS, MPI_RANKS
    (0, 8),  # 0 switches to using [ps]opt
    (1, 8),
    (2, 4),
    (4, 2),
    (8, 1),
    (8, 0),  # 0 switches to using s[smp,opt]
]


def cp2k_version(omp_num_threads, mpi_ranks):
    """Generates the correct CP2K version string based on number of threads and ranks"""
    return f"{'p' if mpi_ranks > 0 else 's'}{'smp' if omp_num_threads > 0 else 'opt'}"


def ofname(omp_num_threads, mpi_ranks, input_file, version):
    """Returns a string of the form 'H2O-32-sopt-1-8.out' with additional steps"""
    return (
        f"{input_file.stem}-{version}-{max(1,mpi_ranks)}-{max(1,omp_num_threads)}.out"
    )


TOTAL_FORCE_EVAL_RE = re.compile(
    r"^\s*ENERGY\| Total FORCE_EVAL.+?:\s*(?P<energy>.+)\n", re.MULTILINE
)


def run_benchmark(cp2k_exe_dir, input_file, mpi_wrapper="mpiexec -np {mpi_ranks}"):
    energies = []

    for omp_num_threads, mpi_ranks in TEST_COMBINATIONS:
        env = {**os.environ, "OMP_NUM_THREADS": "1"}  # override OMP_NUM_THREADS
        args = []

        if mpi_ranks > 0:
            args = mpi_wrapper.format(mpi_ranks=mpi_ranks).split()

        if omp_num_threads > 0:
            env["OMP_NUM_THREADS"] = str(omp_num_threads)

        version = cp2k_version(omp_num_threads, mpi_ranks)
        outfile = SCRIPT_DIR / ofname(omp_num_threads, mpi_ranks, input_file, version)

        args += [cp2k_exe_dir / f"cp2k.{version}", input_file]
        with tempfile.TemporaryDirectory(
            prefix="release-comparison."
        ) as tmpdir, outfile.open("w+", encoding="utf-8") as fhandle:
            print(
                f"Running benchmark '{input_file.name}' with {args[-2]} using {mpi_ranks} MPI ranks and {omp_num_threads} OpenMP threads"
            )
            subprocess.run(
                args,
                stdout=fhandle,
                stderr=subprocess.STDOUT,
                check=True,
                encoding="utf-8",
                cwd=tmpdir,
                env=env,
            )

            fhandle.seek(0)  # be kind, rewind
            for match in TOTAL_FORCE_EVAL_RE.finditer(fhandle.read()):
                pass  # get the last match (it's an MD, we have multiple total energies)
            energies += [float(match["energy"])]
            print(f".. run was successful, final total energy: {energies[-1]:20.14f}")

    mean = sum(energies) / len(energies)
    if len(energies) > 1:
        var = sum((e - mean) ** 2 for e in energies) / len(energies)
    else:
        var = float("NaN")

    allok = True

    print(
        """
CP2K | # Threads | # Ranks | Energy
-----+-----------+---------+------------------------"""
    )

    for (omp_num_threads, mpi_ranks), energy in zip(TEST_COMBINATIONS, energies):
        version = cp2k_version(omp_num_threads, mpi_ranks)
        mark = "✔"
        if abs(energy - mean) >= 1e-10:
            mark = "!"
            allok = False

        print(
            f"{version} | {omp_num_threads:9d} | {mpi_ranks:7d} | {energy:20.14f} {mark}"
        )

    print(
        f"""\
----------------------------------------------------
mean: {mean:20.14f} +/- {var:16.12e}"""
    )

    return allok


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--cp2k-exe-dir",
        metavar="<CP2K-exe-directory>",
        type=pathlib.Path,
        help="Path to where cp2k.* can be found",
        default=(SCRIPT_DIR / "../../../exe/Linux-x86-64-gfortran"),
    )
    parser.add_argument(
        "--input",
        metavar="<CP2K-input to use>",
        type=pathlib.Path,
        help="Path to a CP2K input file to run the benchmark for",
        default=(SCRIPT_DIR / "H2O-32.inp"),
    )
    args = parser.parse_args()

    # resolve any relative path absolute in respect to the current working directory
    cp2k_exe_dir = (pathlib.Path.cwd() / args.cp2k_exe_dir).resolve()
    input_file = (pathlib.Path.cwd() / args.input).resolve()

    if not cp2k_exe_dir.exists():
        print(f"ERROR: '{cp2k_exe_dir}' could not be found")
        sys.exit(1)

    if not run_benchmark(cp2k_exe_dir, input_file):
        sys.exit(1)
