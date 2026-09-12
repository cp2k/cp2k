#!/usr/bin/env python3
"""Check the geometry cache with actual CP2K types and interpolation routines.

Requires a completed GNU/OpenMP CMake build. No SCF calculation or Torch model
is used. All generated files and the machine-readable report stay in work-dir.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    build, work = args.build_dir.resolve(), args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    production = (root / "src/qs_vxc_atom.F").read_text()
    routines = []
    for name in (
        "create_native_grid_interpolation_stencil",
        "native_grid_lagrange_weights",
    ):
        matches = re.findall(
            rf"^   SUBROUTINE {name}\b.*?^   END SUBROUTINE {name}\b",
            production,
            re.M | re.S,
        )
        if len(matches) != 1:
            raise ValueError(f"Expected one production routine: {name}")
        routines.append(matches[0])
    kernel = work / "stencil_kernel.F90"
    kernel.write_text(
        "#define CPASSERT(cond) IF (.NOT. (cond)) ERROR STOP 'CPASSERT'\n"
        "MODULE stencil_kernel\n USE kinds, ONLY: dp\n USE cell_types, ONLY: cell_type\n"
        " USE pw_grid_types, ONLY: pw_grid_type\n USE qs_native_grid_cache\n"
        " IMPLICIT NONE\n CONTAINS\n" + "\n".join(routines) + "\nEND MODULE\n"
    )
    sources = [
        root / "src/qs_native_grid_cache.F",
        kernel,
        Path(__file__).with_name("native_grid_cache_fixture.F90"),
    ]
    commands = subprocess.check_output(
        [
            "ninja",
            "-C",
            str(build),
            "-t",
            "commands",
            "orbital_transformation_matrices_unittest",
        ],
        text=True,
    )
    link = shlex.split(commands.strip().splitlines()[-1])
    if link[:2] != [":", "&&"] or link[-2:] != ["&&", ":"]:
        raise ValueError("Unsupported CMake link layout")
    link = link[2:-2]
    objects = [x for x in link if x.endswith(".F.o")]
    if len(objects) != 1:
        raise ValueError("Expected one driver object")
    flags = [
        "-cpp",
        "-fopenmp",
        "-ffree-form",
        "-ffree-line-length-none",
        "-I",
        str(work),
        "-I",
        str(build / "src/mod_files"),
        "-J",
        str(work),
    ]
    flags += (
        ["-O1", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow"]
        if args.check
        else ["-O3"]
    )
    compiled = []
    for source in sources:
        target = work / (source.stem + ".o")
        subprocess.run(
            [link[0], *flags, "-c", str(source), "-o", str(target)],
            cwd=work,
            check=True,
        )
        compiled.append(str(target))
    executable = work / "cache_test"
    pos = link.index(objects[0])
    link[pos : pos + 1] = compiled
    link[link.index("-o") + 1] = str(executable)
    subprocess.run(link, cwd=build, check=True)
    result = subprocess.run(
        [str(executable)],
        cwd=work,
        text=True,
        capture_output=True,
        env=dict(os.environ, OMP_DYNAMIC="FALSE", OPENBLAS_NUM_THREADS="1"),
    )
    print(result.stdout, end="")
    print(result.stderr, end="")
    (work / "results.json").write_text(
        json.dumps(
            {
                "exit_code": result.returncode,
                "output": result.stdout,
                "stderr": result.stderr,
                "flags": flags,
                "sources": {
                    str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources
                },
                "scope": "Geometry-cache unit tests and extracted unchanged production stencil factory; no SCF/model timing",
            },
            indent=2,
        )
        + "\n"
    )
    result.check_returncode()


if __name__ == "__main__":
    main()
