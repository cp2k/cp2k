#!/usr/bin/env python3
"""Compare the GAPW interpolation/adjoint kernels with an unmodified Git revision.

Uses a completed GNU/OpenMP CMake CP2K build (no Torch model or SCF calculation).
Extracts the actual private routines and caller loops, not rewritten numerical
surrogates. Generated sources, executable and JSON results stay in --work-dir.
The shared CP2K library supplies the real grids, harmonics and orbital tables.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

USES = """
 USE kinds, ONLY: dp
 USE qs_grid_atom, ONLY: grid_atom_type
 USE qs_harmonics_atom, ONLY: harmonics_atom_type
 USE orbital_pointers, ONLY: indco, indso, nco, ncoset, nsoset
 USE orbital_transformation_matrices, ONLY: orbtramat
 USE spherical_harmonics, ONLY: y_lm
 USE cell_types, ONLY: cell_type
 USE particle_types, ONLY: particle_type
 IMPLICIT NONE
 CONTAINS
"""

WRAPPER = """
 SUBROUTINE project(cell, particle_set, grid_atom, harmonics, nspins, flags, &
   composite_grid_coords, composite_grid_atom, composite_local_atoms, &
   composite_atom_start, composite_atom_end, composite_density_grad, &
   composite_grad_grad, composite_kin_grad, cross_cutoff, &
   vxc_h, vxc_s, vxg_h, vxg_s, vtau_h, vtau_s)
 TYPE(cell_type), POINTER :: cell
 TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
 TYPE(grid_atom_type), POINTER :: grid_atom
 TYPE(harmonics_atom_type), POINTER :: harmonics
 INTEGER, INTENT(IN) :: nspins
 LOGICAL, INTENT(IN) :: flags(3)
 REAL(dp), INTENT(IN) :: composite_grid_coords(:, :), composite_density_grad(:, :)
 REAL(dp), INTENT(IN) :: composite_grad_grad(:, :, :), composite_kin_grad(:, :)
 INTEGER, INTENT(IN) :: composite_grid_atom(:), composite_local_atoms(:)
 INTEGER, INTENT(IN) :: composite_atom_start(:), composite_atom_end(:)
 REAL(dp), INTENT(IN) :: cross_cutoff
 REAL(dp), POINTER :: vxc_h(:, :, :), vxc_s(:, :, :), vtau_h(:, :, :), vtau_s(:, :, :)
 REAL(dp), POINTER :: vxg_h(:, :, :, :), vxg_s(:, :, :, :)
 REAL(dp), ALLOCATABLE :: vxc_h_local(:, :, :), vxc_s_local(:, :, :), vtau_h_local(:, :, :), vtau_s_local(:, :, :)
 REAL(dp), ALLOCATABLE :: vxg_h_local(:, :, :, :), vxg_s_local(:, :, :, :)
 INTEGER :: base_shift(3), composite_row, idir, image_i1, image_i2, image_i3
 INTEGER :: image_lower(3), image_upper(3)
 INTEGER :: image_shift(3), image_shell(3), jdir, target_atom, iatom
 INTEGER :: composite_local_atom, composite_local_natom, composite_nflat
 REAL(dp) :: cross_density_adjoint(2), cross_grad_adjoint(3, 2), cross_kin_adjoint(2)
 REAL(dp) :: cross_displacement(3), fractional(3), image_translation(3)
 LOGICAL :: lsd, use_atom_composite_density, use_atom_composite_gradient, use_atom_composite_tau
 iatom = 1
 composite_nflat = SIZE(composite_grid_atom)
 composite_local_natom = SIZE(composite_local_atoms)
 lsd = nspins == 2
 use_atom_composite_density = flags(1)
 use_atom_composite_gradient = flags(2)
 use_atom_composite_tau = flags(3)
 image_shell = 0
 DO idir = 1, 3
   IF (cell%perd(idir) == 1) image_shell(idir) = &
     CEILING(cross_cutoff*SQRT(SUM(cell%h_inv(idir, :)**2))) + 1
 END DO
{loop}
 END SUBROUTINE project
"""


def routine(source, name):
    matches = re.findall(
        rf"^   SUBROUTINE {name}\b.*?^   END SUBROUTINE {name}\b",
        source,
        re.M | re.S,
    )
    if len(matches) != 1:
        raise ValueError(f"Expected one definition of {name}, got {len(matches)}")
    return matches[0]


def caller_loop(source, parallel):
    start = source.index(
        "cross_cutoff = gapw_atom_grid_support_radius(",
        source.index("SUBROUTINE calculate_vxc_atom("),
    )
    # The first occurrence is in the forward path; the adjoint is the last one
    # before the target-owner MPI reduction.
    end = source.index("! Model rows are distributed by target atom")
    region = source[start:end]
    start = region.rindex("cross_cutoff = gapw_atom_grid_support_radius(")
    region = region[start:]
    if parallel is None:
        parallel = "!$OMP PARALLEL" in region
    if parallel:
        start = region.index("!$OMP PARALLEL")
        end = region.index("!$OMP END PARALLEL") + len("!$OMP END PARALLEL")
    else:
        start = region.index("DO composite_local_atom = 1, composite_local_natom")
        end = region.rindex("END IF")
    return region[start:end].strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument(
        "--baseline",
        required=True,
        help="Unmodified Git revision containing the serial or parallel adjoint",
    )
    parser.add_argument(
        "--check", action="store_true", help="Bounds/FPE checks, no timing"
    )
    parser.add_argument("--rows", type=int, default=4000)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    build, work = args.build_dir.resolve(), args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    baseline = subprocess.check_output(
        ["git", "show", f"{args.baseline}:src/qs_vxc_atom.F"], cwd=root, text=True
    )
    patched = (root / "src/qs_vxc_atom.F").read_text()
    helpers = [
        "radial_node_derivative_coefficients",
        "atom_grid_interpolation_weights",
        "interpolate_gapw_atom_grid_fields",
        "add_gapw_atom_grid_interpolation_adjoint",
    ]
    old_loop, new_loop = caller_loop(baseline, None), caller_loop(patched, True)
    sources = []
    for name, source, loop in [
        ("original", baseline, old_loop),
        ("values_only", patched, old_loop),
        ("screened", patched, re.sub(r"^!\$OMP.*\n?", "", new_loop, flags=re.M)),
        ("parallel", patched, new_loop),
    ]:
        path = work / f"{name}.F90"
        content = (
            "#define CPASSERT(cond) IF (.NOT. (cond)) ERROR STOP 'CPASSERT'\n"
            f"MODULE {name}_kernel\n"
            + USES
            + "\n".join(routine(source, n) for n in helpers)
            + "\n"
            + routine(patched, "atom_grid_image_bounds")
            + WRAPPER.format(loop=loop)
            + f"\nEND MODULE {name}_kernel\n"
        )
        path.write_text(content)
        sources.append(path)
    # Reuse the build's resolved compiler and link dependencies without guessing
    # site-specific library paths. Invoke argv directly, never through a shell.
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
        raise ValueError("Unsupported CMake link command layout")
    link = link[2:-2]
    objects = [x for x in link if x.endswith(".F.o")]
    if len(objects) != 1:
        raise ValueError("Expected one unit-test driver object")
    compiler = link[0]
    executable = work / ("adjoint_checked" if args.check else "adjoint_release")
    flags = [
        "-cpp",
        "-fopenmp",
        "-ffree-line-length-none",
        "-I",
        str(build / "src/mod_files"),
    ]
    flags += (
        ["-O1", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow"]
        if args.check
        else ["-O3", "-funroll-loops"]
    )
    sources.append(Path(__file__).with_name("gapw_atom_adjoint_fixture.F90"))
    compiled = []
    for source in sources:
        target = work / (source.stem + ".o")
        subprocess.run(
            [compiler, *flags, "-c", str(source), "-o", str(target)],
            cwd=work,
            check=True,
        )
        compiled.append(str(target))
    pos = link.index(objects[0])
    link[pos : pos + 1] = compiled
    link[link.index("-o") + 1] = str(executable)
    subprocess.run(link, cwd=build, check=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    run_command = [
        str(executable),
        "check" if args.check else "benchmark",
        str(args.rows),
    ]
    if Path("/usr/bin/time").exists():
        run_command = [
            "/usr/bin/time",
            "-l" if sys.platform == "darwin" else "-v",
            *run_command,
        ]
    completed = subprocess.run(
        run_command, cwd=work, env=env, text=True, capture_output=True
    )
    print(completed.stdout, end="")
    print(completed.stderr, end="")
    report = {
        "baseline_revision": subprocess.check_output(
            ["git", "rev-parse", args.baseline], cwd=root, text=True
        ).strip(),
        "baseline_sha256": hashlib.sha256(baseline.encode()).hexdigest(),
        "patched_sha256": hashlib.sha256(patched.encode()).hexdigest(),
        "compiler": subprocess.check_output(
            [compiler, "--version"], text=True
        ).splitlines()[0],
        "flags": flags,
        "exit_code": completed.returncode,
        "output": completed.stdout,
        "stderr": completed.stderr,
        "generated_sources_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources
        },
        "scope": "Extracted forward/adjoint kernels and exact caller loop, real CP2K types/harmonics; no SCF/model test",
    }
    (
        work / ("checked-results.json" if args.check else "benchmark-results.json")
    ).write_text(json.dumps(report, indent=2) + "\n")
    completed.check_returncode()


if __name__ == "__main__":
    main()
