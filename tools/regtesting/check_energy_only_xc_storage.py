#!/usr/bin/env python3
"""Check energy-only XC buffer ownership using a completed GNU/Ninja CP2K build.

Extracts the actual atom-composite allocation and energy-only cleanup blocks.
The existing CP2K library supplies the real PW grid and pool implementations.
This is an ownership regression, not a substitute for a Skala SCF energy test.

Example (use the unpatched parent revision as BASE):
  python3 tools/regtesting/check_energy_only_xc_storage.py \
    --build-dir build --work-dir /tmp/xc-storage-check --baseline BASE
"""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess

FIXTURE = """
PROGRAM xc_storage_test
 USE libcp2k, ONLY: cp2k_init, cp2k_finalize
 USE pw_grid_types, ONLY: pw_grid_type
 USE pw_grids, ONLY: pw_grid_create, pw_grid_release
 USE pw_pool_types, ONLY: pw_pool_type, pw_pool_create, pw_pool_release
 USE pw_types, ONLY: pw_r3d_rs_type
 USE pw_methods, ONLY: pw_zero
 IMPLICIT NONE
 TYPE retained_buffers
   TYPE(pw_r3d_rs_type), POINTER :: rho(:) => NULL(), tau(:) => NULL()
 END TYPE retained_buffers
 TYPE(retained_buffers) :: retained(100)
 TYPE(pw_grid_type), POINTER :: grid => NULL()
 TYPE(pw_pool_type), POINTER :: auxbas_pw_pool => NULL(), xc_pw_pool => NULL()
 TYPE(pw_r3d_rs_type), POINTER :: my_vxc_rho(:) => NULL(), my_vxc_tau(:) => NULL()
 TYPE(pw_r3d_rs_type), POINTER :: v_rspace_new(:) => NULL(), v_tau_rspace(:) => NULL()
 INTEGER :: bounds(2, 3), mspin, ispin, iteration, cache, variant, expected
 LOGICAL, PARAMETER :: fixed = {fixed}

 CALL cp2k_init()
 bounds(1, :) = 0
 bounds(2, :) = 7
 DO cache = 0, 8, 8
   DO mspin = 1, 2
     DO variant = 0, 3
       CALL pw_grid_create(grid, bounds)
       CALL pw_pool_create(auxbas_pw_pool, grid, max_cache=cache)
       xc_pw_pool => auxbas_pw_pool
       DO iteration = 1, 100
         NULLIFY(my_vxc_rho, my_vxc_tau, v_rspace_new, v_tau_rspace)
         IF (variant /= 0) THEN
{allocation}
           v_rspace_new => my_vxc_rho
           v_tau_rspace => my_vxc_tau
           NULLIFY(my_vxc_rho, my_vxc_tau)
           ! Also exercise callers that return only one kind of XC potential.
           IF (variant == 1) THEN
             DO ispin = 1, mspin
               CALL auxbas_pw_pool%give_back_pw(v_tau_rspace(ispin))
             END DO
             DEALLOCATE(v_tau_rspace)
           ELSE IF (variant == 2) THEN
             DO ispin = 1, mspin
               CALL auxbas_pw_pool%give_back_pw(v_rspace_new(ispin))
             END DO
             DEALLOCATE(v_rspace_new)
           END IF
         END IF
{cleanup}
         expected = 2
         IF (.NOT. fixed .AND. variant /= 0) THEN
           expected = expected + iteration*mspin
           IF (variant == 3) expected = expected + iteration*mspin
         END IF
         IF (grid%ref_count /= expected) ERROR STOP 'Unexpected grid reference count'
         IF (fixed) THEN
           IF (ASSOCIATED(v_rspace_new)) ERROR STOP 'Density potential still associated'
           IF (ASSOCIATED(v_tau_rspace)) ERROR STOP 'Tau potential still associated'
         ELSE
           retained(iteration)%rho => v_rspace_new
           retained(iteration)%tau => v_tau_rspace
         END IF
       END DO
       PRINT *, 'cache/spins/variant/references:', cache, mspin, variant, grid%ref_count
       ! Recover deliberately leaked baseline objects after measuring their growth,
       ! so CP2K can finalize its grid-owned MPI communicators without an abort.
       IF (.NOT. fixed) THEN
         DO iteration = 1, 100
           IF (ASSOCIATED(retained(iteration)%rho)) THEN
             DO ispin = 1, mspin
               CALL auxbas_pw_pool%give_back_pw(retained(iteration)%rho(ispin))
             END DO
             DEALLOCATE(retained(iteration)%rho)
           END IF
           IF (ASSOCIATED(retained(iteration)%tau)) THEN
             DO ispin = 1, mspin
               CALL auxbas_pw_pool%give_back_pw(retained(iteration)%tau(ispin))
             END DO
             DEALLOCATE(retained(iteration)%tau)
           END IF
         END DO
       END IF
       IF (grid%ref_count /= 2) ERROR STOP 'Fixture teardown leaked a grid reference'
       CALL pw_pool_release(auxbas_pw_pool)
       NULLIFY(xc_pw_pool)
       CALL pw_grid_release(grid)
     END DO
   END DO
 END DO
 CALL cp2k_finalize()
END PROGRAM xc_storage_test
"""


def cleanup_block(source):
    end = source.index("      END IF ! .NOT. just energy")
    start = source.rindex("\n      ELSE\n", 0, end)
    end = source.index("         IF (do_hfx) THEN", start, end)
    return source[start + len("\n      ELSE\n") : end]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--baseline", required=True, help="Unfixed Git revision")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    build, work = args.build_dir.resolve(), args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    vxc = (root / "src/qs_vxc.F").read_text()
    start = vxc.index("            ALLOCATE (my_vxc_rho(mspin), my_vxc_tau(mspin))")
    end = vxc.index("         ELSE IF (native_skala_grid) THEN", start)
    allocation = vxc[start:end]
    original = subprocess.check_output(
        ["git", "show", f"{args.baseline}:src/qs_ks_methods.F"], cwd=root, text=True
    )
    patched = (root / "src/qs_ks_methods.F").read_text()
    if cleanup_block(original).strip():
        raise ValueError("The baseline is not the expected unfixed energy-only path")
    if "give_back_pw" not in cleanup_block(patched):
        raise ValueError("No energy-only PW cleanup found in the patched source")
    commands = subprocess.check_output(
        ["ninja", "-C", str(build), "-t", "commands", "kpsym_unittest"], text=True
    )
    link = shlex.split(commands.strip().splitlines()[-1])
    if link[:2] != [":", "&&"] or link[-2:] != ["&&", ":"]:
        raise ValueError("Unsupported CMake link command layout")
    link = link[2:-2]
    link = [arg for arg in link if not arg.startswith("-Wl,--dependency-file=")]
    objects = [x for x in link if x.endswith(".F.o")]
    if len(objects) != 1:
        raise ValueError("Expected a single unit-test driver object")
    compiler = link[0]
    flags = [
        "-fopenmp",
        "-ffree-line-length-none",
        "-fcheck=all",
        "-O1",
        "-g",
        "-I",
        str(build / "src/mod_files"),
    ]
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    # Assert exact CP2K ownership counts instead of library-wide LSan reports.
    env["LSAN_OPTIONS"] = "detect_leaks=0"
    results = {}
    for name, source in [("original", original), ("fixed", patched)]:
        src = work / f"{name}.F90"
        src.write_text(
            FIXTURE.format(
                fixed=".TRUE." if name == "fixed" else ".FALSE.",
                allocation=allocation,
                cleanup=cleanup_block(source),
            )
        )
        obj, exe = work / f"{name}.o", work / name
        subprocess.run(
            [compiler, *flags, "-c", str(src), "-o", str(obj)], check=True, cwd=work
        )
        command = list(link)
        command[command.index(objects[0])] = str(obj)
        command[command.index("-o") + 1] = str(exe)
        subprocess.run(command, check=True, cwd=build)
        run = subprocess.run(
            [str(exe)], env=env, capture_output=True, text=True, timeout=60
        )
        results[name] = {
            "returncode": run.returncode,
            "checks": [
                line.strip()
                for line in run.stdout.splitlines()
                if "cache/spins/variant" in line
            ],
            "cache_disabled_warnings": run.stdout.count("hit max_cache"),
            "stderr": run.stderr,
        }
        if run.returncode:
            results[name]["output_tail"] = run.stdout[-4000:]
            print(json.dumps(results, indent=2))
            raise SystemExit(run.returncode)
    (work / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
