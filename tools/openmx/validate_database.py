#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Validate every generated entry, optionally with CP2K's real basis reader.

CP2K RUN_TYPE NONE and ghost kinds test parsing/initialization without choosing
physical pseudopotentials or claiming an OpenMX-equivalent Hamiltonian.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np

from convert_basis import ELEMENTS, gaussian_overlap


def read_basis(path):
    lines = iter(
        line.split("#", 1)[0].strip() for line in path.read_text().splitlines()
    )
    lines = iter(line for line in lines if line)
    entries = {}
    for line in lines:
        element, name = line.split()
        if element not in ELEMENTS[1:] or name in entries:
            raise ValueError("Unknown element or duplicate basis name")
        nsets = int(next(lines))
        if nsets < 1:
            raise ValueError("Empty basis")
        nao, primitives, norm_error = 0, 0, 0.0
        for _ in range(nsets):
            header = [int(value) for value in next(lines).split()]
            if len(header) != 5:
                raise ValueError("Expected one angular channel per set")
            principal, lmin, lmax, nprimitive, nradial = header
            if (
                lmin != lmax
                or not 0 <= lmin <= 4
                or principal < 1
                or min(nprimitive, nradial) < 1
            ):
                raise ValueError("Invalid set header")
            table = np.array(
                [
                    [float(value) for value in next(lines).split()]
                    for _ in range(nprimitive)
                ]
            )
            if (
                table.shape != (nprimitive, 1 + nradial)
                or not np.isfinite(table).all()
                or np.any(table[:, 0] <= 0)
            ):
                raise ValueError("Invalid exponent/coefficient table")
            alpha, coeff = table[:, 0], table[:, 1:]
            gram = coeff.T @ gaussian_overlap(alpha, lmin) @ coeff
            if np.linalg.eigvalsh((gram + gram.T) / 2)[0] <= 0:
                raise ValueError(f"{name}: linearly dependent/invalid contractions")
            norm_error = max(norm_error, float(np.max(np.abs(np.diag(gram) - 1))))
            nao += (2 * lmin + 1) * nradial
            primitives += nprimitive
        if norm_error > 1e-6:
            raise ValueError(
                f"{name}: unstable contraction normalization ({norm_error})"
            )
        entries[name] = {
            "element": element,
            "nao": nao,
            "nsets": nsets,
            "primitive_rows": primitives,
            "normalization_max_error": norm_error,
        }
    if not entries:
        raise ValueError("Empty basis file")
    return entries


def initialization_input(basis, entries):
    lines = [
        "&GLOBAL",
        "  PROJECT basis_read_test",
        "  RUN_TYPE NONE",
        "  PRINT_LEVEL MEDIUM",
        "&END GLOBAL",
        "&FORCE_EVAL",
        "  METHOD Quickstep",
        "  &DFT",
        f"    BASIS_SET_FILE_NAME {basis}",
        "    &QS",
        "      EPS_DEFAULT 1.0E-10",
        "    &END QS",
        "    &MGRID",
        "      CUTOFF 30",
        "      REL_CUTOFF 10",
        "      NGRIDS 1",
        "    &END MGRID",
        "    &SCF",
        "      MAX_SCF 1",
        "    &END SCF",
        "    &XC",
        "      &XC_FUNCTIONAL PBE",
        "      &END XC_FUNCTIONAL",
        "    &END XC",
        "  &END DFT",
        "  &SUBSYS",
        "    &CELL",
        "      ABC 8 8 8",
        "    &END CELL",
        "    &COORD",
    ]
    # Coordinates merely make each KIND active in initialization. Keep ghosts
    # distinct so CP2K's topology sanity checks also succeed.
    for i, _ in enumerate(entries):
        lines.append(f"      basis_{i} {1+2*(i%4)} {1+2*((i//4)%4)} {1+2*(i//16)}")
    lines.extend(["    &END COORD"])
    for i, (name, metadata) in enumerate(entries):
        lines.extend(
            [
                f"    &KIND basis_{i}",
                f"      ELEMENT {metadata['element']}",
                f"      BASIS_SET {name}",
                "      GHOST T",
                "    &END KIND",
            ]
        )
    lines.extend(["  &END SUBSYS", "&END FORCE_EVAL"])
    return "\n".join(lines) + "\n"


def cp2k_batch(task):
    executable, basis, entries = task
    with tempfile.TemporaryDirectory(prefix="basis-reader-test-") as temporary:
        root = Path(temporary)
        (root / "read.inp").write_text(initialization_input(basis, entries))
        run = subprocess.run(
            [executable, "-i", "read.inp", "-o", "read.out"],
            cwd=root,
            env=dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1"),
            capture_output=True,
            text=True,
            timeout=120,
        )
        output = (root / "read.out").read_text() if (root / "read.out").exists() else ""
        # Verify actual basis loading, not just a successful input syntax check.
        kinds = len(re.findall(r"Atomic kind:", output))
        ok = run.returncode == 0 and "PROGRAM ENDED" in output and kinds == len(entries)
        return {
            "names": [name for name, _ in entries],
            "passed": ok,
            "basis_kinds_printed": kinds,
            "returncode": run.returncode,
            **({} if ok else {"diagnostic": (run.stderr + output)[-6000:]}),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("basis", type=Path)
    parser.add_argument("--cp2k", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if (
        args.output.exists()
        or not 1 <= args.batch_size <= 32
        or not 1 <= args.workers <= 8
    ):
        parser.error("Need a new output path, batch size 1..32 and workers 1..8")
    entries = read_basis(args.basis)
    report_path = args.basis.with_name(args.basis.name + ".json")
    report = json.loads(report_path.read_text())
    sha = hashlib.sha256(args.basis.read_bytes()).hexdigest()
    if sha != report["basis_sha256"]:
        raise ValueError("Basis/report hash mismatch")
    expected = {
        b["cp2k_name"]: b for b in report["basis_sets"] if b.get("exported", True)
    }
    if entries.keys() != expected.keys() or any(
        v["nao"] != expected[k]["nao"] or v["element"] != expected[k]["element"]
        for k, v in entries.items()
    ):
        raise ValueError("Basis/report coverage or AO-count mismatch")
    result = {
        "schema": "openmx-cp2k-database-validation-v1",
        "basis_sha256": sha,
        "parsed_entries": len(entries),
        "elements": len({e["element"] for e in entries.values()}),
        "max_normalization_error": max(
            e["normalization_max_error"] for e in entries.values()
        ),
        "cp2k_basis_reader_tested": bool(args.cp2k),
        "batches": [],
    }
    if args.cp2k:
        values = list(entries.items())
        tasks = [
            (
                str(args.cp2k.resolve()),
                str(args.basis.resolve()),
                values[i : i + args.batch_size],
            )
            for i in range(0, len(values), args.batch_size)
        ]
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for i, batch in enumerate(pool.map(cp2k_batch, tasks)):
                result["batches"].append(batch)
                print(
                    f"CP2K batch {i+1}/{len(tasks)}: {'PASS' if batch['passed'] else 'FAIL'}",
                    flush=True,
                )
    result["passed"] = all(b["passed"] for b in result["batches"])
    with args.output.open("x") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k != "batches"}, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
