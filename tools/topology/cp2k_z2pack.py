"""CP2K overlap-loop adapter for Z2Pack (requires the NNKP export mode).

Each call uses an isolated directory, retains its log, and converges the same
SCF problem before evaluating the requested loop at that potential. A copied
SCF restart can be supplied in input_files; the SCF mesh must not be the loop.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from z2pack.system import OverlapSystem


def nnkp_text(kpoints, lattice):
    """Write one directed closed loop, preserving unwrapped reduced coordinates.

    lattice contains direct lattice vectors as ROWS, in Angstrom. The final
    point must equal the first plus an integer reciprocal lattice vector.
    """
    points = np.asarray(kpoints, dtype=float)
    cell = np.asarray(lattice, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 3:
        raise ValueError("A loop requires at least two points plus its endpoint")
    if not np.isfinite(points).all():
        raise ValueError("Nonfinite loop coordinates")
    if cell.shape != (3, 3) or not np.isfinite(cell).all():
        raise ValueError("lattice must be a finite 3x3 matrix")
    if abs(np.linalg.det(cell)) < 1e-12:
        raise ValueError("Singular lattice")
    delta = points[-1] - points[0]
    shift = np.rint(delta).astype(int)
    if not np.allclose(delta, shift, rtol=0, atol=1e-9):
        raise ValueError(
            "Loop endpoint must equal its start modulo a reciprocal vector"
        )
    if np.any(np.linalg.norm(np.diff(points, axis=0), axis=1) < 1e-12):
        raise ValueError("Adjacent loop points must differ")
    reciprocal = 2 * np.pi * np.linalg.inv(cell).T
    lines = ["# Explicit overlap loop for CP2K"]
    for name, matrix in (("real_lattice", cell), ("recip_lattice", reciprocal)):
        lines += [f"begin {name}"]
        lines += [" ".join(f"{x:.17g}" for x in row) for row in matrix]
        lines += [f"end {name}"]
    count = len(points) - 1
    lines += ["begin kpoints", str(count)]
    lines += [" ".join(f"{x:.17g}" for x in row) for row in points[:-1]]
    lines += ["end kpoints", "begin nnkpts", "1"]
    for i in range(1, count):
        lines.append(f"{i} {i + 1} 0 0 0")
    lines += [f"{count} 1 {' '.join(str(x) for x in shift)}", "end nnkpts"]
    lines += ["begin exclude_bands", "0", "end exclude_bands", ""]
    return "\n".join(lines)


def read_loop_mmn(path, kpoints, num_bands=None):
    """Read exactly the requested directed loop, validating dimensions and shifts."""
    lines = iter(Path(path).read_text().splitlines())
    try:
        next(lines)
        nb, nk, nn = map(int, next(lines).split())
        if nb < 1 or nk != len(kpoints) - 1 or nn != 1:
            raise ValueError("Unexpected overlap dimensions")
        if num_bands is not None and nb != num_bands:
            raise ValueError(f"Expected {num_bands} bands, received {nb}")
        closure = tuple(np.rint(np.asarray(kpoints[-1]) - kpoints[0]).astype(int))
        matrices = {}
        for _ in range(nk):
            source, target, gx, gy, gz = map(int, next(lines).split())
            if not 1 <= source <= nk or source in matrices:
                raise ValueError("Duplicate or invalid source point")
            expected_target = source % nk + 1
            expected_shift = closure if source == nk else (0, 0, 0)
            if target != expected_target or (gx, gy, gz) != expected_shift:
                raise ValueError("Overlap connection does not match the requested loop")
            values = []
            for _ in range(nb * nb):
                real, imag = map(float, next(lines).replace("D", "E").split())
                values.append(complex(real, imag))
            matrix = np.asarray(values).reshape(nb, nb).T
            if not np.isfinite(matrix).all():
                raise ValueError("Nonfinite overlap")
            matrices[source] = matrix
        if any(line.strip() for line in lines):
            raise ValueError("Unexpected trailing overlap blocks")
    except (StopIteration, TypeError) as error:
        raise ValueError("Truncated or malformed overlap file") from error
    return [matrices[i] for i in range(1, nk + 1)]


class CP2KSystem(OverlapSystem):
    """Z2Pack system with retained calculation directories and no shell execution.

    The CP2K input must use KPOINTS_SOURCE NNKP, NNKP_FILE loop.nnkp, and
    SEED_NAME loop. Relative dependencies (including restart files) belong in
    input_files. For MPI, command may be ['mpiexec', '-n', '2', '/path/cp2k.psmp'].
    """

    def __init__(
        self,
        *,
        input_file,
        lattice,
        command,
        workdir,
        input_files=(),
        num_bands=None,
        polar=True,
        singular_tol=1e-10,
        timeout=None,
        env=None,
    ):
        self.input_file = Path(input_file).resolve(strict=True)
        self.input_files = [Path(p).resolve(strict=True) for p in input_files]
        names = [p.name for p in self.input_files]
        if len(set(names)) != len(names) or set(names) & {
            "input.inp",
            "loop.nnkp",
            "loop.mmn",
        }:
            raise ValueError("Input dependencies have conflicting filenames")
        self.lattice = np.asarray(lattice, dtype=float)
        if isinstance(command, (str, bytes)) or not command:
            raise ValueError("command must be a nonempty argument list")
        self.command = [str(arg) for arg in command]
        self.workdir = Path(workdir).resolve()
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.num_bands = num_bands
        self.polar = polar
        self.singular_tol = singular_tol
        self.timeout = timeout
        self.env = dict(os.environ)
        self.env.setdefault("OMP_NUM_THREADS", "1")
        self.env.setdefault("OPENBLAS_NUM_THREADS", "1")
        self.env.update(env or {})
        self.last_run = None

    def get_mmn(self, kpt):
        text = nnkp_text(kpt, self.lattice)
        directory = Path(tempfile.mkdtemp(prefix="loop-", dir=self.workdir))
        self.last_run = directory
        shutil.copyfile(self.input_file, directory / "input.inp")
        for path in self.input_files:
            shutil.copyfile(path, directory / path.name)
        (directory / "loop.nnkp").write_text(text)
        argv = self.command + ["-i", "input.inp"]
        (directory / "request.json").write_text(
            json.dumps(
                {
                    "command": argv,
                    "points": np.asarray(kpt).tolist(),
                    "polar": self.polar,
                },
                indent=2,
            )
        )
        with (directory / "run.log").open("w") as log:
            result = subprocess.run(
                argv,
                cwd=directory,
                env=self.env,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=self.timeout,
                check=False,
            )
        if result.returncode:
            raise RuntimeError(
                f"CP2K failed ({result.returncode}); inspect {directory / 'run.log'}"
            )
        matrices = read_loop_mmn(directory / "loop.mmn", kpt, self.num_bands)
        singular_values = []
        for i, matrix in enumerate(matrices):
            u, s, vh = np.linalg.svd(matrix)
            singular_values.append(float(s.min()))
            if s.min() <= self.singular_tol:
                raise ValueError(
                    f"Singular subspace overlap on link {i + 1}; refine the loop or check band selection"
                )
            if self.polar:
                matrices[i] = u @ vh
        (directory / "diagnostics.json").write_text(
            json.dumps(
                {
                    "minimum_link_singular_value": min(singular_values),
                    "polar": self.polar,
                    "bands": len(matrices[0]),
                },
                indent=2,
            )
        )
        return matrices
