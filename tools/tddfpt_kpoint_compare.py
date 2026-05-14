#!/usr/bin/env python3
"""Compare small TDDFPT k-point calculations against Gamma supercells.

This is a research/debugging helper for the TDDFPT k-point kernel branch.  It
generates a tiny H2 chain input, runs KERNEL NONE/FULL with a k-point mesh and
with the corresponding Gamma supercell, then prints the first TDDFPT states.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

TDDFPT_LINE = re.compile(
    r"^\s*TDDFPT\|\s+"
    r"(?P<state>\d+)\s+"
    r"(?P<energy>[-+0-9.EeDd]+)\s+"
    r"(?P<dx>[-+0-9.EeDd]+)\s+"
    r"(?P<dy>[-+0-9.EeDd]+)\s+"
    r"(?P<dz>[-+0-9.EeDd]+)\s+"
    r"(?P<osc>[-+0-9.EeDd]+)"
)
CHECKSUM_LINE = re.compile(r"TDDFPT\s*:\s*CheckSum\s*=\s*(?P<value>[-+0-9.EeDd]+)")
METRIC_LINE = re.compile(r"Metric symmetry error:\s*(?P<value>[-+0-9.EeDd]+)")


@dataclass(frozen=True)
class State:
    index: int
    energy_ev: float
    dipole: tuple[float, float, float]
    oscillator: float


@dataclass(frozen=True)
class Result:
    label: str
    output: Path
    states: tuple[State, ...]
    checksum: float | None
    metric_error: float | None


def as_float(text: str) -> float:
    return float(text.replace("D", "E"))


def parse_output(label: str, output: Path) -> Result:
    states: list[State] = []
    checksum: float | None = None
    metric_error: float | None = None

    for line in output.read_text(errors="replace").splitlines():
        match = TDDFPT_LINE.match(line)
        if match:
            states.append(
                State(
                    index=int(match.group("state")),
                    energy_ev=as_float(match.group("energy")),
                    dipole=(
                        as_float(match.group("dx")),
                        as_float(match.group("dy")),
                        as_float(match.group("dz")),
                    ),
                    oscillator=as_float(match.group("osc")),
                )
            )
            continue

        match = CHECKSUM_LINE.search(line)
        if match:
            checksum = as_float(match.group("value"))
            continue

        match = METRIC_LINE.search(line)
        if match:
            metric_error = as_float(match.group("value"))

    if not states:
        raise RuntimeError(f"No TDDFPT states found in {output}")

    return Result(
        label=label,
        output=output,
        states=tuple(states),
        checksum=checksum,
        metric_error=metric_error,
    )


def xc_block(functional: str) -> str:
    functional = functional.upper()
    if functional == "PBE":
        return "      &XC_FUNCTIONAL PBE\n      &END XC_FUNCTIONAL"
    if functional in {"PADE", "LDA"}:
        return "      &XC_FUNCTIONAL PADE\n      &END XC_FUNCTIONAL"
    raise ValueError(f"Unsupported functional for this helper: {functional}")


def h2_coords(nx: int, ny: int, nz: int) -> str:
    lines: list[str] = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                x0 = 4.0 * ix
                y0 = 4.0 * iy
                z0 = 4.0 * iz
                lines.append(
                    f"      H {x0 + 1.65:8.4f} {y0 + 2.00:8.4f} {z0 + 2.00:8.4f}"
                )
                lines.append(
                    f"      H {x0 + 2.35:8.4f} {y0 + 2.00:8.4f} {z0 + 2.00:8.4f}"
                )
    return "\n".join(lines)


def make_input(
    *,
    project: str,
    kernel: str,
    functional: str,
    nstates: int,
    grid: tuple[int, int, int],
    gamma_centered: bool,
    supercell: bool,
    triplet: bool,
) -> str:
    nx, ny, nz = grid if supercell else (1, 1, 1)
    kpoints = ""
    if not supercell:
        gamma_line = "      GAMMA_CENTERED T\n" if gamma_centered else ""
        kpoints = f"""    &KPOINTS
      FULL_GRID ON
{gamma_line}      PARALLEL_GROUP_SIZE -1
      SCHEME MONKHORST-PACK {grid[0]} {grid[1]} {grid[2]}
      WAVEFUNCTIONS COMPLEX
    &END KPOINTS
"""

    triplet_line = "      RKS_TRIPLETS T\n" if triplet else ""
    return f"""&GLOBAL
  PRINT_LEVEL LOW
  PROJECT {project}
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
{kpoints}    &MGRID
      CUTOFF 100
      REL_CUTOFF 30
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-10
      METHOD GPW
    &END QS
    &SCF
      ADDED_MOS 4
      CHOLESKY OFF
      EPS_EIGVAL 1.0E-8
      EPS_SCF 1.0E-8
      MAX_SCF 50
      SCF_GUESS ATOMIC
      &MIXING
        ALPHA 0.35
        METHOD BROYDEN_MIXING
      &END MIXING
      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
    &XC
{xc_block(functional)}
    &END XC
  &END DFT
  &PROPERTIES
    &TDDFPT
      KERNEL {kernel}
      NSTATES {nstates}
{triplet_line}    &END TDDFPT
  &END PROPERTIES
  &SUBSYS
    &CELL
      ABC {4.0 * nx:.1f} {4.0 * ny:.1f} {4.0 * nz:.1f}
    &END CELL
    &COORD
{h2_coords(nx, ny, nz)}
    &END COORD
    &KIND H
      BASIS_SET ORB DZVP-MOLOPT-GGA-GTH-q1
      POTENTIAL GTH-GGA-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
"""


def run_case(
    *,
    cp2k: Path,
    data_dir: Path,
    work_dir: Path,
    label: str,
    input_text: str,
    keep_going: bool,
) -> Result | None:
    inp = work_dir / f"{label}.inp"
    out = work_dir / f"{label}.out"
    inp.write_text(input_text)
    env = os.environ.copy()
    env["CP2K_DATA_DIR"] = str(data_dir)
    env.setdefault("OMP_NUM_THREADS", "1")
    print(f"Running {label} ...", flush=True)
    completed = subprocess.run(
        [str(cp2k), "-i", str(inp), "-o", str(out)],
        cwd=work_dir,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if completed.returncode != 0:
        print(completed.stdout, file=sys.stderr)
        print(f"FAILED {label}; output: {out}", file=sys.stderr)
        if keep_going:
            return None
        raise SystemExit(completed.returncode)
    return parse_output(label, out)


def print_result(result: Result, limit: int, bright_threshold: float) -> None:
    metric = "n/a" if result.metric_error is None else f"{result.metric_error:.3e}"
    checksum = "n/a" if result.checksum is None else f"{result.checksum:.6e}"
    print(f"\n{result.label}")
    print(f"  output:   {result.output}")
    print(f"  checksum: {checksum}")
    print(f"  metric:   {metric}")
    print("  state       energy/eV       oscillator        |dipole|")
    for state in result.states[:limit]:
        dip_norm = sum(x * x for x in state.dipole) ** 0.5
        print(
            f"  {state.index:5d}  {state.energy_ev:14.6f}"
            f"  {state.oscillator:15.6e}  {dip_norm:14.6e}"
        )
    bright = [state for state in result.states if state.oscillator > bright_threshold]
    if bright:
        print(f"  bright states, f > {bright_threshold:.1e}:")
        for state in bright[:limit]:
            print(
                f"    {state.index:5d}  {state.energy_ev:14.6f}  {state.oscillator:15.6e}"
            )
    else:
        print(f"  bright states, f > {bright_threshold:.1e}: none")


def compare_pair(kpoint: Result, supercell: Result, limit: int) -> None:
    print(f"\nPairwise energy deltas: {kpoint.label} - {supercell.label}")
    print("  state       dk/eV      osc(k)        osc(sc)")
    for left, right in zip(kpoint.states[:limit], supercell.states[:limit]):
        print(
            f"  {left.index:5d}  {left.energy_ev - right.energy_ev:10.6f}"
            f"  {left.oscillator:12.5e}  {right.oscillator:12.5e}"
        )

    print(f"\nNearest-energy matches: {kpoint.label} -> {supercell.label}")
    print("  kp      k/eV   sc      sc/eV       dk/eV      osc(k)        osc(sc)")
    for left in kpoint.states[:limit]:
        right = min(
            supercell.states, key=lambda state: abs(left.energy_ev - state.energy_ev)
        )
        print(
            f"  {left.index:2d}  {left.energy_ev:10.6f}"
            f"  {right.index:3d}  {right.energy_ev:10.6f}"
            f"  {left.energy_ev - right.energy_ev:10.6f}"
            f"  {left.oscillator:12.5e}  {right.oscillator:12.5e}"
        )


def parse_grid(values: Iterable[str]) -> tuple[int, int, int]:
    grid = tuple(int(v) for v in values)
    if len(grid) != 3 or any(v < 1 for v in grid):
        raise argparse.ArgumentTypeError("grid needs three positive integers")
    return grid  # type: ignore[return-value]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cp2k", type=Path, default=Path("build-serial/bin/cp2k.ssmp"))
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument(
        "--grid", nargs=3, default=("2", "1", "1"), metavar=("NX", "NY", "NZ")
    )
    parser.add_argument("--functional", default="PBE", choices=("PBE", "PADE", "LDA"))
    parser.add_argument("--nstates", type=int, default=4)
    parser.add_argument(
        "--kernels", nargs="+", default=("NONE", "FULL"), choices=("NONE", "FULL")
    )
    parser.add_argument("--no-gamma-centered", action="store_true")
    parser.add_argument("--triplet", action="store_true")
    parser.add_argument("--bright-threshold", type=float, default=1.0e-8)
    parser.add_argument("--keep-going", action="store_true")
    args = parser.parse_args()

    repo = Path.cwd()
    cp2k = (repo / args.cp2k).resolve() if not args.cp2k.is_absolute() else args.cp2k
    data_dir = (
        (repo / args.data_dir).resolve()
        if not args.data_dir.is_absolute()
        else args.data_dir
    )
    grid = parse_grid(args.grid)
    work_dir = args.work_dir or Path(tempfile.mkdtemp(prefix="cp2k_tddfpt_kp_compare_"))
    work_dir = work_dir.resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    if not cp2k.exists():
        raise SystemExit(f"CP2K binary not found: {cp2k}")
    if not data_dir.exists():
        raise SystemExit(f"CP2K data directory not found: {data_dir}")

    print(f"Work dir: {work_dir}")
    print(f"Grid:     {grid[0]} {grid[1]} {grid[2]}")
    print(f"XC:       {args.functional}")
    print(f"Triplet:  {args.triplet}")

    results: dict[str, Result] = {}
    for kernel in args.kernels:
        for supercell in (False, True):
            label = (
                f"H2_{'gamma_supercell' if supercell else 'kpoint'}_{kernel.lower()}"
            )
            input_text = make_input(
                project=label,
                kernel=kernel,
                functional=args.functional,
                nstates=args.nstates,
                grid=grid,
                gamma_centered=not args.no_gamma_centered,
                supercell=supercell,
                triplet=args.triplet,
            )
            result = run_case(
                cp2k=cp2k,
                data_dir=data_dir,
                work_dir=work_dir,
                label=label,
                input_text=input_text,
                keep_going=args.keep_going,
            )
            if result is not None:
                results[label] = result

    for result in results.values():
        print_result(result, args.nstates, args.bright_threshold)

    for kernel in args.kernels:
        left = results.get(f"H2_kpoint_{kernel.lower()}")
        right = results.get(f"H2_gamma_supercell_{kernel.lower()}")
        if left is not None and right is not None:
            compare_pair(left, right, args.nstates)

    print(f"\nKept inputs/outputs in {work_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
