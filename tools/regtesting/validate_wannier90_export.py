#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Run and validate CP2K-generated Wannier90 exports.

CP2K regtests can validate the generated ``.mmn`` and ``.eig`` files internally,
but the external Wannier90 executable has a few additional failure modes. This
helper is intentionally small: it can run ``wannier90.x`` for one seed, parse the
resulting ``.wout`` file, and compare final spreads/centres between two runs.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

WF_RE = re.compile(
    r"WF centre and spread\s+(\d+)\s+\(\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?),\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?),\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*\)\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
)
OMEGA_RE = re.compile(
    r"Omega\s+(I|D|OD|Total)\s*=\s*" r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
)
ERROR_RE = re.compile(r"(^|\s)(ERROR:|Error in|Exiting\.\.\.\.\.\.\.)")


@dataclass(frozen=True)
class WannierFunction:
    index: int
    center: tuple[float, float, float]
    spread: float


@dataclass(frozen=True)
class WoutSummary:
    seed: str
    path: Path
    functions: tuple[WannierFunction, ...]
    omega_i: float
    omega_d: float
    omega_od: float
    omega_total: float
    all_done: bool


def parse_wout(seed: str, workdir: Path) -> WoutSummary:
    path = workdir / f"{seed}.wout"
    text = path.read_text(encoding="utf-8", errors="replace")
    if ERROR_RE.search(text):
        raise ValueError(f"{path}: Wannier90 reported an error")
    final_pos = text.rfind("Final State")
    if final_pos < 0:
        raise ValueError(f"{path}: missing final Wannier90 state")
    final_text = text[final_pos:]

    functions: list[WannierFunction] = []
    for match in WF_RE.finditer(final_text):
        functions.append(
            WannierFunction(
                index=int(match.group(1)),
                center=(
                    float(match.group(2)),
                    float(match.group(3)),
                    float(match.group(4)),
                ),
                spread=float(match.group(5)),
            )
        )
    if not functions:
        raise ValueError(f"{path}: missing final Wannier-function centres")

    omegas: dict[str, float] = {}
    for match in OMEGA_RE.finditer(final_text):
        omegas[match.group(1)] = float(match.group(2))
    required = {"I", "D", "OD", "Total"}
    missing = required - omegas.keys()
    if missing:
        raise ValueError(
            f"{path}: missing omega field(s): {', '.join(sorted(missing))}"
        )

    return WoutSummary(
        seed=seed,
        path=path,
        functions=tuple(functions),
        omega_i=omegas["I"],
        omega_d=omegas["D"],
        omega_od=omegas["OD"],
        omega_total=omegas["Total"],
        all_done="All done: wannier90 exiting" in text,
    )


def run_wannier90(executable: str, seed: str, workdir: Path, timeout: float) -> None:
    result = subprocess.run(
        [executable, seed],
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        check=False,
    )
    (workdir / f"{seed}.wannier90.stdout").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(
            f"{executable} {seed} failed with exit code {result.returncode}"
        )


def print_summary(summary: WoutSummary) -> None:
    done = "yes" if summary.all_done else "no"
    print(f"seed = {summary.seed}")
    print(f"wout = {summary.path}")
    print(f"all_done = {done}")
    print(f"num_wann = {len(summary.functions)}")
    print(f"omega_i = {summary.omega_i:.12e}")
    print(f"omega_d = {summary.omega_d:.12e}")
    print(f"omega_od = {summary.omega_od:.12e}")
    print(f"omega_total = {summary.omega_total:.12e}")
    for function in summary.functions:
        cx, cy, cz = function.center
        print(
            f"wf {function.index:5d}: center = "
            f"({cx:.12e}, {cy:.12e}, {cz:.12e}), "
            f"spread = {function.spread:.12e}"
        )


def max_abs_delta(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right):
        raise ValueError(f"length mismatch: {len(left)} != {len(right)}")
    return max((abs(x - y) for x, y in zip(left, right)), default=0.0)


def compare_summaries(
    left: WoutSummary, right: WoutSummary
) -> tuple[float, float, float]:
    if len(left.functions) != len(right.functions):
        raise ValueError(
            f"number of Wannier functions differs: "
            f"{len(left.functions)} != {len(right.functions)}"
        )
    left_spreads = [function.spread for function in left.functions]
    right_spreads = [function.spread for function in right.functions]
    spread_delta = max_abs_delta(left_spreads, right_spreads)

    left_centers = [value for function in left.functions for value in function.center]
    right_centers = [value for function in right.functions for value in function.center]
    center_delta = max_abs_delta(left_centers, right_centers)
    omega_delta = abs(left.omega_total - right.omega_total)
    return omega_delta, spread_delta, center_delta


def command_run(args: argparse.Namespace) -> int:
    workdir = Path(args.workdir)
    if args.run:
        run_wannier90(args.executable, args.seed, workdir, args.timeout)
    summary = parse_wout(args.seed, workdir)
    print_summary(summary)
    if args.require_all_done and not summary.all_done:
        print("Wannier90 did not report normal completion.", file=sys.stderr)
        return 1
    if args.max_omega_total is not None and summary.omega_total > args.max_omega_total:
        print(
            f"Omega total {summary.omega_total:.12e} exceeds "
            f"{args.max_omega_total:.12e}.",
            file=sys.stderr,
        )
        return 1
    return 0


def command_compare(args: argparse.Namespace) -> int:
    left = parse_wout(args.left_seed, Path(args.left_workdir))
    right = parse_wout(args.right_seed, Path(args.right_workdir))
    omega_delta, spread_delta, center_delta = compare_summaries(left, right)
    print(f"omega_total_delta = {omega_delta:.12e}")
    print(f"spread_delta = {spread_delta:.12e}")
    print(f"center_delta = {center_delta:.12e}")

    failed = False
    if args.omega_total_tol is not None and omega_delta > args.omega_total_tol:
        failed = True
    if args.spread_tol is not None and spread_delta > args.spread_tol:
        failed = True
    if args.center_tol is not None and center_delta > args.center_tol:
        failed = True
    return 1 if failed else 0


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run", help="run and/or parse one Wannier90 seed"
    )
    run_parser.add_argument("seed", help="Wannier90 seed name")
    run_parser.add_argument(
        "--workdir", default=".", help="directory containing the seed files"
    )
    run_parser.add_argument(
        "--run", action="store_true", help="run wannier90.x before parsing"
    )
    run_parser.add_argument(
        "--executable", default="wannier90.x", help="Wannier90 executable"
    )
    run_parser.add_argument(
        "--timeout", type=float, default=300.0, help="execution timeout in seconds"
    )
    run_parser.add_argument("--max-omega-total", type=float, default=None)
    run_parser.add_argument("--require-all-done", action="store_true")
    run_parser.set_defaults(func=command_run)

    compare_parser = subparsers.add_parser(
        "compare", help="compare two parsed .wout files"
    )
    compare_parser.add_argument("left_seed")
    compare_parser.add_argument("right_seed")
    compare_parser.add_argument("--left-workdir", default=".")
    compare_parser.add_argument("--right-workdir", default=".")
    compare_parser.add_argument("--omega-total-tol", type=float, default=None)
    compare_parser.add_argument("--spread-tol", type=float, default=None)
    compare_parser.add_argument("--center-tol", type=float, default=None)
    compare_parser.set_defaults(func=command_compare)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (OSError, RuntimeError, subprocess.TimeoutExpired, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
