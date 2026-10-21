#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Literal, Sequence, Tuple, cast

import numpy as np
import numpy.typing as npt

FIT_FACTOR = 1.02
CERT_FACTOR = 1.05
ABS_TOL = 1.0e-15
FIT_RC_POINTS = 3
FIT_X_POINTS = 512
CERT_RC_POINTS = 25
CERT_X_POINTS = 4096
DENSE_INTERP_DEGREE = 7
SPARSE_INTERP_DEGREE = 4
DENSE_SOURCE_MIN = 20
GW_K = (26, 28, 30, 32, 34)
FloatArray = npt.NDArray[np.float64]


# ======================================================================================
def format_fortran_source(source: str) -> str:
    """Format generated Fortran with CP2K's precommit formatter."""
    repo_root = Path(__file__).resolve().parents[2]
    formatter = repo_root / "tools" / "precommit" / "format_fortran.py"
    result = subprocess.run(
        [sys.executable, str(formatter), "--stdout"],
        cwd=repo_root,
        input=source,
        capture_output=True,
        check=False,
        text=True,
    )
    if result.returncode or result.stderr:
        details = result.stderr.strip() or f"exit status {result.returncode}"
        raise RuntimeError(f"CP2K Fortran formatter failed:\n{details}")
    return result.stdout


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="Generates minimax_exp_data.F file.")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    input_path = script_dir / "1_xData"
    output_path = script_dir.parent.parent / "src" / "minimax" / "minimax_exp_data.F"
    assert input_path.exists()

    raw = parse_data_set(input_path)
    source = select_source_data(raw)
    by_k = group_by_k(source)
    sparse = fit_sparse_data(by_k)
    certify_sparse_data(by_k, sparse)
    code = format_fortran_source(generate_fortran_code(by_k, sparse))

    n_source = sum(len(v) for v in by_k.values())
    n_knots = sum(len(v) for v in sparse.values())
    n_coeff = sum(2 * k * len(v) for k, v in sparse.items())
    print(
        f"Sparse minimax data: {n_source} source points -> {n_knots} knots ({n_coeff} coefficient doubles)."
    )

    if args.check:
        assert output_path.exists(), f"Missing generated file {output_path}."
        assert output_path.read_text(encoding="utf8") == code
        print(f"File {output_path} is consistent with generator script.")
    else:
        output_path.write_text(code, encoding="utf8")
        print(f"Wrote {output_path}.")


# ======================================================================================
@dataclass
class Approximation:
    k: int
    Rc: float
    err: float
    omega: List[float]
    alpha: List[float]


ApproxByK = Dict[int, List[Approximation]]


# ======================================================================================
def parse_data_set(input_path: Path) -> List[Approximation]:
    with open(input_path / "tabelle", "r", encoding="utf8") as fh:
        lines = [line.strip() for line in fh.readlines()]

    ExpectEnum = Literal["HEADER", "SEPERATOR", "BODY"]
    expect: ExpectEnum = "HEADER"
    approx: List[Approximation] = []

    for line in lines[41:540]:
        if expect == "HEADER":
            if line == "":
                break
            assert line.startswith("k =  |")
            kvals = [int(x) for x in line[6:].split()]
            expect = "SEPERATOR"
            continue

        if expect == "SEPERATOR":
            assert line == "-" * 76
            expect = "BODY"
            continue

        if expect == "BODY":
            if line == "-" * 76:
                expect = "HEADER"
                continue

            Rc_str, errors_str = line.split("|")
            errors_str = errors_str.replace("--", "0.0")
            errors_str = errors_str.replace("**", "1e-18")
            Rc = float(Rc_str)
            errors = [float(x) for x in errors_str.split()]

            for k, err in zip(kvals, errors):
                if err == 0.0 or k >= 58:
                    continue
                Rc_filename = re.sub(r"E0?(?=\d)", "_", Rc_str.strip())
                filename = input_path / f"1_xk{k:02d}.{Rc_filename}"
                assert filename.exists()

                coeff_lines = filename.read_text(encoding="utf8").strip().split("\n")
                assert len(coeff_lines) == 2 * k
                ww = [float(x.split()[0]) for x in coeff_lines[:k]]
                aa = [float(x.split()[0]) for x in coeff_lines[k:]]
                approx.append(Approximation(k=k, Rc=Rc, err=err, omega=ww, alpha=aa))

    assert len(approx) == 2761
    return approx


# ======================================================================================
def format_array(
    array: Sequence[Any],
    indent: int,
    fmt: str = "{}",
    wrap: int = 5,
    begin: str = "[",
    end: str = "]",
    continuation_indent: int | None = None,
) -> List[str]:
    lines: List[str] = [" " * indent + begin]
    for i, x in enumerate(array):
        if i % wrap == 0 and i > 0:
            lines[-1] += "&"
            if continuation_indent is None:
                lines.append(" " * (indent + 1))
            else:
                lines.append(" " * continuation_indent + "& ")
        lines[-1] += fmt.format(x)
        if i + 1 != len(array):
            lines[-1] += ", "
    lines[-1] += end
    return lines


# ======================================================================================
def filter_orig(approx_in: List[Approximation]) -> List[Approximation]:
    """Filter to the ranges supported by the current exponential-minimax API."""

    Rc_range = [
        (2e00, 1e01),
        (2e00, 5e01),
        (2e00, 2e02),
        (2e00, 5e02),
        (2e00, 2e03),
        (2e00, 3e03),
        (2e00, 7e03),
        (1e01, 2e04),
        (1e01, 3e04),
        (1e01, 1e05),
        (1e01, 2e05),
        (1e01, 3e05),
        (1e01, 4e05),
        (1e01, 7e05),
        (1e01, 2e06),
        (1e01, 3e06),
        (2e01, 4e06),
        (4e01, 7e06),
        (5e01, 1e07),
        (6e01, 2e07),
        (9e01, 3e07),
        (1e02, 4e07),
        (2e02, 7e07),
        (2e02, 1e08),
        (2e02, 2e08),
        (4e02, 3e08),
        (1e04, 4e08),
        (1e08, 7e08),
        (7e08, 1e09),
        (7e08, 2e09),
        (7e08, 2e09),
        (7e08, 3e09),
        (7e08, 4e09),
        (7e08, 7e09),
        (7e08, 1e10),
        (1e10, 2e10),
        (1e10, 2e10),
        (1e10, 3e10),
        (1e10, 4e10),
        (1e10, 5e10),
        (1e10, 7e10),
        (1e10, 1e11),
        (1e11, 2e11),
        (1e11, 2e11),
        (1e11, 3e11),
        (1e11, 4e11),
        (1e11, 5e11),
        (1e11, 7e11),
        (1e11, 1e12),
        (2e08, 2e12),
        (1e12, 2e12),
        (1e12, 3e12),
        (4e11, 4e12),
    ]

    approx_out = []
    for a in approx_in:
        if a.k <= 53:
            Rc_min, Rc_max = Rc_range[a.k - 1]
            if Rc_min <= a.Rc <= Rc_max:
                approx_out.append(a)
    return approx_out


# ======================================================================================
def select_source_data(raw: List[Approximation]) -> List[Approximation]:
    """Select source ranges used by the unified exponential-minimax backend."""

    source = filter_orig(raw)
    # The historical GW grids use much wider low/intermediate-Rc ranges for
    # these five k values. Keep all published source points so the same
    # continuous interpolation backend can serve both the general and GW APIs.
    source = [a for a in source if a.k not in GW_K]
    source.extend(a for a in raw if a.k in GW_K)
    return source


# ======================================================================================
def group_by_k(approx: Iterable[Approximation]) -> ApproxByK:
    by_k: ApproxByK = {}
    for a in approx:
        by_k.setdefault(a.k, []).append(a)
    for values in by_k.values():
        values.sort(key=lambda a: a.Rc)
    assert sorted(by_k) == list(range(1, 54))
    return by_k


# ======================================================================================
def interpolation_degree(source: Sequence[Approximation]) -> int:
    """Choose a stable local interpolation degree for one k block."""

    k = source[0].k
    if k <= 15:
        return min(5, len(source) - 1)
    if len(source) >= DENSE_SOURCE_MIN:
        return DENSE_INTERP_DEGREE
    return min(SPARSE_INTERP_DEGREE, len(source) - 1)


# ======================================================================================
def lagrange_weights(nodes: FloatArray, x: float) -> FloatArray:
    weights: FloatArray = np.ones(len(nodes), dtype=np.float64)
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if i != j:
                weights[i] *= (x - nodes[j]) / (nodes[i] - nodes[j])
    return weights


# ======================================================================================
def sampled_error(Rc: float, coeff: FloatArray, npoints: int) -> float:
    k = len(coeff) // 2
    alpha = coeff[:k]
    omega = coeff[k:]
    x = np.exp(np.linspace(0.0, math.log(Rc), npoints))
    error = np.exp(-np.outer(x, alpha)) @ omega - 1.0 / x
    return float(np.max(np.abs(error)))


# ======================================================================================
def historical_choice(knots: Sequence[Approximation], Rc: float) -> FloatArray:
    """Preserve the current behavior for the few k values with only two source points."""

    assert len(knots) == 2
    coeff: List[FloatArray] = [cast(FloatArray, np.r_[a.alpha, a.omega]) for a in knots]
    errors = [sampled_error(Rc, c, 100) for c in coeff]
    selected: FloatArray = coeff[int(errors[1] <= errors[0])]
    return selected


# ======================================================================================
def interpolate_coefficients(
    knots: Sequence[Approximation], Rc: float, degree: int
) -> FloatArray:
    ranges: FloatArray = np.asarray([a.Rc for a in knots], dtype=np.float64)
    if Rc <= ranges[0]:
        return cast(FloatArray, np.r_[knots[0].alpha, knots[0].omega])
    if Rc >= ranges[-1]:
        return cast(FloatArray, np.r_[knots[-1].alpha, knots[-1].omega])
    if len(knots) == 2:
        return historical_choice(knots, Rc)

    x_nodes: FloatArray = np.log(ranges)
    x = math.log(Rc)
    i_right = int(np.searchsorted(x_nodes, x))
    n_full = degree + 1

    if len(knots) <= n_full:
        n_interp = len(knots)
        start = 0
    elif i_right == 1:
        n_interp = degree
        start = 0
    elif i_right == len(knots) - 1:
        n_interp = degree
        start = len(knots) - n_interp
    else:
        n_interp = n_full
        start = i_right - n_interp // 2

    start = max(0, min(start, len(knots) - n_interp))
    i_stop = start + n_interp
    weights = lagrange_weights(x_nodes[start:i_stop], x)
    values: FloatArray = cast(
        FloatArray,
        np.stack([np.r_[knots[i].alpha, knots[i].omega] for i in range(start, i_stop)]),
    )
    return weights @ values


# ======================================================================================
def coefficients_are_valid(coeff: FloatArray) -> bool:
    k = len(coeff) // 2
    alpha = coeff[:k]
    omega = coeff[k:]
    return bool(
        np.all(alpha > 0.0)
        and np.all(omega > 0.0)
        and np.all(np.diff(alpha) > 0.0)
        and np.all(np.diff(omega) > 0.0)
    )


# ======================================================================================
def make_cases(
    source: Sequence[Approximation], points_per_interval: int
) -> List[Tuple[float, float]]:
    # Keep source nodes first.  Besides being easier to inspect, this gives a
    # deterministic tie break when several trial interpolants violate ordering.
    cases: List[Tuple[float, float]] = [(a.Rc, a.err) for a in source]
    fractions = np.arange(1, points_per_interval + 1) / (points_per_interval + 1)
    for i in range(len(source) - 1):
        left = math.log(source[i].Rc)
        right = math.log(source[i + 1].Rc)
        for fraction in fractions:
            Rc = math.exp(left + fraction * (right - left))
            cases.append((Rc, source[i + 1].err))
    return cases


# ======================================================================================
def fit_knots(source: Sequence[Approximation]) -> List[Approximation]:
    """Adaptively add source knots until the minimax error envelope is satisfied."""

    degree = interpolation_degree(source)
    n = len(source)
    if n <= degree + 1:
        return list(source)

    n_initial = max(5, degree + 1)
    selected = sorted(
        set(np.rint(np.linspace(0, n - 1, n_initial)).astype(int).tolist())
    )
    cases = make_cases(source, FIT_RC_POINTS)

    while True:
        knots = [source[i] for i in selected]
        worst_ratio = -1.0
        worst_Rc = source[0].Rc
        for Rc, ref_error in cases:
            coeff = interpolate_coefficients(knots, Rc, degree)
            if coefficients_are_valid(coeff):
                error = sampled_error(Rc, coeff, FIT_X_POINTS)
                ratio = error / (FIT_FACTOR * ref_error + ABS_TOL)
            else:
                ratio = math.inf
            if ratio > worst_ratio:
                worst_ratio = ratio
                worst_Rc = Rc

        if worst_ratio <= 1.0:
            return knots

        candidates = [i for i in range(n) if i not in selected]
        if not candidates:
            raise RuntimeError(
                f"Could not fit k={source[0].k} within the minimax error envelope."
            )
        add = min(
            candidates,
            key=lambda i: abs(math.log(source[i].Rc) - math.log(worst_Rc)),
        )
        selected.append(add)
        selected.sort()


# ======================================================================================
def fit_sparse_data(by_k: ApproxByK) -> ApproxByK:
    sparse: ApproxByK = {}
    for k, source in sorted(by_k.items()):
        sparse[k] = fit_knots(source)
        print(
            f"k={k:2d}: {len(source):2d} -> {len(sparse[k]):2d} knots "
            f"(degree {interpolation_degree(source)})"
        )
    return sparse


# ======================================================================================
def certify_sparse_data(by_k: ApproxByK, sparse: ApproxByK) -> None:
    """Perform one dense certification after adaptive construction."""

    worst_ratio = -1.0
    worst_case: Tuple[int, float, float, float] | None = None
    for k, source in sorted(by_k.items()):
        knots = sparse[k]
        degree = interpolation_degree(source)
        for Rc, ref_error in make_cases(source, CERT_RC_POINTS):
            coeff = interpolate_coefficients(knots, Rc, degree)
            if not coefficients_are_valid(coeff):
                raise RuntimeError(
                    f"Invalid interpolated coefficients for k={k}, Rc={Rc:g}."
                )
            error = sampled_error(Rc, coeff, CERT_X_POINTS)
            ratio = error / (CERT_FACTOR * ref_error + ABS_TOL)
            if ratio > worst_ratio:
                worst_ratio = ratio
                worst_case = (k, Rc, error, ref_error)
            if ratio > 1.0:
                raise RuntimeError(
                    f"Certification failed for k={k}, Rc={Rc:g}: "
                    f"error={error:.16e}, reference={ref_error:.16e}."
                )

    assert worst_case is not None
    k, Rc, error, ref_error = worst_case
    print(
        f"Certification: max ratio={worst_ratio:.12f} "
        f"at k={k}, Rc={Rc:g} "
        f"(error={error:.6e}, reference={ref_error:.6e})."
    )


# ======================================================================================
def pointers(blocks: Sequence[Sequence[Any]]) -> List[int]:
    result = [1]
    for block in blocks:
        result.append(result[-1] + len(block))
    return result


# ======================================================================================
def generate_fortran_code(by_k: ApproxByK, sparse: ApproxByK) -> str:
    source_blocks = [by_k[k] for k in range(1, 54)]
    knot_blocks = [sparse[k] for k in range(1, 54)]
    k_p = pointers(knot_blocks)
    ref_p = pointers(source_blocks)
    knots = [a for block in knot_blocks for a in block]
    refs = [a for block in source_blocks for a in block]
    degrees = [interpolation_degree(by_k[k]) for k in range(1, 54)]

    output: List[str] = [
        "!--------------------------------------------------------------------------------------------------!",
        "!   CP2K: A general program to perform molecular dynamics simulations                              !",
        "!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !",
        "!                                                                                                  !",
        "!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !",
        "!--------------------------------------------------------------------------------------------------!",
        "",
        r"! **************************************************************************************************",
        r"!> \brief Sparse data for exponential minimax approximations of 1/x.",
        r"!>        Coefficients are interpolated locally as functions of log(Rc).",
        r"!>        The retained knots are selected adaptively against the published minimax-error envelope.",
        r"!>        Generated from http://www.mis.mpg.de/scicomp/EXP_SUM/1_x",
        r"!>        See also https://doi.org/10.1007/s00791-018-00308-4",
        r"!>        This file is generated by tools/minimax_tools/minimax_to_fortran_source.py.",
        r"! **************************************************************************************************",
        "MODULE minimax_exp_data",
        "   USE kinds,                           ONLY: dp",
        "",
        "   IMPLICIT NONE",
        "   PRIVATE",
        "   PUBLIC :: R_max, R_mm, err_ref, get_minimax_coeff_low, interp_degree, k_max, k_p, n_approx, &",
        "             n_ref, R_ref, ref_p",
        "",
        f"   INTEGER, PARAMETER :: k_max = 53",
        f"   INTEGER, PARAMETER :: n_approx = {len(knots)}",
        f"   INTEGER, PARAMETER :: n_ref = {len(refs)}",
        f"   REAL(KIND=dp), PARAMETER :: R_max = {max(a.Rc for a in refs):.1E}_dp",
        "",
        "   ! The generator uses a 2% construction margin and certifies against the existing 5% bound.",
        f"   ! Adaptive construction: {FIT_RC_POINTS} Rc probes/interval, {FIT_X_POINTS} x points.",
        f"   ! Dense certification: {CERT_RC_POINTS} Rc probes/interval, {CERT_X_POINTS} x points.",
        "",
        "   INTEGER, PARAMETER, DIMENSION(k_max + 1) :: k_p = &",
    ]
    output += format_array(k_p, indent=53, wrap=8)
    output += [
        "",
        "   INTEGER, PARAMETER, DIMENSION(k_max + 1) :: ref_p = &",
    ]
    output += format_array(ref_p, indent=55, wrap=8)
    output += [
        "",
        "   INTEGER, PARAMETER, DIMENSION(k_max) :: interp_degree = &",
    ]
    output += format_array(degrees, indent=58, wrap=12)
    output += [
        "",
        "   REAL(KIND=dp), PARAMETER, DIMENSION(n_approx) :: R_mm = &",
    ]
    output += format_array([a.Rc for a in knots], indent=52, fmt="{:.1E}_dp", wrap=4)
    output += [
        "",
        "   REAL(KIND=dp), PARAMETER, DIMENSION(n_ref) :: R_ref = &",
    ]
    output += format_array(
        [a.Rc for a in refs],
        indent=49,
        fmt="{:.1E}_dp",
        wrap=6,
        continuation_indent=3,
    )
    output += [
        "",
        "   REAL(KIND=dp), PARAMETER, DIMENSION(n_ref) :: err_ref = &",
    ]
    output += format_array(
        [a.err for a in refs],
        indent=51,
        fmt="{:.3E}_dp",
        wrap=6,
        continuation_indent=3,
    )
    output += [
        "",
        "CONTAINS",
        "",
        r"! **************************************************************************************************",
        r"!> \brief Return one retained exponential-minimax coefficient vector.",
        r"!> \param i Global retained-knot index.",
        r"!> \param aw Exponents followed by weights.",
        r"! **************************************************************************************************",
        "   PURE SUBROUTINE get_minimax_coeff_low(i, aw)",
        "      INTEGER, INTENT(IN)                                :: i",
        "      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)          :: aw",
        "",
        "      SELECT CASE (i)",
        "",
    ]

    for i, a in enumerate(knots, start=1):
        output += [f"      CASE ({i}) ! k={a.k}, Rc={a.Rc:g}", "         aw(:) = &"]
        output += format_array(
            [*a.alpha, *a.omega],
            indent=12,
            fmt="{:.17E}_dp",
            wrap=3,
        )
        output += [""]

    output += [
        "      END SELECT",
        "",
        "   END SUBROUTINE get_minimax_coeff_low",
        "END MODULE minimax_exp_data",
        "",
    ]
    return "\n".join(output)


# ======================================================================================
if __name__ == "__main__":
    main()

# EOF
