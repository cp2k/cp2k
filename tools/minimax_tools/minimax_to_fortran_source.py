#!/usr/bin/python3

import re
import argparse
from pathlib import Path
from typing import Any, List, Literal
from dataclasses import dataclass


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="Generates minimax_exp_k53.F file.")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    input_path = script_dir / "1_xData"
    output_path = script_dir.parent.parent / "src" / "minimax" / "minimax_exp_k53.F"
    assert input_path.exists()
    assert output_path.exists()

    approx = parse_data_set(input_path)
    approx = filter_orig(approx)
    code = generate_fortran_code(approx)

    if args.check:
        assert output_path.read_text(encoding="utf8") == code
        print(f"File {output_path} is consisted with generator script.")
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
                break  # we are done
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
            errors_str = errors_str.replace("**", "1e-18")  # an educated guess
            Rc = float(Rc_str)
            errors = [float(x) for x in errors_str.split()]

            for k, err in zip(kvals, errors):
                if err == 0.0 or k >= 58:
                    continue
                # slight change of notation e.g. 6E09 --> 6_9
                Rc_filename = re.sub("E0?(?=\d)", "_", Rc_str.strip())
                filename = input_path / f"1_xk{k:02d}.{Rc_filename}"
                assert filename.exists()

                # read coefficients
                coeff_lines = filename.read_text(encoding="utf8").strip().split("\n")
                assert len(coeff_lines) == 2 * k
                ww = [float(x.split()[0]) for x in coeff_lines[:k]]
                aa = [float(x.split()[0]) for x in coeff_lines[k:]]
                approx.append(Approximation(k=k, Rc=Rc, err=err, omega=ww, alpha=aa))

    assert len(approx) == 2761  # sanity check
    return approx


# ======================================================================================
def format_array(
    array: List[Any],
    indent: int,
    fmt: str = "{}",
    wrap: int = 5,
    begin: str = "[",
    end: str = "]",
) -> List[str]:

    lines: List[str] = [" " * indent + begin]
    for i, x in enumerate(array):
        if i % wrap == 0 and i > 0:
            lines[-1] += "&"
            lines.append(" " * (indent + 1))
        lines[-1] += fmt.format(x)
        if not i + 1 == len(array):
            lines[-1] += ", "
    lines[-1] += end

    return lines


# ======================================================================================
def filter_orig(approx_in: List[Approximation]) -> List[Approximation]:
    """Filters to roughly match the first version of the data set."""

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
            if Rc_min <= a.Rc and a.Rc <= Rc_max:
                approx_out.append(a)

    return approx_out


# ======================================================================================
def generate_fortran_code(approx: List[Approximation]) -> str:
    approx.sort(key=lambda a: (a.k, a.Rc))

    # pointers to k
    k_p: List[int] = []
    my_k = 0
    for i, a in enumerate(approx):
        if my_k < a.k:
            my_k = a.k
            k_p.append(i + 1)
    k_p.append(len(approx) + 1)

    output: List[str] = []
    output += [
        "!--------------------------------------------------------------------------------------------------!",
        "!   CP2K: A general program to perform molecular dynamics simulations                              !",
        "!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !",
        "!                                                                                                  !",
        "!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !",
        "!--------------------------------------------------------------------------------------------------!",
        "",
        r"! **************************************************************************************************",
        r"!> \brief Routines to calculate the minimax coefficients in order to",
        r"!>        approximate 1/x as a sum over exponential functions",
        r"!>        1/x ~ SUM_{i}^{K} w_i EXP(-a_i * x) for x belonging to [1:Rc].",
        r"!>        This module contains coefficients for minimax approximations with 1 <= k <= 53.",
        r"!>        Generated from data from",
        r"!>        http://www.mis.mpg.de/scicomp/EXP_SUM/1_x",
        r"!>        See also https://doi.org/10.1007/s00791-018-00308-4",
        r"!>        This module should not be modified manually and should not be used anywhere",
        r"!>        except in main minimax module.",
        r"!>        This file was created using the scripts in cp2k/tools/minimax_tools.",
        r"! **************************************************************************************************",
        "MODULE minimax_exp_k53",
        "   USE kinds,                           ONLY: dp",
        "",
        "   IMPLICIT NONE",
        "   PRIVATE",
        "   PUBLIC :: R_max, R_mm, err_mm, get_minimax_coeff_low, k_max, k_mm, k_p, n_approx",
        "",
        "   INTEGER, PARAMETER :: n_approx = {}".format(len(approx)),
        "   INTEGER, PARAMETER :: n_k = {}".format(len(k_p) - 1),
        "   INTEGER, PARAMETER :: k_max = {}".format((max(a.k for a in approx))),
        "   REAL(KIND=dp), PARAMETER :: R_max = {:.1E}_dp".format(
            max(a.Rc for a in approx)
        ),
        "",
    ]

    output += ["   INTEGER, PARAMETER, DIMENSION(n_k + 1) :: k_p = &"]
    output += format_array(k_p, indent=45)
    output += [""]

    output += ["   INTEGER, PARAMETER, DIMENSION(n_approx) :: k_mm = &"]
    output += format_array([a.k for a in approx], indent=46)
    output += [""]

    output += ["   REAL(KIND=dp), PARAMETER, DIMENSION(n_approx) :: R_mm = &"]
    output += format_array([a.Rc for a in approx], indent=52, fmt="{:.1E}_dp")
    output += [""]

    output += ["   REAL(KIND=dp), PARAMETER, DIMENSION(n_approx) :: err_mm = &"]
    output += format_array([a.err for a in approx], indent=52, fmt="{:.3E}_dp")
    output += [""]

    output += [
        "CONTAINS",
        "",
        r"! **************************************************************************************************",
        r"!> \brief ...",
        r"!> \param i ...",
        r"!> \param aw ...",
        r"! **************************************************************************************************",
        "   SUBROUTINE get_minimax_coeff_low(i, aw)",
        "      INTEGER, INTENT(IN)                                :: i",
        "      REAL(KIND=dp), DIMENSION(k_mm(i)*2), INTENT(OUT)   :: aw",
        "",
        "      SELECT CASE (i)",
        "",
    ]

    for i, a in enumerate(approx):
        output += [f"      CASE ({i+1})"]
        output += [""]
        output += ["         aw(:) = & ! a"]
        output += format_array(
            a.alpha, indent=12, fmt="{}_dp", wrap=3, begin="[", end=", &"
        )
        output += ["             ! w"]
        output += format_array(
            a.omega, indent=12, fmt="{}_dp", wrap=3, begin=" ", end="]"
        )
        output += [""]

    output += [
        "      END SELECT",
        "",
        "   END SUBROUTINE get_minimax_coeff_low",
        "END MODULE minimax_exp_k53",
        "",
    ]

    return "\n".join(output)


# ======================================================================================
if __name__ == "__main__":
    main()

# EOF
