#!/usr/bin/env python3
"""Generate CP2K BASIS_TB and gxTB_BASIS data files.

The generated file contains static Gaussian contractions converted from
q-vSZP, PTB vDZP/vDZP auxiliary, and tblite GFN1/GFN2 Slater exponents.
The gxTB_BASIS file also stores the q/CN-dependent coefficient column used by
the experimental BASIS_SET ORB TB reader.
"""

from __future__ import annotations

import datetime as _dt
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_FILE = ROOT / "data" / "BASIS_TB"
GXTB_OUT_FILE = ROOT / "data" / "gxTB_BASIS"
CP2K_STO_NG = ROOT / "src" / "aobasis" / "sto_ng.F"

REPOS = {
    "qvSZP": "https://github.com/grimme-lab/qvSZP",
    "ptb": "https://github.com/grimme-lab/ptb",
    "save_tblite": "https://github.com/thfroitzheim/save_tblite",
}

SYMBOLS = [
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]

L_MAP = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5, "i": 6}
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")


@dataclass
class Shell:
    l: int
    exps: list[float]
    coeffs: list[float]
    n_label: int | None = None
    coeffs_env: list[float] | None = None


def run(cmd: list[str], cwd: Path | None = None) -> str:
    return subprocess.check_output(cmd, cwd=cwd, text=True).strip()


def clone_sources(tmp: Path) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for name, url in REPOS.items():
        dest = tmp / name
        subprocess.check_call(
            ["git", "clone", "--depth", "1", url, str(dest)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        paths[name] = dest
    return paths


def repo_ref(path: Path) -> str:
    return run(["git", "rev-parse", "--short=12", "HEAD"], cwd=path)


def parse_float(value: str) -> float:
    value = value.strip().replace("_dp", "").replace("_wp", "")
    return float(value.replace("D", "E").replace("d", "e"))


def fortran_numbers(text: str, as_int: bool = False) -> list[int] | list[float]:
    cleaned: list[str] = []
    for line in text.splitlines():
        cleaned.append(line.split("!")[0])
    body = "\n".join(cleaned).replace("&", " ")
    body = body.replace("_wp", "").replace("_dp", "")
    values = FLOAT_RE.findall(body)
    if as_int:
        return [int(float(v.replace("D", "E").replace("d", "e"))) for v in values]
    return [parse_float(v) for v in values]


def extract_fortran_array(
    text: str, name: str, as_int: bool = False
) -> list[int] | list[float]:
    pat = re.compile(
        r"::\s*" + re.escape(name) + r"\b[^=]*=\s*(?:reshape\(\s*)?\[&?(.*?)\]",
        re.IGNORECASE | re.DOTALL,
    )
    match = pat.search(text)
    if not match:
        raise ValueError(f"Could not find Fortran array {name}")
    return fortran_numbers(match.group(1), as_int=as_int)


def extract_fortran_int_parameter(text: str, name: str) -> int:
    match = re.search(
        r"integer\s*,\s*parameter\s*::\s*" + re.escape(name) + r"\s*=\s*(\d+)",
        text,
        re.IGNORECASE,
    )
    if not match:
        raise ValueError(f"Could not find integer parameter {name}")
    return int(match.group(1))


def parse_cp2k_sto_ng(
    path: Path,
) -> dict[tuple[int, int], tuple[list[float], list[float]]]:
    """Parse CP2K's STO-NG table as keyed by (m, ng).

    CP2K computes m = nq*(nq-1)/2 + lq + 1.
    """
    table: dict[tuple[int, int], tuple[list[float], list[float]]] = {}
    outer_m: int | None = None
    ng: int | None = None
    alphas: dict[int, float] = {}
    coeffs: dict[int, float] = {}

    outer_case = re.compile(r"^      CASE \((\d+)\)\s*!\s*\d+[A-Za-z]")
    inner_case = re.compile(r"^         CASE \((\d+)\)")
    assignment = re.compile(
        r"alpha\((\d+)\)\s*=\s*([^;]+);\s*coef\(\1\)\s*=\s*([^!\n]+)"
    )

    for line in path.read_text().splitlines():
        outer = outer_case.match(line)
        if outer:
            outer_m = int(outer.group(1))
            ng = None
            alphas = {}
            coeffs = {}
            continue

        inner = inner_case.match(line)
        if inner and outer_m is not None:
            ng = int(inner.group(1))
            alphas = {}
            coeffs = {}
            continue

        if outer_m is None or ng is None:
            continue

        match = assignment.search(line)
        if not match:
            continue
        idx = int(match.group(1))
        alphas[idx] = parse_float(match.group(2))
        coeffs[idx] = parse_float(match.group(3))
        if len(alphas) == ng and len(coeffs) == ng:
            table[(outer_m, ng)] = (
                [alphas[i] for i in range(1, ng + 1)],
                [coeffs[i] for i in range(1, ng + 1)],
            )

    return table


def sto_ng_shell(
    sto_table: dict[tuple[int, int], tuple[list[float], list[float]]],
    principal_n: int,
    l: int,
    zeta: float,
    ng: int,
) -> Shell:
    m = principal_n * (principal_n - 1) // 2 + l + 1
    try:
        alphas, coeffs = sto_table[(m, ng)]
    except KeyError as exc:
        raise ValueError(
            f"Missing STO-{ng}G table for n={principal_n}, l={l}, m={m}"
        ) from exc
    scale = zeta * zeta
    return Shell(l=l, exps=[a * scale for a in alphas], coeffs=coeffs, n_label=l + 1)


def parse_star_basis(path: Path) -> dict[int, list[Shell]]:
    """Parse grimme-lab qvSZP/PTB star-separated basis files.

    q-vSZP primitive lines carry an extra environment coefficient; CP2K's
    static basis format cannot express it, so this parser keeps only the
    exponent and the base contraction coefficient.
    """
    records: list[list[str]] = []
    current: list[str] = []

    for line in path.read_text().splitlines():
        if line.strip().startswith("*"):
            if current:
                records.append(current)
                current = []
            continue
        current.append(line)
    if current:
        records.append(current)

    parsed: dict[int, list[Shell]] = {}
    shell_line = re.compile(r"^\s*(\d+)\s+([A-Za-z])\s*$")

    for record in records:
        lines = [
            ln.strip() for ln in record if ln.strip() and not ln.strip().startswith("#")
        ]
        if not lines:
            continue
        first = lines[0].split()
        if not first or not first[0].isdigit():
            continue
        z = int(first[0])
        shells: list[Shell] = []
        i = 1
        while i < len(lines):
            match = shell_line.match(lines[i])
            if not match:
                i += 1
                continue
            nprim = int(match.group(1))
            l_char = match.group(2).lower()
            if l_char not in L_MAP:
                raise ValueError(f"Unsupported angular momentum {l_char!r} in {path}")
            l = L_MAP[l_char]
            exps: list[float] = []
            coeffs: list[float] = []
            for prim in lines[i + 1 : i + 1 + nprim]:
                nums = fortran_numbers(prim)
                if len(nums) < 2:
                    raise ValueError(f"Bad primitive line in {path}: {prim!r}")
                exps.append(float(nums[0]))
                coeffs.append(float(nums[1]))
            if len(exps) != nprim:
                raise ValueError(f"Short shell for Z={z} in {path}")
            shells.append(Shell(l=l, exps=exps, coeffs=coeffs, n_label=l + 1))
            i += nprim + 1
        if shells:
            parsed[z] = shells

    return parsed


def parse_qvszp_tb_basis(
    path: Path,
) -> dict[int, tuple[tuple[float, float, float, float], list[Shell]]]:
    """Parse q-vSZP-style basis files including environment coefficients.

    The record header contains k2, k1, and k3. k0 is one for the current
    upstream parameterization.
    """
    records: list[list[str]] = []
    current: list[str] = []

    for line in path.read_text().splitlines():
        if line.strip().startswith("*"):
            if current:
                records.append(current)
                current = []
            continue
        current.append(line)
    if current:
        records.append(current)

    parsed: dict[int, tuple[tuple[float, float, float, float], list[Shell]]] = {}
    shell_line = re.compile(r"^\s*(\d+)\s+([A-Za-z])\s*$")

    for record in records:
        lines = [
            ln.strip() for ln in record if ln.strip() and not ln.strip().startswith("#")
        ]
        if not lines:
            continue
        first = lines[0].split()
        if len(first) < 4 or not first[0].isdigit():
            continue
        z = int(first[0])
        k2 = float(first[1])
        k1 = float(first[2])
        k3 = float(first[3])
        shells: list[Shell] = []
        i = 1
        while i < len(lines):
            match = shell_line.match(lines[i])
            if not match:
                i += 1
                continue
            nprim = int(match.group(1))
            l_char = match.group(2).lower()
            if l_char not in L_MAP:
                raise ValueError(f"Unsupported angular momentum {l_char!r} in {path}")
            l = L_MAP[l_char]
            exps: list[float] = []
            coeffs: list[float] = []
            coeffs_env: list[float] = []
            for prim in lines[i + 1 : i + 1 + nprim]:
                nums = fortran_numbers(prim)
                if len(nums) < 3:
                    raise ValueError(f"Bad TB primitive line in {path}: {prim!r}")
                exps.append(float(nums[0]))
                coeffs.append(float(nums[1]))
                coeffs_env.append(float(nums[2]))
            if len(exps) != nprim:
                raise ValueError(f"Short shell for Z={z} in {path}")
            shells.append(
                Shell(
                    l=l,
                    exps=exps,
                    coeffs=coeffs,
                    coeffs_env=coeffs_env,
                    n_label=l + 1,
                )
            )
            i += nprim + 1
        if shells:
            parsed[z] = ((1.0, k1, k2, k3), shells)

    return parsed


def parse_tblite_gfn(
    path: Path, sto_table: dict[tuple[int, int], tuple[list[float], list[float]]]
) -> dict[int, list[Shell]]:
    text = path.read_text()
    max_elem = extract_fortran_int_parameter(text, "max_elem")
    max_shell = extract_fortran_int_parameter(text, "max_shell")
    nshell = extract_fortran_array(text, "nshell", as_int=True)
    angular = extract_fortran_array(text, "ang_shell", as_int=True)
    principal = extract_fortran_array(text, "principal_quantum_number", as_int=True)
    zetas = extract_fortran_array(text, "slater_exponent", as_int=False)

    if len(nshell) != max_elem:
        raise ValueError(f"Unexpected nshell length in {path}")
    expected_2d = max_elem * max_shell
    if (
        len(angular) != expected_2d
        or len(principal) != expected_2d
        or len(zetas) != expected_2d
    ):
        raise ValueError(f"Unexpected 2D array length in {path}")

    parsed: dict[int, list[Shell]] = {}
    for z in range(1, max_elem + 1):
        shells: list[Shell] = []
        for ish in range(nshell[z - 1]):
            idx = (z - 1) * max_shell + ish
            l = int(angular[idx])
            n = int(principal[idx])
            zeta = float(zetas[idx])
            shells.append(sto_ng_shell(sto_table, n, l, zeta, 6))
        if shells:
            parsed[z] = shells
    return parsed


def write_basis_block(
    lines: list[str], z: int, names: list[str], shells: list[Shell]
) -> None:
    symbol = SYMBOLS[z]
    lines.append(f"{symbol} {' '.join(names)}")
    lines.append(f"  {len(shells)}")
    for shell in shells:
        nprim = len(shell.exps)
        n_label = shell.n_label if shell.n_label is not None else shell.l + 1
        lines.append(f"  {n_label:2d} {shell.l:d} {shell.l:d} {nprim:3d} 1")
        for exp, coeff in zip(shell.exps, shell.coeffs):
            lines.append(f"    {exp:24.14f} {coeff:24.14f}")
    lines.append("#")


def append_family(
    lines: list[str],
    title: str,
    names: list[str],
    basis: dict[int, list[Shell]],
) -> None:
    lines.append("")
    lines.append(f"# ---- {title} ----")
    for z in sorted(basis):
        if z <= 0 or z >= len(SYMBOLS):
            raise ValueError(f"Unsupported atomic number {z}")
        write_basis_block(lines, z, names, basis[z])


def write_tb_basis_block(
    lines: list[str],
    z: int,
    names: list[str],
    params: tuple[float, float, float, float],
    shells: list[Shell],
) -> None:
    symbol = SYMBOLS[z]
    k0, k1, k2, k3 = params
    fixed_q = 0.0
    fixed_cn = 0.0
    lines.append(f"{symbol} {' '.join(names)}")
    lines.append(
        f"  {k0:20.12f} {k1:20.12f} {k2:20.12f} {k3:20.12f}"
        f" {fixed_q:20.12f} {fixed_cn:20.12f}"
    )
    lines.append(f"  {len(shells)}")
    for shell in shells:
        if shell.coeffs_env is None:
            raise ValueError(f"TB shell for Z={z} is missing environment coefficients")
        nprim = len(shell.exps)
        n_label = shell.n_label if shell.n_label is not None else shell.l + 1
        lines.append(f"  {n_label:2d} {shell.l:d} {shell.l:d} {nprim:3d} 1")
        for exp, coeff, coeff_env in zip(shell.exps, shell.coeffs, shell.coeffs_env):
            lines.append(f"    {exp:24.14f} {coeff:24.14f} {coeff_env:24.14f}")
    lines.append("#")


def append_tb_family(
    lines: list[str],
    title: str,
    names: list[str],
    basis: dict[int, tuple[tuple[float, float, float, float], list[Shell]]],
) -> None:
    lines.append("")
    lines.append(f"# ---- {title} ----")
    for z in sorted(basis):
        if z <= 0 or z >= len(SYMBOLS):
            raise ValueError(f"Unsupported atomic number {z}")
        params, shells = basis[z]
        write_tb_basis_block(lines, z, names, params, shells)


def validate_cp2k_basis_file(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if not stripped or stripped.startswith("#"):
            i += 1
            continue
        header = stripped.split()
        if len(header) < 2:
            raise ValueError(f"Bad basis header at line {i + 1}: {lines[i]!r}")
        zsym = header[0]
        names = header[1:]
        i += 1
        if i >= len(lines):
            raise ValueError(f"Missing nset after {zsym} {' '.join(names)}")
        nset = int(lines[i].strip())
        i += 1
        for _ in range(nset):
            fields = lines[i].split()
            if len(fields) < 5:
                raise ValueError(f"Bad set line at line {i + 1}: {lines[i]!r}")
            nexp = int(fields[3])
            nshell = sum(int(v) for v in fields[4:])
            if nshell <= 0:
                raise ValueError(f"No shell declared at line {i + 1}")
            i += 1 + nexp
        counts[names[0]] = counts.get(names[0], 0) + 1
    return counts


def validate_cp2k_tb_basis_file(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if not stripped or stripped.startswith("#"):
            i += 1
            continue
        header = stripped.split()
        if len(header) < 2:
            raise ValueError(f"Bad TB basis header at line {i + 1}: {lines[i]!r}")
        names = header[1:]
        i += 1
        params = lines[i].split()
        if len(params) != 6:
            raise ValueError(f"Bad TB parameter line at line {i + 1}: {lines[i]!r}")
        [float(value) for value in params]
        i += 1
        nset = int(lines[i].strip())
        i += 1
        for _ in range(nset):
            fields = lines[i].split()
            if len(fields) != 5 or fields[1] != fields[2] or fields[4] != "1":
                raise ValueError(
                    f"Unsupported TB set line at line {i + 1}: {lines[i]!r}"
                )
            nexp = int(fields[3])
            i += 1
            for _ in range(nexp):
                values = lines[i].split()
                if len(values) != 3:
                    raise ValueError(
                        f"Bad TB primitive line at line {i + 1}: {lines[i]!r}"
                    )
                if float(values[0]) <= 0.0:
                    raise ValueError(
                        f"Nonpositive exponent at line {i + 1}: {lines[i]!r}"
                    )
                i += 1
        counts[names[0]] = counts.get(names[0], 0) + 1
    return counts


def build() -> tuple[dict[str, int], dict[str, int]]:
    sto_table = parse_cp2k_sto_ng(CP2K_STO_NG)
    if not sto_table:
        raise ValueError(f"No STO-NG table parsed from {CP2K_STO_NG}")

    with tempfile.TemporaryDirectory(prefix="cp2k_tb_basis_") as tmp_dir:
        sources = clone_sources(Path(tmp_dir))
        refs = {name: repo_ref(path) for name, path in sources.items()}

        qvszp = parse_star_basis(sources["qvSZP"] / "q-vSZP_basis" / "basisq")
        qvszp_gxtb = parse_star_basis(sources["qvSZP"] / "q-vSZP_basis" / "basisq_gxtb")
        qvszp_gxtb_tb = parse_qvszp_tb_basis(
            sources["qvSZP"] / "q-vSZP_basis" / "basisq_gxtb"
        )
        ptb_vdzp = parse_star_basis(sources["ptb"] / ".basis_vDZP")
        ptb_aux_vdzp = parse_star_basis(sources["ptb"] / ".auxbasis_vDZP")
        gfn1 = parse_tblite_gfn(
            sources["save_tblite"] / "src" / "tblite" / "xtb" / "gfn1.f90", sto_table
        )
        gfn2 = parse_tblite_gfn(
            sources["save_tblite"] / "src" / "tblite" / "xtb" / "gfn2.f90", sto_table
        )

        date = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        lines = [
            "# CP2K BASIS_TB",
            f"# Generated by tools/generate_tb_basis.py on {date}",
            "#",
            "# This is a CP2K BASIS_SET-format file, intended for use via",
            "#   BASIS_SET_FILE_NAME BASIS_TB",
            "#",
            "# Sources:",
            f"# - qvSZP:       {REPOS['qvSZP']} @ {refs['qvSZP']}",
            f"# - PTB:         {REPOS['ptb']} @ {refs['ptb']}",
            f"# - save_tblite: {REPOS['save_tblite']} @ {refs['save_tblite']}",
            f"# - STO-NG:      CP2K {CP2K_STO_NG.relative_to(ROOT)}",
            "#",
            "# Limitations:",
            "# - The q-vSZP source basis carries charge/environment coefficients.",
            "#   CP2K BASIS_SET files are static, so these coefficients are not",
            "#   represented here. The first contraction coefficient column is used.",
            "# - qvSZP-gxTB is the g-xTB/q-xTB-adapted basis from basisq_gxtb,",
            "#   included here as the static Gaussian contractions provided upstream.",
            "# - GFN1/GFN2 entries are static STO-6G expansions of tblite Slater exponents.",
            "#   They are useful as Gaussian basis definitions, not as a replacement for",
            "#   CP2K's native xTB/tblite Hamiltonian and parameter handling.",
            "# - PTB ECP data are not included; this file contains orbital/auxiliary bases.",
        ]

        append_family(lines, "q-vSZP static orbital basis", ["qvSZP", "q-vSZP"], qvszp)
        append_family(
            lines,
            "q-vSZP g-xTB/q-xTB static orbital basis",
            ["qvSZP-gxTB", "q-vSZP-gxTB"],
            qvszp_gxtb,
        )
        append_family(lines, "PTB vDZP orbital basis", ["vDZP", "vDZP-PTB"], ptb_vdzp)
        append_family(
            lines,
            "PTB vDZP auxiliary basis",
            ["vDZP-AUX", "vDZP-AUX-PTB"],
            ptb_aux_vdzp,
        )
        append_family(
            lines,
            "tblite GFN1-xTB static STO-6G Gaussian basis",
            ["GFN1-xTB-STO6G", "GFN1-xTB", "GFN1"],
            gfn1,
        )
        append_family(
            lines,
            "tblite GFN2-xTB static STO-6G Gaussian basis",
            ["GFN2-xTB-STO6G", "GFN2-xTB", "GFN2"],
            gfn2,
        )

        OUT_FILE.write_text("\n".join(lines) + "\n")
        gxtb_lines = [
            "# CP2K gxTB_BASIS",
            f"# Generated by tools/generate_tb_basis.py on {date}",
            "#",
            "# This file is read by BASIS_SET ORB TB <name>.",
            "# Format:",
            "#   Element symbol  Name of the basis set  Alias names",
            "#   k0 k1 k2 k3 fixed_q fixed_cn",
            "#   nset",
            "#   n lmin lmax nexp 1",
            "#   exponent coefficient coefficient_env",
            "#",
            "# CP2K initializes the q/CN-dependent contractions at fixed_q and",
            "# fixed_cn. With TB_BASIS_UPDATE ITERATIVE the ORB TB reader updates",
            "# them from SCF Mulliken charges and geometry coordination numbers.",
            "#",
            "# Sources:",
            f"# - qvSZP: {REPOS['qvSZP']} @ {refs['qvSZP']}",
        ]
        append_tb_family(
            gxtb_lines,
            "q-vSZP g-xTB/q-xTB TB orbital basis",
            ["qvSZP-gxTB", "q-vSZP-gxTB"],
            qvszp_gxtb_tb,
        )
        GXTB_OUT_FILE.write_text("\n".join(gxtb_lines) + "\n")

    return validate_cp2k_basis_file(OUT_FILE), validate_cp2k_tb_basis_file(
        GXTB_OUT_FILE
    )


def main() -> int:
    if not shutil.which("git"):
        print("git is required to fetch source basis repositories", file=sys.stderr)
        return 2
    try:
        counts, tb_counts = build()
    except Exception as exc:  # noqa: BLE001 - command-line tool should report context
        print(f"Failed to generate {OUT_FILE}: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {OUT_FILE}")
    for name in sorted(counts):
        print(f"{name:16s} {counts[name]:4d} element blocks")
    print(f"Wrote {GXTB_OUT_FILE}")
    for name in sorted(tb_counts):
        print(f"{name:16s} {tb_counts[name]:4d} TB element blocks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
