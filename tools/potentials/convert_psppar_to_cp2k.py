#!/usr/bin/env python3
"""Convert HGH/GTH psppar files with optional NLCC data to CP2K format.

Example:
    python3 convert_psppar_to_cp2k.py psppar_dir MPC_POTENTIALS
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import math
import re


NUM_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?")
PROJECTOR_RE = re.compile(r"\b([spdfg])-projector\b")


ATOMIC_NUMBERS = {
    "H": 1,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Y": 39,
    "Au": 79,
}


ORBITALS = {
    # Only elements present in Yike.tar are needed here. Occupancies are neutral
    # atomic ground states; core electrons are stripped from the innermost
    # orbitals to recover CP2K's n_elec(s), n_elec(p), ... line.
    "H": [(1, "s", 1)],
    "Li": [(1, "s", 2), (2, "s", 1)],
    "Be": [(1, "s", 2), (2, "s", 2)],
    "B": [(1, "s", 2), (2, "s", 2), (2, "p", 1)],
    "C": [(1, "s", 2), (2, "s", 2), (2, "p", 2)],
    "N": [(1, "s", 2), (2, "s", 2), (2, "p", 3)],
    "O": [(1, "s", 2), (2, "s", 2), (2, "p", 4)],
    "F": [(1, "s", 2), (2, "s", 2), (2, "p", 5)],
    "Na": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 1)],
    "Mg": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2)],
    "Al": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 1)],
    "Si": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 2)],
    "P": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 3)],
    "S": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 4)],
    "Cl": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 5)],
    "K": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (4, "s", 1)],
    "Ca": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (4, "s", 2)],
    "Sc": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 1), (4, "s", 2)],
    "Ti": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 2), (4, "s", 2)],
    "V": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 3), (4, "s", 2)],
    "Cr": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 5), (4, "s", 1)],
    "Mn": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 5), (4, "s", 2)],
    "Fe": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 6), (4, "s", 2)],
    "Co": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 7), (4, "s", 2)],
    "Ni": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 8), (4, "s", 2)],
    "Cu": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 10), (4, "s", 1)],
    "Zn": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 10), (4, "s", 2)],
    "Ga": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 10), (4, "s", 2), (4, "p", 1)],
    "Ge": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 10), (4, "s", 2), (4, "p", 2)],
    "As": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 10), (4, "s", 2), (4, "p", 3)],
    "Se": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 10), (4, "s", 2), (4, "p", 4)],
    "Br": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 10), (4, "s", 2), (4, "p", 5)],
    "Y": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6), (3, "d", 10), (4, "s", 2), (4, "p", 6), (4, "d", 1), (5, "s", 2)],
    "Au": [
        (1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6),
        (3, "d", 10), (4, "s", 2), (4, "p", 6), (4, "d", 10), (4, "f", 14),
        (5, "s", 2), (5, "p", 6), (5, "d", 10), (6, "s", 1),
    ],
}

L_ORDER = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}


@dataclass
class Projector:
    l: str
    radius: float
    nproj: int
    h: list[float]
    soc: list[float]


@dataclass
class PspparPotential:
    source_name: str
    element: str
    zatom: int
    zion: int
    date: str
    pspcod: int
    ixc: int
    rloc: float
    nloc: int
    local_coeffs: list[float]
    projectors: list[Projector]
    rcore: float | None
    qcore: float | None

    @property
    def has_nlcc(self) -> bool:
        return self.rcore is not None and self.qcore is not None and self.zatom > self.zion

    @property
    def c_core(self) -> float | None:
        if not self.has_nlcc:
            return None
        return (
            self.qcore
            * 4.0
            * math.pi
            * (self.zatom - self.zion)
            / ((math.sqrt(2.0 * math.pi) * self.rcore) ** 3)
        )


def numbers(line: str) -> list[float]:
    return [float(match.group(0)) for match in NUM_RE.finditer(line)]


def element_from_path(path: Path) -> str:
    element = path.name.rsplit(".", 1)[-1]
    if element not in ATOMIC_NUMBERS:
        raise ValueError(f"Cannot determine supported element from {path}")
    return element


def valence_counts(element: str, zion: int) -> list[int]:
    zatom = ATOMIC_NUMBERS[element]
    core = zatom - zion
    if core < 0:
        raise ValueError(f"{element}: zion={zion} exceeds Z={zatom}")
    remaining_core = core
    counts = {l: 0 for l in L_ORDER}
    for _n, l, occupancy in sorted(ORBITALS[element], key=lambda item: (item[0], L_ORDER[item[1]])):
        removed = min(occupancy, remaining_core)
        remaining_core -= removed
        counts[l] += occupancy - removed
    if remaining_core:
        raise ValueError(f"{element}: could not strip {core} core electrons")
    highest_l = max((L_ORDER[l] for l, occ in counts.items() if occ), default=0)
    ordered = ["s", "p", "d", "f", "g"][: highest_l + 1]
    result = [counts[l] for l in ordered]
    if sum(result) != zion:
        raise ValueError(f"{element}: valence counts {result} do not sum to zion={zion}")
    return result


def parse_psppar(path: Path) -> PspparPotential:
    lines = path.read_text().splitlines()
    element = element_from_path(path)

    zline = numbers(lines[1])
    zatom, zion, date = int(round(zline[0])), int(round(zline[1])), str(int(round(zline[2])))
    psp_line = numbers(lines[2])
    pspcod, ixc = int(round(psp_line[0])), int(round(psp_line[1]))
    local = numbers(lines[3])
    rloc, nloc = local[0], int(round(local[1]))
    local_coeffs = local[2 : 2 + nloc]
    if len(local_coeffs) != nloc:
        raise ValueError(f"{path}: expected {nloc} local coefficients")
    nsep = int(round(numbers(lines[4])[0]))

    idx = 5
    projectors: list[Projector] = []
    for _ in range(nsep):
        while idx < len(lines) and PROJECTOR_RE.search(lines[idx]) is None:
            idx += 1
        if idx >= len(lines):
            raise ValueError(f"{path}: missing projector block")

        match = PROJECTOR_RE.search(lines[idx])
        assert match is not None
        label_l = match.group(1)
        vals = numbers(lines[idx])
        radius, nproj = vals[0], int(round(vals[1]))
        needed = nproj * (nproj + 1) // 2
        h_values = vals[2:]
        soc_values: list[float] = []
        idx += 1

        while idx < len(lines) and PROJECTOR_RE.search(lines[idx]) is None and "rcore" not in lines[idx]:
            vals = numbers(lines[idx])
            missing = needed - len(h_values)
            if missing > 0:
                h_values.extend(vals[:missing])
                soc_values.extend(vals[missing:])
            else:
                soc_values.extend(vals)
            idx += 1

        if len(h_values) != needed:
            raise ValueError(f"{path}: {label_l}-projector expected {needed} h values, found {len(h_values)}")
        projectors.append(Projector(label_l, radius, nproj, h_values, soc_values))

    rcore = qcore = None
    while idx < len(lines):
        if "rcore" in lines[idx]:
            vals = numbers(lines[idx])
            rcore, qcore = vals[0], vals[1]
            break
        idx += 1

    return PspparPotential(
        source_name=path.name,
        element=element,
        zatom=zatom,
        zion=zion,
        date=date,
        pspcod=pspcod,
        ixc=ixc,
        rloc=rloc,
        nloc=nloc,
        local_coeffs=local_coeffs,
        projectors=projectors,
        rcore=rcore,
        qcore=qcore,
    )


def format_float(value: float) -> str:
    return f"{value:20.14f}"


def format_projector(projector: Projector) -> list[str]:
    lines: list[str] = []
    cursor = 0
    for row in range(projector.nproj):
        row_values = projector.h[cursor : cursor + projector.nproj - row]
        cursor += projector.nproj - row
        values = "".join(format_float(value) for value in row_values)
        if row == 0:
            lines.append(f"    {projector.radius:20.14f} {projector.nproj:2d}{values}")
        else:
            lines.append(f"    {'':23s}{values}")
    return lines


def format_cp2k_potential(potential: PspparPotential) -> str:
    lines = [
        f"# source: {potential.source_name}; zatom={potential.zatom}; zion={potential.zion}; date={potential.date}; pspcod={potential.pspcod}; ixc={potential.ixc}",
        f"{potential.element} GTH-NLCC-PBE-q{potential.zion}",
        "    " + " ".join(f"{value:4d}" for value in valence_counts(potential.element, potential.zion)),
        f"    {potential.rloc:20.14f} {potential.nloc:2d}" + "".join(format_float(value) for value in potential.local_coeffs),
    ]
    if potential.has_nlcc:
        lines.extend(
            [
                "    NLCC  1",
                f"    {potential.rcore:20.14f}  1{format_float(potential.c_core or 0.0)}",
            ]
        )
    lines.append(f"    {len(potential.projectors):2d}")
    for projector in potential.projectors:
        lines.extend(format_projector(projector))
    return "\n".join(lines)


def parse_cp2k_entries(path: Path) -> dict[tuple[str, str], dict[str, object]]:
    raw_lines = path.read_text().splitlines()
    lines = [line for line in raw_lines if line.strip() and not line.lstrip().startswith("#")]
    entries: dict[tuple[str, str], dict[str, object]] = {}
    i = 0
    while i < len(lines):
        header = lines[i].split()
        if len(header) < 2 or not re.fullmatch(r"[A-Z][a-z]?", header[0]):
            i += 1
            continue
        element, name = header[0], header[1]
        i += 1
        config = [int(round(value)) for value in numbers(lines[i])]
        i += 1
        local_values = numbers(lines[i])
        rloc, nloc = local_values[0], int(round(local_values[1]))
        local_coeffs = local_values[2 : 2 + nloc]
        i += 1
        nlcc = None
        if lines[i].strip().startswith("NLCC"):
            i += 1
            vals = numbers(lines[i])
            nlcc = vals[:3]
            i += 1
        nprj = int(round(numbers(lines[i])[0]))
        i += 1
        projectors = []
        for _ in range(nprj):
            vals = numbers(lines[i])
            radius, nproj = vals[0], int(round(vals[1]))
            needed = nproj * (nproj + 1) // 2
            h_values = vals[2:]
            i += 1
            while len(h_values) < needed:
                h_values.extend(numbers(lines[i]))
                i += 1
            projectors.append([radius, nproj, h_values[:needed]])
        entries[(element, name)] = {
            "element": element,
            "name": name,
            "config": config,
            "local": [rloc, nloc, local_coeffs],
            "nlcc": nlcc,
            "projectors": projectors,
        }
    return entries


def numeric_vector(entry: dict[str, object]) -> list[float]:
    vector: list[float] = []
    vector.extend(float(value) for value in entry["config"])  # type: ignore[arg-type]
    rloc, nloc, coeffs = entry["local"]  # type: ignore[misc]
    vector.extend([float(rloc), float(nloc)])
    vector.extend(float(value) for value in coeffs)
    nlcc = entry["nlcc"]
    if nlcc is not None:
        vector.extend(float(value) for value in nlcc)  # type: ignore[arg-type]
    else:
        vector.append(-1.0)
    for radius, nproj, values in entry["projectors"]:  # type: ignore[union-attr]
        vector.extend([float(radius), float(nproj)])
        vector.extend(float(value) for value in values)
    return vector


def max_abs_delta(a: dict[str, object], b: dict[str, object]) -> float | None:
    va = numeric_vector(a)
    vb = numeric_vector(b)
    if len(va) != len(vb):
        return None
    return max(abs(x - y) for x, y in zip(va, vb))


def write_comparison(converted_path: Path, reference_path: Path, out_path: Path) -> None:
    converted = parse_cp2k_entries(converted_path)
    reference = parse_cp2k_entries(reference_path)
    ref_by_element: dict[str, list[str]] = {}
    for element, name in reference:
        ref_by_element.setdefault(element, []).append(name)

    rows = [
        "| Element | Converted potential | Valence | NLCC | CP2K reference | Comparison |",
        "|---|---:|---:|---:|---|---|",
    ]
    for element, name in sorted(converted, key=lambda item: (ATOMIC_NUMBERS[item[0]], item[1])):
        entry = converted[(element, name)]
        valence = sum(entry["config"])  # type: ignore[arg-type]
        has_nlcc = "yes" if entry["nlcc"] is not None else "no"
        ref_names = sorted(ref_by_element.get(element, []))
        same_name = (element, name) in reference
        if same_name:
            delta = max_abs_delta(entry, reference[(element, name)])
            if delta is not None and delta < 1e-10:
                comparison = "identical to same named CP2K entry"
            elif delta is None:
                comparison = "same name, different structure"
            else:
                comparison = f"same name, parameters differ (max |delta| {delta:.3g})"
        elif ref_names:
            comparison = "same element, different q/name"
        else:
            comparison = "new element for this reference file"
        rows.append(
            f"| {element} | `{name}` | {valence} | {has_nlcc} | "
            f"{', '.join(f'`{item}`' for item in ref_names) if ref_names else '-'} | {comparison} |"
        )
    out_path.write_text("\n".join(rows) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_file", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--comparison", type=Path)
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob("psppar*"), key=lambda path: ATOMIC_NUMBERS[element_from_path(path)])
    potentials = [parse_psppar(path) for path in paths]

    nonzero_soc = [
        (potential.element, projector.l, value)
        for potential in potentials
        for projector in potential.projectors
        for value in projector.soc
        if abs(value) > 1e-12
    ]
    if nonzero_soc:
        detail = ", ".join(f"{el}/{l}={value:g}" for el, l, value in nonzero_soc[:10])
        raise ValueError(f"Non-zero spin-orbit projector data found; scalar conversion would lose data: {detail}")

    header = """################################################################################
#
# MPC GTH pseudopotentials converted from HGH/GTH psppar files.
# Source format: HGH/GTH psppar, PBE (ixc=-101130).
# CP2K c_core conversion:
#   c_core = qcore * 4*pi*(Z-Zion)/(sqrt(2*pi)*r_core)^3
#
# Entries intentionally omit the generic GTH-NLCC-PBE alias to avoid ambiguity
# with CP2K's existing NLCC_POTENTIALS entries for the same elements.
#
################################################################################
"""
    args.output_file.write_text(header + "\n# PBE functional\n#\n" + "\n#\n".join(format_cp2k_potential(p) for p in potentials) + "\n")

    if args.reference and args.comparison:
        write_comparison(args.output_file, args.reference, args.comparison)


if __name__ == "__main__":
    main()
