#!/usr/bin/env python3
"""Compare TDDFPT k-point calculations against Gamma supercells.

This is a research/debugging helper for the TDDFPT k-point kernel branch.  It
generates small molecular and crystalline inputs, runs KERNEL NONE/FULL with a
k-point mesh and with the corresponding Gamma supercell, then prints the first
TDDFPT states.
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


@dataclass(frozen=True)
class Atom:
    element: str
    fractional: tuple[float, float, float]


@dataclass(frozen=True)
class SystemDefinition:
    name: str
    description: str
    cell: tuple[
        tuple[float, float, float],
        tuple[float, float, float],
        tuple[float, float, float],
    ]
    atoms: tuple[Atom, ...]
    basis: dict[str, str]
    potential: dict[str, str]
    cutoff: float = 100.0
    rel_cutoff: float = 30.0


@dataclass(frozen=True)
class EnergyGroup:
    first_state: int
    last_state: int
    count: int
    energy_min: float
    energy_max: float
    energy_avg: float
    oscillator_sum: float


@dataclass(frozen=True)
class SpectrumBlock:
    kpoint_groups: tuple[EnergyGroup, ...]
    supercell_groups: tuple[EnergyGroup, ...]


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


def fcc_primitive_cell(
    a: float,
) -> tuple[
    tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]
]:
    return ((0.0, 0.5 * a, 0.5 * a), (0.5 * a, 0.0, 0.5 * a), (0.5 * a, 0.5 * a, 0.0))


SYSTEMS: dict[str, SystemDefinition] = {
    "h2": SystemDefinition(
        name="h2",
        description="H2 molecule in a periodically repeated large box",
        cell=((4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)),
        atoms=(
            Atom("H", (1.65 / 4.0, 0.5, 0.5)),
            Atom("H", (2.35 / 4.0, 0.5, 0.5)),
        ),
        basis={"H": "DZVP-MOLOPT-GGA-GTH-q1"},
        potential={"H": "GTH-GGA-q1"},
    ),
    "si": SystemDefinition(
        name="si",
        description="primitive diamond silicon cell",
        cell=fcc_primitive_cell(5.43),
        atoms=(Atom("Si", (0.0, 0.0, 0.0)), Atom("Si", (0.25, 0.25, 0.25))),
        basis={"Si": "DZVP-MOLOPT-GGA-GTH-q4"},
        potential={"Si": "GTH-GGA-q4"},
    ),
    "c": SystemDefinition(
        name="c",
        description="primitive diamond carbon cell",
        cell=fcc_primitive_cell(3.57),
        atoms=(Atom("C", (0.0, 0.0, 0.0)), Atom("C", (0.25, 0.25, 0.25))),
        basis={"C": "DZVP-MOLOPT-GGA-GTH-q4"},
        potential={"C": "GTH-GGA-q4"},
    ),
    "bn": SystemDefinition(
        name="bn",
        description="primitive zincblende boron nitride cell",
        cell=fcc_primitive_cell(3.62),
        atoms=(Atom("B", (0.0, 0.0, 0.0)), Atom("N", (0.25, 0.25, 0.25))),
        basis={"B": "DZVP-MOLOPT-GGA-GTH-q3", "N": "DZVP-MOLOPT-GGA-GTH-q5"},
        potential={"B": "GTH-GGA-q3", "N": "GTH-GGA-q5"},
    ),
    "lif": SystemDefinition(
        name="lif",
        description="primitive rocksalt lithium fluoride cell",
        cell=fcc_primitive_cell(4.03),
        atoms=(Atom("Li", (0.0, 0.0, 0.0)), Atom("F", (0.5, 0.5, 0.5))),
        basis={"Li": "DZVP-MOLOPT-GGA-GTH-q1", "F": "DZVP-MOLOPT-GGA-GTH-q7"},
        potential={"Li": "GTH-GGA-q1", "F": "GTH-GGA-q7"},
    ),
    "nacl": SystemDefinition(
        name="nacl",
        description="primitive rocksalt sodium chloride cell",
        cell=fcc_primitive_cell(5.64),
        atoms=(Atom("Na", (0.0, 0.0, 0.0)), Atom("Cl", (0.5, 0.5, 0.5))),
        basis={"Na": "DZVP-MOLOPT-GGA-GTH-q1", "Cl": "DZVP-MOLOPT-GGA-GTH-q7"},
        potential={"Na": "GTH-GGA-q1", "Cl": "GTH-GGA-q7"},
    ),
}


def scaled_vector(
    scale: int, vector: tuple[float, float, float]
) -> tuple[float, float, float]:
    return (scale * vector[0], scale * vector[1], scale * vector[2])


def frac_to_cart(
    cell: tuple[
        tuple[float, float, float],
        tuple[float, float, float],
        tuple[float, float, float],
    ],
    fractional: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        sum(fractional[idir] * cell[idir][icomponent] for idir in range(3))
        for icomponent in range(3)
    )


def format_cell(
    cell: tuple[
        tuple[float, float, float],
        tuple[float, float, float],
        tuple[float, float, float],
    ],
) -> str:
    labels = ("A", "B", "C")
    return "\n".join(
        f"      {label} {vector[0]:12.7f} {vector[1]:12.7f} {vector[2]:12.7f}"
        for label, vector in zip(labels, cell)
    )


def system_cell(
    system: SystemDefinition, grid: tuple[int, int, int], supercell: bool
) -> tuple[tuple[float, float, float], ...]:
    if not supercell:
        return system.cell
    return tuple(scaled_vector(grid[idir], system.cell[idir]) for idir in range(3))


def system_coords(
    system: SystemDefinition, grid: tuple[int, int, int], supercell: bool
) -> str:
    nx, ny, nz = grid if supercell else (1, 1, 1)
    lines: list[str] = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for atom in system.atoms:
                    fractional = (
                        atom.fractional[0] + ix,
                        atom.fractional[1] + iy,
                        atom.fractional[2] + iz,
                    )
                    x, y, z = frac_to_cart(system.cell, fractional)
                    lines.append(
                        f"      {atom.element:2s} {x:12.7f} {y:12.7f} {z:12.7f}"
                    )
    return "\n".join(lines)


def kind_sections(system: SystemDefinition) -> str:
    lines: list[str] = []
    for element in sorted(system.basis):
        lines.extend(
            [
                f"    &KIND {element}",
                f"      BASIS_SET ORB {system.basis[element]}",
                f"      POTENTIAL {system.potential[element]}",
                "    &END KIND",
            ]
        )
    return "\n".join(lines)


def make_input(
    *,
    system: SystemDefinition,
    project: str,
    kernel: str,
    functional: str,
    nstates: int,
    added_mos: int,
    grid: tuple[int, int, int],
    gamma_centered: bool,
    supercell: bool,
    triplet: bool,
    dipole_operator: str,
) -> str:
    supercell_scale = grid[0] * grid[1] * grid[2] if supercell else 1
    run_added_mos = added_mos * supercell_scale
    run_nstates = nstates * supercell_scale
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
    dipole_section = ""
    if not supercell and dipole_operator != "VELOCITY":
        dipole_section = f"""      &DIPOLE_MOMENTS
        DIPOLE_FORM {dipole_operator}
      &END DIPOLE_MOMENTS
"""
    cell = system_cell(system, grid, supercell)
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
      CUTOFF {system.cutoff:.0f}
      REL_CUTOFF {system.rel_cutoff:.0f}
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-10
      METHOD GPW
    &END QS
    &SCF
      ADDED_MOS {run_added_mos}
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
      NSTATES {run_nstates}
{triplet_line}{dipole_section}    &END TDDFPT
  &END PROPERTIES
  &SUBSYS
    &CELL
{format_cell(cell)}
    &END CELL
    &COORD
{system_coords(system, grid, supercell)}
    &END COORD
{kind_sections(system)}
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
    env_overrides: dict[str, str] | None,
    keep_going: bool,
) -> Result | None:
    inp = work_dir / f"{label}.inp"
    out = work_dir / f"{label}.out"
    inp.write_text(input_text)
    env = os.environ.copy()
    env["CP2K_DATA_DIR"] = str(data_dir)
    env.setdefault("OMP_NUM_THREADS", "1")
    if env_overrides:
        env.update(env_overrides)
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


def group_states(
    states: tuple[State, ...], tolerance_ev: float
) -> tuple[EnergyGroup, ...]:
    if not states:
        return ()

    groups: list[EnergyGroup] = []
    current: list[State] = [states[0]]
    group_start = states[0].energy_ev
    for state in states[1:]:
        if abs(state.energy_ev - group_start) <= tolerance_ev:
            current.append(state)
            continue
        groups.append(make_group(current))
        current = [state]
        group_start = state.energy_ev
    groups.append(make_group(current))
    return tuple(groups)


def make_group(states: list[State]) -> EnergyGroup:
    energies = [state.energy_ev for state in states]
    return EnergyGroup(
        first_state=states[0].index,
        last_state=states[-1].index,
        count=len(states),
        energy_min=min(energies),
        energy_max=max(energies),
        energy_avg=sum(energies) / len(energies),
        oscillator_sum=sum(state.oscillator for state in states),
    )


def group_range(group: EnergyGroup) -> str:
    if group.first_state == group.last_state:
        return f"{group.first_state}"
    return f"{group.first_state}-{group.last_state}"


def groups_range(groups: tuple[EnergyGroup, ...]) -> str:
    if not groups:
        return "-"
    return ",".join(group_range(group) for group in groups)


def groups_count(groups: tuple[EnergyGroup, ...]) -> int:
    return sum(group.count for group in groups)


def groups_energy_avg(groups: tuple[EnergyGroup, ...]) -> float | None:
    count = groups_count(groups)
    if count == 0:
        return None
    return sum(group.energy_avg * group.count for group in groups) / count


def groups_oscillator_sum(groups: tuple[EnergyGroup, ...]) -> float:
    return sum(group.oscillator_sum for group in groups)


def print_group_spectrum(result: Result, tolerance_ev: float) -> None:
    groups = group_states(result.states, tolerance_ev)
    print(f"\nGrouped spectrum: {result.label}  (tol={tolerance_ev:.1e} eV)")
    print("  states         n       e_min/eV       e_max/eV       osc_sum")
    for group in groups:
        print(
            f"  {group_range(group):>8s}  {group.count:4d}"
            f"  {group.energy_min:13.6f}  {group.energy_max:13.6f}"
            f"  {group.oscillator_sum:12.5e}"
        )
    print(f"  total oscillator over printed states: {sum_group_osc(groups):.8e}")


def sum_group_osc(groups: tuple[EnergyGroup, ...]) -> float:
    return sum(group.oscillator_sum for group in groups)


def folded_spectrum_blocks(
    kpoint: Result,
    supercell: Result,
    group_tolerance_ev: float,
    match_tolerance_ev: float,
) -> tuple[SpectrumBlock, ...]:
    items: list[tuple[float, str, EnergyGroup]] = []
    for group in group_states(kpoint.states, group_tolerance_ev):
        items.append((group.energy_avg, "kpoint", group))
    for group in group_states(supercell.states, group_tolerance_ev):
        items.append((group.energy_avg, "supercell", group))
    if not items:
        return ()

    blocks: list[SpectrumBlock] = []
    current: list[tuple[float, str, EnergyGroup]] = []
    current_max = -float("inf")

    def flush(block_items: list[tuple[float, str, EnergyGroup]]) -> None:
        kpoint_groups = tuple(
            group for _energy, side, group in block_items if side == "kpoint"
        )
        supercell_groups = tuple(
            group for _energy, side, group in block_items if side == "supercell"
        )
        blocks.append(SpectrumBlock(kpoint_groups, supercell_groups))

    for item in sorted(items, key=lambda entry: entry[0]):
        energy = item[0]
        if not current or energy - current_max <= match_tolerance_ev:
            current.append(item)
            current_max = max(current_max, energy)
            continue
        flush(current)
        current = [item]
        current_max = energy
    flush(current)
    return tuple(blocks)


def format_optional(value: float | None, width: int = 8, precision: int = 4) -> str:
    if value is None:
        return f"{'-':>{width}s}"
    return f"{value:{width}.{precision}f}"


def compare_groups(
    kpoint: Result,
    supercell: Result,
    group_tolerance_ev: float,
    match_tolerance_ev: float,
    supercell_scale: int,
) -> None:
    blocks = folded_spectrum_blocks(
        kpoint, supercell, group_tolerance_ev, match_tolerance_ev
    )
    print(f"\nDegeneracy/folding-aware spectra: {kpoint.label} -> {supercell.label}")
    if supercell_scale > 1:
        print(
            "  kp states    nk     k/eV   sc states    nsc    sc/eV"
            "      dk/eV    nsc/nk     osc(k)      osc(sc)    osc(sc/N)"
            "   k/(sc/N)"
        )
    else:
        print(
            "  kp states    nk     k/eV   sc states    nsc    sc/eV"
            "      dk/eV    nsc/nk     osc(k)      osc(sc)"
        )
    for block in blocks:
        left_energy = groups_energy_avg(block.kpoint_groups)
        right_energy = groups_energy_avg(block.supercell_groups)
        dk = None
        if left_energy is not None and right_energy is not None:
            dk = left_energy - right_energy
        left_count = groups_count(block.kpoint_groups)
        right_count = groups_count(block.supercell_groups)
        fold = None
        if left_count > 0 and right_count > 0:
            fold = right_count / left_count
        left_osc = groups_oscillator_sum(block.kpoint_groups)
        right_osc = groups_oscillator_sum(block.supercell_groups)
        line = (
            f"  {groups_range(block.kpoint_groups):>9s}  {left_count:4d}"
            f"  {format_optional(left_energy)}"
            f"  {groups_range(block.supercell_groups):>9s}  {right_count:5d}"
            f"  {format_optional(right_energy)}  {format_optional(dk, 9, 5)}"
            f"  {format_optional(fold, 8, 3)}"
            f"  {left_osc:10.3e}  {right_osc:10.3e}"
        )
        if supercell_scale > 1:
            right_scaled = right_osc / supercell_scale
            ratio = None
            if abs(right_scaled) > 1.0e-12:
                ratio = left_osc / right_scaled
            line += f"  {right_scaled:10.3e}  {format_optional(ratio, 9, 3)}"
        print(line)
    left_total = sum(groups_oscillator_sum(block.kpoint_groups) for block in blocks)
    right_total = sum(groups_oscillator_sum(block.supercell_groups) for block in blocks)
    matched_left_total = sum(
        groups_oscillator_sum(block.kpoint_groups)
        for block in blocks
        if block.kpoint_groups and block.supercell_groups
    )
    matched_right_total = sum(
        groups_oscillator_sum(block.supercell_groups)
        for block in blocks
        if block.kpoint_groups and block.supercell_groups
    )
    if supercell_scale > 1:
        right_total_scaled = right_total / supercell_scale
        matched_right_scaled = matched_right_total / supercell_scale
        print(
            f"  matched oscillator difference vs sc/N:"
            f" {matched_left_total - matched_right_scaled:.8e}"
        )
        print(
            f"  total oscillator difference vs sc/N:"
            f" {left_total - right_total_scaled:.8e}"
        )
    print(
        "  total oscillator difference over printed states:"
        f" {left_total - right_total:.8e}"
    )


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

    print(f"\nNearest-energy matches: {supercell.label} -> {kpoint.label}")
    print("  sc     sc/eV   kp       k/eV       dk/eV      osc(sc)       osc(k)")
    for right in supercell.states[:limit]:
        left = min(
            kpoint.states, key=lambda state: abs(right.energy_ev - state.energy_ev)
        )
        print(
            f"  {right.index:2d}  {right.energy_ev:10.6f}"
            f"  {left.index:3d}  {left.energy_ev:10.6f}"
            f"  {left.energy_ev - right.energy_ev:10.6f}"
            f"  {right.oscillator:12.5e}  {left.oscillator:12.5e}"
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
        "--system",
        default="h2",
        choices=tuple(SYSTEMS),
        help="Small molecular or crystalline system to generate.",
    )
    parser.add_argument(
        "--grid", nargs=3, default=("2", "1", "1"), metavar=("NX", "NY", "NZ")
    )
    parser.add_argument("--functional", default="PBE", choices=("PBE", "PADE", "LDA"))
    parser.add_argument("--nstates", type=int, default=4)
    parser.add_argument("--added-mos", type=int, default=4)
    parser.add_argument(
        "--kernels", nargs="+", default=("NONE", "FULL"), choices=("NONE", "FULL")
    )
    parser.add_argument(
        "--kernel-parts",
        nargs="+",
        default=(),
        choices=("FULL", "GAP", "XC", "HARTREE"),
        help="Run k-point KERNEL FULL with internal kernel-part diagnostics.",
    )
    parser.add_argument("--density-weight-power", type=float)
    parser.add_argument("--sigma-weight-power", type=float)
    parser.add_argument("--density-scale", type=float)
    parser.add_argument("--sigma-scale", type=float)
    parser.add_argument("--no-gamma-centered", action="store_true")
    parser.add_argument("--triplet", action="store_true")
    parser.add_argument("--bright-threshold", type=float, default=1.0e-8)
    parser.add_argument(
        "--group-tol-ev",
        type=float,
        default=1.0e-4,
        help="Energy tolerance for grouping near-degenerate TDDFPT states.",
    )
    parser.add_argument(
        "--match-tol-ev",
        type=float,
        default=5.0e-2,
        help="Energy tolerance for matching k-point blocks to folded Gamma blocks.",
    )
    parser.add_argument(
        "--dipole-operator",
        default="VELOCITY",
        choices=("VELOCITY", "SCF_MOMENT"),
        help="Diagnostic k-point transition-dipole operator for the k-point run.",
    )
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
    system = SYSTEMS[args.system]
    work_dir = args.work_dir or Path(tempfile.mkdtemp(prefix="cp2k_tddfpt_kp_compare_"))
    work_dir = work_dir.resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    if not cp2k.exists():
        raise SystemExit(f"CP2K binary not found: {cp2k}")
    if not data_dir.exists():
        raise SystemExit(f"CP2K data directory not found: {data_dir}")

    print(f"Work dir: {work_dir}")
    print(f"System:   {system.name} ({system.description})")
    print(f"Grid:     {grid[0]} {grid[1]} {grid[2]}")
    print(f"Fold N:   {grid[0] * grid[1] * grid[2]}")
    print(f"XC:       {args.functional}")
    print(f"Triplet:  {args.triplet}")
    print(f"Dipoles:  {args.dipole_operator}")
    print(f"Match:   {args.match_tol_ev:.1e} eV")
    if args.kernel_parts:
        print(f"Parts:    {' '.join(args.kernel_parts)}")

    debug_env: dict[str, str] = {}
    if args.density_weight_power is not None:
        debug_env["CP2K_TDDFPT_KP_DENSITY_WEIGHT_POWER"] = str(
            args.density_weight_power
        )
    if args.sigma_weight_power is not None:
        debug_env["CP2K_TDDFPT_KP_SIGMA_WEIGHT_POWER"] = str(args.sigma_weight_power)
    if args.density_scale is not None:
        debug_env["CP2K_TDDFPT_KP_DENSITY_SCALE"] = str(args.density_scale)
    if args.sigma_scale is not None:
        debug_env["CP2K_TDDFPT_KP_SIGMA_SCALE"] = str(args.sigma_scale)

    results: dict[str, Result] = {}
    for kernel in args.kernels:
        for supercell in (False, True):
            if kernel == "FULL" and not supercell:
                parts = args.kernel_parts or ("",)
            else:
                parts = ("",)
            for part in parts:
                suffix = f"_{part.lower()}" if part else ""
                label = (
                    f"{system.name}_{'gamma_supercell' if supercell else 'kpoint'}_"
                    f"{kernel.lower()}{suffix}"
                )
                input_text = make_input(
                    system=system,
                    project=label,
                    kernel=kernel,
                    functional=args.functional,
                    nstates=args.nstates,
                    added_mos=args.added_mos,
                    grid=grid,
                    gamma_centered=not args.no_gamma_centered,
                    supercell=supercell,
                    triplet=args.triplet,
                    dipole_operator=args.dipole_operator,
                )
                env_overrides = dict(debug_env)
                if part:
                    env_overrides["CP2K_TDDFPT_KP_KERNEL_PART"] = part
                result = run_case(
                    cp2k=cp2k,
                    data_dir=data_dir,
                    work_dir=work_dir,
                    label=label,
                    input_text=input_text,
                    env_overrides=env_overrides,
                    keep_going=args.keep_going,
                )
                if result is not None:
                    results[label] = result

    for result in results.values():
        print_result(result, len(result.states), args.bright_threshold)
        print_group_spectrum(result, args.group_tol_ev)

    for kernel in args.kernels:
        left = results.get(f"{system.name}_kpoint_{kernel.lower()}")
        right = results.get(f"{system.name}_gamma_supercell_{kernel.lower()}")
        if left is not None and right is not None:
            compare_pair(left, right, len(left.states))
            compare_groups(
                left,
                right,
                args.group_tol_ev,
                args.match_tol_ev,
                grid[0] * grid[1] * grid[2],
            )
        if kernel == "FULL" and args.kernel_parts and right is not None:
            for part in args.kernel_parts:
                part_result = results.get(f"{system.name}_kpoint_full_{part.lower()}")
                if part_result is not None:
                    compare_pair(part_result, right, len(part_result.states))
                    compare_groups(
                        part_result,
                        right,
                        args.group_tol_ev,
                        args.match_tol_ev,
                        grid[0] * grid[1] * grid[2],
                    )

    print(f"\nKept inputs/outputs in {work_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
