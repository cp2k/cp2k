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
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

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
SPRD_RE = re.compile(
    r"O_D=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s+"
    r"O_OD=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s+"
    r"O_TOT=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
)
ERROR_RE = re.compile(r"(^|\s)(ERROR:|Error in|Exiting\.\.\.\.\.\.\.)")
ITERATION_LABEL_RE = re.compile(r"^(?P<variant>.+)-iter-(?P<iteration>\d+)$")
CONTROLLED_WIN_KEYS = {
    "auto_projections",
    "guiding_centres",
    "guiding_centers",
    "num_iter",
    "use_bloch_phases",
}
VARIANT_CHOICES = (
    "identity",
    "zero-iter",
    "identity-amn",
    "identity-amn-zero-iter",
    "phase-amn",
    "phase-amn-zero-iter",
    "smooth-amn",
    "smooth-amn-zero-iter",
    "neighbor-smooth-amn",
    "neighbor-smooth-amn-zero-iter",
    "global-smooth-amn",
    "global-smooth-amn-zero-iter",
    "unitary-smooth-amn",
    "unitary-smooth-amn-zero-iter",
    "transport-unitary-amn",
    "transport-unitary-amn-zero-iter",
    "transport-polar-amn",
    "transport-polar-amn-zero-iter",
    "transport-unitary-guiding-centres",
    "transport-polar-guiding-centres",
    "scdm-polar-amn",
    "scdm-polar-amn-zero-iter",
    "qr-pivot-amn",
    "qr-pivot-amn-zero-iter",
    "native-auto-projections",
    "native-auto-projections-zero-iter",
    "external-bloch",
    "external-bloch-zero-iter",
    "guiding-centres",
    "guiding-centres-zero-iter",
    "identity-zero-iter",
)
DEFAULT_SCAN_VARIANTS = (
    "zero-iter",
    "identity-amn",
    "identity-amn-zero-iter",
    "scdm-polar-amn",
    "scdm-polar-amn-zero-iter",
    "transport-unitary-amn",
    "transport-unitary-amn-zero-iter",
    "transport-polar-amn",
    "transport-polar-amn-zero-iter",
)
DEFAULT_ITERATION_SCAN_VARIANTS = (
    "scdm-polar-amn",
    "transport-unitary-amn",
    "transport-polar-amn",
)
EXTENDED_SCAN_VARIANTS = (
    "zero-iter",
    "identity-amn",
    "identity-amn-zero-iter",
    "phase-amn",
    "phase-amn-zero-iter",
    "smooth-amn",
    "smooth-amn-zero-iter",
    "neighbor-smooth-amn",
    "neighbor-smooth-amn-zero-iter",
    "global-smooth-amn",
    "global-smooth-amn-zero-iter",
    "unitary-smooth-amn",
    "unitary-smooth-amn-zero-iter",
    "transport-unitary-amn",
    "transport-unitary-amn-zero-iter",
    "transport-polar-amn",
    "transport-polar-amn-zero-iter",
    "scdm-polar-amn",
    "scdm-polar-amn-zero-iter",
    "qr-pivot-amn",
    "qr-pivot-amn-zero-iter",
)
CP2K_INITIAL_PROJECTIONS = (
    "IDENTITY",
    "AO_FIRST",
    "AO_RANKED",
    "AO_SCDM",
    "AO_FIRST_ORTHO",
    "AO_RANKED_ORTHO",
    "AO_SCDM_ORTHO",
    "AO_RANKED_RAW",
    "AO_SCDM_RAW",
    "AO_PAIR_SCDM_RAW",
    "AO_PAIR_SIGNED_SCDM_RAW",
    "CENTER_SCDM_RAW",
    "AUTO_SCDM_RAW",
)
CP2K_AUTO_PROJECTION_SCORES = (
    "CONDITIONING",
    "SPREAD_PROXY",
    "BALANCED",
    "MMN_PROXY",
    "MLWF_PROXY",
)
CP2K_PROJECTION_SEED_SUFFIXES = {
    "IDENTITY": "BID",
    "AO_FIRST": "AOF",
    "AO_RANKED": "AOR",
    "AO_SCDM": "AOS",
    "AO_FIRST_ORTHO": "AOFO",
    "AO_RANKED_ORTHO": "AORO",
    "AO_SCDM_ORTHO": "AOSO",
    "AO_RANKED_RAW": "AORW",
    "AO_SCDM_RAW": "AOSW",
    "AO_PAIR_SCDM_RAW": "AOPW",
    "AO_PAIR_SIGNED_SCDM_RAW": "AOPS",
    "CENTER_SCDM_RAW": "CENW",
    "AUTO_SCDM_RAW": "AUTW",
}
CP2K_AUTO_PROJECTION_SCORE_SEED_SUFFIXES = {
    "CONDITIONING": "CON",
    "SPREAD_PROXY": "SPR",
    "BALANCED": "BAL",
    "MMN_PROXY": "MMN",
    "MLWF_PROXY": "MLW",
}
CP2K_PROJECTION_SMOOTHINGS = (
    "NONE",
    "SEQUENTIAL",
    "NEIGHBOR",
    "UNITARY_NEIGHBOR",
    "MMN_TRANSPORT",
)
CP2K_PROJECTION_SMOOTHING_SEED_SUFFIXES = {
    "NONE": "N",
    "SEQUENTIAL": "S",
    "NEIGHBOR": "G",
    "UNITARY_NEIGHBOR": "U",
    "MMN_TRANSPORT": "M",
}
DEFAULT_CP2K_PROJECTION_MATRIX_PROJECTIONS = (
    "AO_SCDM_RAW",
    "AO_PAIR_SCDM_RAW",
    "AO_PAIR_SIGNED_SCDM_RAW",
    "CENTER_SCDM_RAW",
    "AUTO_SCDM_RAW",
)
MAX_WANNIER90_SEED_LENGTH = 48
STRATEGY_MATRIX_HEADER = (
    "case",
    "seed",
    "status",
    "selected_label",
    "selected_classification",
    "selected_best_iteration",
    "selected_best_omega_total",
    "selected_final_omega_total",
    "selected_final_minus_best",
    "recommended_action",
    "num_candidates",
    "num_failed_candidates",
    "cap_label",
    "cap_minus_target",
    "cap_omega_d_delta",
    "cap_omega_od_delta",
    "cap_spread_delta",
    "cap_center_delta",
    "cap_classification",
    "case_workdir",
    "error",
)
STRATEGY_MATRIX_CAP_TOL_FIELDS = (
    "cap_minus_target",
    "cap_omega_d_delta",
    "cap_omega_od_delta",
    "cap_spread_delta",
    "cap_center_delta",
)


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


@dataclass(frozen=True)
class WoutIteration:
    iteration: int
    delta_spread: float
    rms_gradient: float
    omega_total: float


@dataclass(frozen=True)
class WoutState:
    iteration: int
    functions: tuple[WannierFunction, ...]
    omega_d: float
    omega_od: float
    omega_total: float


@dataclass(frozen=True)
class WoutHistorySummary:
    classification: str
    n_iterations: int
    start: WoutIteration
    best: WoutIteration
    final: WoutIteration
    final_omega_total: float
    final_minus_best: float
    first_uphill_iteration: int | None
    max_step_increase: float


@dataclass(frozen=True)
class WoutHistoryRecord:
    spec: RunSpec
    summary: WoutSummary
    history: WoutHistorySummary


@dataclass(frozen=True)
class StrategyCapRun:
    source_record: WoutHistoryRecord
    cap_spec: RunSpec


@dataclass(frozen=True)
class StrategyCapComparison:
    source_label: str
    cap_label: str
    target_num_iter: int
    target_best_omega_total: float
    target_omega_d: float
    target_omega_od: float
    cap_final_omega_total: float
    cap_omega_d: float
    cap_omega_od: float
    cap_minus_target: float
    omega_d_delta: float
    omega_od_delta: float
    spread_delta: float
    center_delta: float
    cap_classification: str


@dataclass(frozen=True)
class StrategyCase:
    label: str
    seed: str
    source_workdir: Path


@dataclass(frozen=True)
class CP2KStrategyCase:
    label: str
    input_path: Path
    seed: str | None


@dataclass(frozen=True)
class W90Cell:
    vectors: tuple[tuple[float, float, float], ...]
    inverse: tuple[tuple[float, float, float], ...]


@dataclass(frozen=True)
class ComparisonSummary:
    omega_total_delta: float
    spread_delta: float
    center_delta: float
    assignment: tuple[int, ...]


@dataclass(frozen=True)
class AmnSummary:
    seed: str
    path: Path
    num_bands: int
    num_kpts: int
    num_wann: int
    max_identity_deviation: float
    max_column_norm_delta: float
    max_offdiag: float
    nonzero_fraction: float
    is_identity: bool


@dataclass(frozen=True)
class AmnMatrix:
    seed: str
    path: Path
    num_bands: int
    num_kpts: int
    num_wann: int
    matrices: tuple[tuple[tuple[complex, ...], ...], ...]


@dataclass(frozen=True)
class AmnComparison:
    max_abs_delta: float
    rms_abs_delta: float
    max_singular_delta: float
    rms_singular_delta: float
    max_subspace_deviation: float
    rms_subspace_deviation: float
    min_subspace_svalue: float
    max_subspace_kpt_index: int
    max_kpt_delta: float
    max_kpt_index: int


@dataclass(frozen=True)
class EigSummary:
    seed: str
    path: Path
    num_bands: int
    num_kpts: int
    values: tuple[tuple[float, ...], ...]


@dataclass(frozen=True)
class MmnBlock:
    ikpt: int
    jkpt: int
    cell: tuple[int, int, int]
    values: tuple[tuple[complex, ...], ...]


@dataclass(frozen=True)
class MmnSummary:
    seed: str
    path: Path
    num_bands: int
    num_kpts: int
    num_neighbors: int
    blocks: tuple[MmnBlock, ...]


@dataclass(frozen=True)
class ExportComparison:
    eig_max_abs_delta: float
    mmn_max_abs_delta: float
    mmn_rms_abs_delta: float
    mmn_max_singular_delta: float
    mmn_rms_singular_delta: float
    mmn_max_block_delta: float
    mmn_max_block_index: int


@dataclass(frozen=True)
class PreparedVariant:
    seed: str
    variant: str
    workdir: Path
    removed_amn: bool


@dataclass(frozen=True)
class RunSpec:
    label: str
    seed: str
    workdir: Path


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


def parse_wout_iterations(seed: str, workdir: Path) -> tuple[WoutIteration, ...]:
    path = workdir / f"{seed}.wout"
    iterations: list[WoutIteration] = []
    in_iteration_table = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if (
            "Delta Spread" in line
            and "RMS Gradient" in line
            and "Spread (Bohr^2)" in line
        ):
            in_iteration_table = True
            continue
        if not in_iteration_table:
            continue
        iteration = _parse_wout_iteration_line(line)
        if iteration is None:
            continue
        iterations.append(
            WoutIteration(
                iteration=iteration[0],
                delta_spread=iteration[1],
                rms_gradient=iteration[2],
                omega_total=iteration[3],
            )
        )
    return tuple(iterations)


def _parse_wout_iteration_line(line: str) -> tuple[int, float, float, float] | None:
    fields = line.replace("<-- CONV", "").split()
    if len(fields) < 4 or not fields[0].isdigit():
        return None
    try:
        return int(fields[0]), float(fields[1]), float(fields[2]), float(fields[3])
    except ValueError:
        return None


def parse_wout_states(seed: str, workdir: Path) -> tuple[WoutState, ...]:
    path = workdir / f"{seed}.wout"
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    states: list[WoutState] = []
    functions: list[WannierFunction] = []
    in_wannierise = False

    for iline, line in enumerate(lines):
        if (
            "Delta Spread" in line
            and "RMS Gradient" in line
            and "Spread (Bohr^2)" in line
        ):
            in_wannierise = True
            continue
        if not in_wannierise:
            continue

        match = WF_RE.search(line)
        if match is not None:
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
            continue

        iteration = _parse_wout_iteration_line(line)
        if iteration is None or not functions:
            continue
        sprd_match = None
        for follow_line in lines[iline + 1 : iline + 4]:
            sprd_match = SPRD_RE.search(follow_line)
            if sprd_match is not None:
                break
        if sprd_match is None:
            functions = []
            continue
        states.append(
            WoutState(
                iteration=iteration[0],
                functions=tuple(functions),
                omega_d=float(sprd_match.group(1)),
                omega_od=float(sprd_match.group(2)),
                omega_total=float(sprd_match.group(3)),
            )
        )
        functions = []
    return tuple(states)


def wout_state_as_summary(seed: str, path: Path, state: WoutState) -> WoutSummary:
    return WoutSummary(
        seed=seed,
        path=path,
        functions=state.functions,
        omega_i=float("nan"),
        omega_d=state.omega_d,
        omega_od=state.omega_od,
        omega_total=state.omega_total,
        all_done=True,
    )


def summarize_wout_history(
    iterations: Sequence[WoutIteration],
    final_omega_total: float,
    *,
    drift_tol: float,
    uphill_tol: float,
) -> WoutHistorySummary:
    if not iterations:
        raise ValueError("missing Wannier90 iteration table")

    start = iterations[0]
    best = min(iterations, key=lambda iteration: iteration.omega_total)
    first_uphill: int | None = None
    max_step_increase = 0.0
    previous = start
    for iteration in iterations[1:]:
        step_increase = iteration.omega_total - previous.omega_total
        max_step_increase = max(max_step_increase, step_increase)
        if first_uphill is None and step_increase > uphill_tol:
            first_uphill = iteration.iteration
        previous = iteration

    final = iterations[-1]
    final_minus_best = final_omega_total - best.omega_total
    if len(iterations) == 1:
        classification = "zero_iter"
    elif final_minus_best > drift_tol and best.iteration < final.iteration:
        classification = "late_drift"
    else:
        classification = "stable"

    return WoutHistorySummary(
        classification=classification,
        n_iterations=len(iterations),
        start=start,
        best=best,
        final=final,
        final_omega_total=final_omega_total,
        final_minus_best=final_minus_best,
        first_uphill_iteration=first_uphill,
        max_step_increase=max_step_increase,
    )


def read_wout_history_record(
    spec: RunSpec, *, drift_tol: float, uphill_tol: float
) -> WoutHistoryRecord:
    summary = parse_wout(spec.seed, spec.workdir)
    history = summarize_wout_history(
        parse_wout_iterations(spec.seed, spec.workdir),
        summary.omega_total,
        drift_tol=drift_tol,
        uphill_tol=uphill_tol,
    )
    return WoutHistoryRecord(spec=spec, summary=summary, history=history)


def _parse_float_triplet(line: str) -> tuple[float, float, float]:
    fields = line.split()
    if len(fields) < 3:
        raise ValueError(f"expected three floating-point values, got: {line!r}")
    return float(fields[0]), float(fields[1]), float(fields[2])


def _det3(matrix: tuple[tuple[float, float, float], ...]) -> float:
    a11, a12, a13 = matrix[0]
    a21, a22, a23 = matrix[1]
    a31, a32, a33 = matrix[2]
    return (
        a11 * (a22 * a33 - a23 * a32)
        - a12 * (a21 * a33 - a23 * a31)
        + a13 * (a21 * a32 - a22 * a31)
    )


def _invert3(
    matrix: tuple[tuple[float, float, float], ...],
) -> tuple[tuple[float, float, float], ...]:
    det = _det3(matrix)
    if abs(det) < 1.0e-14:
        raise ValueError("unit cell is singular")
    a11, a12, a13 = matrix[0]
    a21, a22, a23 = matrix[1]
    a31, a32, a33 = matrix[2]
    return (
        (
            (a22 * a33 - a23 * a32) / det,
            (a13 * a32 - a12 * a33) / det,
            (a12 * a23 - a13 * a22) / det,
        ),
        (
            (a23 * a31 - a21 * a33) / det,
            (a11 * a33 - a13 * a31) / det,
            (a13 * a21 - a11 * a23) / det,
        ),
        (
            (a21 * a32 - a22 * a31) / det,
            (a12 * a31 - a11 * a32) / det,
            (a11 * a22 - a12 * a21) / det,
        ),
    )


def _matvec(
    matrix: tuple[tuple[float, float, float], ...], vector: tuple[float, float, float]
) -> tuple[float, float, float]:
    return (
        sum(matrix[0][j] * vector[j] for j in range(3)),
        sum(matrix[1][j] * vector[j] for j in range(3)),
        sum(matrix[2][j] * vector[j] for j in range(3)),
    )


def parse_win_cell(seed: str, workdir: Path) -> W90Cell:
    path = workdir / f"{seed}.win"
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = -1
    for iline, line in enumerate(lines):
        if line.strip().lower() == "begin unit_cell_cart":
            start = iline + 1
            break
    if start < 0:
        raise ValueError(f"{path}: missing unit_cell_cart block")

    vectors: list[tuple[float, float, float]] = []
    iline = start
    while iline < len(lines) and len(vectors) < 3:
        line = lines[iline].strip()
        iline += 1
        if not line:
            continue
        if line.lower() in {"bohr", "ang", "angstrom"}:
            continue
        if line.lower() == "end unit_cell_cart":
            break
        vectors.append(_parse_float_triplet(line))
    if len(vectors) != 3:
        raise ValueError(f"{path}: incomplete unit_cell_cart block")
    cell_vectors = tuple(vectors)
    return W90Cell(vectors=cell_vectors, inverse=_invert3(cell_vectors))


def parse_amn_matrix(seed: str, workdir: Path) -> AmnMatrix:
    path = workdir / f"{seed}.amn"
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header_index = -1
    num_bands = 0
    num_kpts = 0
    num_wann = 0
    for iline, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("!"):
            continue
        fields = stripped.split()
        if len(fields) >= 3:
            num_bands = int(fields[0])
            num_kpts = int(fields[1])
            num_wann = int(fields[2])
            header_index = iline
            break
    if header_index < 0:
        raise ValueError(f"{path}: missing amn dimensions")
    if num_bands <= 0 or num_kpts <= 0 or num_wann <= 0:
        raise ValueError(f"{path}: invalid amn dimensions")

    matrices = [
        [[0.0j for _ in range(num_wann)] for _ in range(num_bands)]
        for _ in range(num_kpts)
    ]
    entries = 0
    for line in lines[header_index + 1 :]:
        stripped = line.strip()
        if not stripped or stripped.startswith("!"):
            continue
        fields = stripped.split()
        if len(fields) < 5:
            raise ValueError(f"{path}: malformed amn entry: {line!r}")
        iband = int(fields[0])
        iwann = int(fields[1])
        ikpt = int(fields[2])
        if not (
            1 <= iband <= num_bands and 1 <= iwann <= num_wann and 1 <= ikpt <= num_kpts
        ):
            raise ValueError(f"{path}: amn index out of range: {line!r}")
        value = complex(float(fields[3]), float(fields[4]))
        matrices[ikpt - 1][iband - 1][iwann - 1] = value
        entries += 1
    expected_entries = num_bands * num_kpts * num_wann
    if entries != expected_entries:
        raise ValueError(
            f"{path}: expected {expected_entries} amn entries, found {entries}"
        )
    return AmnMatrix(
        seed=seed,
        path=path,
        num_bands=num_bands,
        num_kpts=num_kpts,
        num_wann=num_wann,
        matrices=tuple(tuple(tuple(row) for row in matrix) for matrix in matrices),
    )


def parse_amn(seed: str, workdir: Path, identity_tol: float) -> AmnSummary:
    amn = parse_amn_matrix(seed, workdir)
    max_identity_deviation = 0.0
    max_column_norm_delta = 0.0
    max_offdiag = 0.0
    nonzero = 0
    expected_entries = amn.num_bands * amn.num_kpts * amn.num_wann
    for matrix in amn.matrices:
        for iband in range(amn.num_bands):
            for iwann in range(amn.num_wann):
                if abs(matrix[iband][iwann]) > identity_tol:
                    nonzero += 1
        for iwann in range(amn.num_wann):
            column_norm = sum(
                abs(matrix[iband][iwann]) ** 2 for iband in range(amn.num_bands)
            )
            max_column_norm_delta = max(max_column_norm_delta, abs(column_norm - 1.0))
            for jwann in range(amn.num_wann):
                overlap = sum(
                    matrix[iband][iwann].conjugate() * matrix[iband][jwann]
                    for iband in range(amn.num_bands)
                )
                target = 1.0 if iwann == jwann else 0.0
                deviation = abs(overlap - target)
                max_identity_deviation = max(max_identity_deviation, deviation)
                if iwann != jwann:
                    max_offdiag = max(max_offdiag, abs(overlap))

    nonzero_fraction = float(nonzero) / float(expected_entries)
    is_identity = (
        amn.num_bands == amn.num_wann
        and max_identity_deviation <= identity_tol
        and abs(nonzero_fraction - 1.0 / float(amn.num_bands)) <= identity_tol
    )
    return AmnSummary(
        seed=amn.seed,
        path=amn.path,
        num_bands=amn.num_bands,
        num_kpts=amn.num_kpts,
        num_wann=amn.num_wann,
        max_identity_deviation=max_identity_deviation,
        max_column_norm_delta=max_column_norm_delta,
        max_offdiag=max_offdiag,
        nonzero_fraction=nonzero_fraction,
        is_identity=is_identity,
    )


def write_amn_matrix(amn: AmnMatrix, path: Path, comment: str) -> None:
    lines = [f"! {comment}", f"{amn.num_bands:8d}{amn.num_kpts:8d}{amn.num_wann:8d}"]
    for ikpt, matrix in enumerate(amn.matrices, start=1):
        for iwann in range(amn.num_wann):
            for iband in range(amn.num_bands):
                value = matrix[iband][iwann]
                lines.append(
                    f"{iband + 1:8d}{iwann + 1:8d}{ikpt:8d}"
                    f"{value.real:30.14E}{value.imag:30.14E}"
                )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _phase_fixed_columns(
    matrix: Sequence[Sequence[complex]], tolerance: float = 1.0e-14
) -> tuple[tuple[complex, ...], ...]:
    nrows = len(matrix)
    if nrows == 0:
        return ()
    ncols = len(matrix[0])
    if any(len(row) != ncols for row in matrix):
        raise ValueError("phase-fixing input must be rectangular")

    columns = [[matrix[irow][icol] for irow in range(nrows)] for icol in range(ncols)]
    fixed_columns: list[list[complex]] = []
    for column in columns:
        pivot = max(column, key=abs)
        if abs(pivot) > tolerance:
            phase = pivot.conjugate() / abs(pivot)
            fixed_columns.append([phase * value for value in column])
        else:
            fixed_columns.append(list(column))
    return tuple(
        tuple(fixed_columns[icol][irow] for icol in range(ncols))
        for irow in range(nrows)
    )


def _column_overlap(
    left_matrix: Sequence[Sequence[complex]],
    right_matrix: Sequence[Sequence[complex]],
    left_column: int,
    right_column: int,
) -> complex:
    return sum(
        left_matrix[irow][left_column].conjugate() * right_matrix[irow][right_column]
        for irow in range(len(left_matrix))
    )


def _align_amn_matrix_to_reference(
    reference: Sequence[Sequence[complex]], matrix: Sequence[Sequence[complex]]
) -> tuple[tuple[complex, ...], ...]:
    nrows = len(matrix)
    if nrows == 0:
        return ()
    ncols = len(matrix[0])
    if any(len(row) != ncols for row in matrix):
        raise ValueError("AMN gauge alignment input must be rectangular")
    if len(reference) != nrows or any(len(row) != ncols for row in reference):
        raise ValueError("AMN gauge alignment inputs must have matching dimensions")

    cost = [
        [
            -abs(_column_overlap(reference, matrix, left_column, right_column))
            for right_column in range(ncols)
        ]
        for left_column in range(ncols)
    ]
    assignment = solve_assignment(cost)
    columns: list[list[complex]] = []
    for left_column, right_column in enumerate(assignment):
        column = [matrix[irow][right_column] for irow in range(nrows)]
        overlap = _column_overlap(reference, matrix, left_column, right_column)
        if abs(overlap) > 1.0e-14:
            phase = overlap.conjugate() / abs(overlap)
            column = [phase * value for value in column]
        columns.append(column)
    return tuple(
        tuple(columns[icol][irow] for icol in range(ncols)) for irow in range(nrows)
    )


def _require_numpy() -> Any:
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError(
            "This Wannier90 gauge diagnostic requires numpy. "
            "Install numpy or use one of the pure-Python variants."
        ) from exc
    return np


def _unitary_align_amn_matrix_to_reference(
    reference: Sequence[Sequence[complex]], matrix: Sequence[Sequence[complex]]
) -> tuple[tuple[complex, ...], ...]:
    nrows = len(matrix)
    if nrows == 0:
        return ()
    ncols = len(matrix[0])
    if any(len(row) != ncols for row in matrix):
        raise ValueError("unitary AMN gauge alignment input must be rectangular")
    if len(reference) != nrows or any(len(row) != ncols for row in reference):
        raise ValueError(
            "unitary AMN gauge alignment inputs must have matching dimensions"
        )

    np = _require_numpy()
    source = np.array(matrix, dtype=np.complex128)
    target = np.array(reference, dtype=np.complex128)
    left, _, right_h = np.linalg.svd(source.conj().T @ target, full_matrices=False)
    aligned = source @ (left @ right_h)
    return tuple(tuple(complex(value) for value in row) for row in aligned.tolist())


def _polar_amn_matrix(
    matrix: Sequence[Sequence[complex]],
) -> tuple[tuple[complex, ...], ...]:
    np = _require_numpy()
    source = np.array(matrix, dtype=np.complex128)
    left, _, right_h = np.linalg.svd(source, full_matrices=False)
    polar = left @ right_h
    return tuple(tuple(complex(value) for value in row) for row in polar.tolist())


def _amn_matrix_delta(
    left: Sequence[Sequence[complex]], right: Sequence[Sequence[complex]]
) -> float:
    if len(left) != len(right):
        raise ValueError("AMN matrix row counts differ")
    max_delta = 0.0
    for left_row, right_row in zip(left, right):
        if len(left_row) != len(right_row):
            raise ValueError("AMN matrix column counts differ")
        for left_value, right_value in zip(left_row, right_row):
            max_delta = max(max_delta, abs(left_value - right_value))
    return max_delta


def _average_amn_matrices(
    matrices: Sequence[Sequence[Sequence[complex]]],
) -> tuple[tuple[complex, ...], ...]:
    if not matrices:
        raise ValueError("cannot average an empty AMN matrix list")
    nrows = len(matrices[0])
    if nrows == 0:
        return ()
    ncols = len(matrices[0][0])
    if any(
        len(matrix) != nrows or any(len(row) != ncols for row in matrix)
        for matrix in matrices
    ):
        raise ValueError("AMN matrix average inputs must have matching dimensions")
    scale = 1.0 / float(len(matrices))
    return tuple(
        tuple(
            scale * sum(matrix[irow][icol] for matrix in matrices)
            for icol in range(ncols)
        )
        for irow in range(nrows)
    )


def _build_mmn_neighbor_graph(num_kpts: int, mmn: MmnSummary) -> list[list[int]]:
    neighbors: list[list[int]] = [[] for _ in range(num_kpts)]
    for block in mmn.blocks:
        ileft = block.ikpt - 1
        iright = block.jkpt - 1
        if not (0 <= ileft < num_kpts) or not (0 <= iright < num_kpts):
            raise ValueError(
                f"mmn neighbor block references invalid k-point pair "
                f"{block.ikpt}, {block.jkpt}"
            )
        if iright not in neighbors[ileft]:
            neighbors[ileft].append(iright)
        if ileft not in neighbors[iright]:
            neighbors[iright].append(ileft)
    for node_neighbors in neighbors:
        node_neighbors.sort()
    return neighbors


def _mmn_outgoing_blocks(num_kpts: int, mmn: MmnSummary) -> list[list[MmnBlock]]:
    outgoing: list[list[MmnBlock]] = [[] for _ in range(num_kpts)]
    for block in mmn.blocks:
        ikpt = block.ikpt - 1
        jkpt = block.jkpt - 1
        if not (0 <= ikpt < num_kpts) or not (0 <= jkpt < num_kpts):
            raise ValueError(
                f"mmn neighbor block references invalid k-point pair "
                f"{block.ikpt}, {block.jkpt}"
            )
        outgoing[ikpt].append(block)
    return outgoing


def _transport_amn_matrix(
    block: MmnBlock, matrix: Sequence[Sequence[complex]]
) -> tuple[tuple[complex, ...], ...]:
    nrows = len(matrix)
    if nrows == 0:
        return ()
    ncols = len(matrix[0])
    if any(len(row) != ncols for row in matrix):
        raise ValueError("AMN transport input must be rectangular")
    if len(block.values) != nrows or any(len(row) != nrows for row in block.values):
        raise ValueError("MMN transport block and AMN matrix dimensions differ")

    return tuple(
        tuple(
            sum(block.values[irow][jrow] * matrix[jrow][icol] for jrow in range(nrows))
            for icol in range(ncols)
        )
        for irow in range(nrows)
    )


def _smooth_amn_columns(amn: AmnMatrix) -> AmnMatrix:
    smoothed: list[tuple[tuple[complex, ...], ...]] = []
    previous = _phase_fixed_columns(amn.matrices[0])
    smoothed.append(previous)
    for matrix in amn.matrices[1:]:
        previous = _align_amn_matrix_to_reference(previous, matrix)
        smoothed.append(previous)
    return AmnMatrix(
        seed=amn.seed,
        path=amn.path,
        num_bands=amn.num_bands,
        num_kpts=amn.num_kpts,
        num_wann=amn.num_wann,
        matrices=tuple(smoothed),
    )


def _neighbor_smooth_amn_columns(amn: AmnMatrix, mmn: MmnSummary) -> AmnMatrix:
    if (amn.num_bands, amn.num_kpts) != (mmn.num_bands, mmn.num_kpts):
        raise ValueError(
            "amn/mmn dimensions differ: "
            f"{amn.num_bands}x{amn.num_kpts} != {mmn.num_bands}x{mmn.num_kpts}"
        )

    neighbors = _build_mmn_neighbor_graph(amn.num_kpts, mmn)
    smoothed: list[tuple[tuple[complex, ...], ...] | None] = [None] * amn.num_kpts
    smoothed[0] = _phase_fixed_columns(amn.matrices[0])
    queue = deque([0])
    while queue:
        parent = queue.popleft()
        parent_matrix = smoothed[parent]
        if parent_matrix is None:
            raise ValueError("internal neighbor smoothing error: missing parent")
        for child in neighbors[parent]:
            if smoothed[child] is not None:
                continue
            smoothed[child] = _align_amn_matrix_to_reference(
                parent_matrix, amn.matrices[child]
            )
            queue.append(child)

    for ikpt, matrix in enumerate(smoothed):
        if matrix is None:
            smoothed[ikpt] = _phase_fixed_columns(amn.matrices[ikpt])

    return AmnMatrix(
        seed=amn.seed,
        path=amn.path,
        num_bands=amn.num_bands,
        num_kpts=amn.num_kpts,
        num_wann=amn.num_wann,
        matrices=tuple(matrix for matrix in smoothed if matrix is not None),
    )


def _global_smooth_amn_columns(
    amn: AmnMatrix,
    mmn: MmnSummary,
    *,
    max_iterations: int = 50,
    tolerance: float = 1.0e-12,
) -> AmnMatrix:
    if (amn.num_bands, amn.num_kpts) != (mmn.num_bands, mmn.num_kpts):
        raise ValueError(
            "amn/mmn dimensions differ: "
            f"{amn.num_bands}x{amn.num_kpts} != {mmn.num_bands}x{mmn.num_kpts}"
        )

    neighbors = _build_mmn_neighbor_graph(amn.num_kpts, mmn)
    current = [tuple(tuple(row) for row in matrix) for matrix in amn.matrices]
    for _ in range(max_iterations):
        updated: list[tuple[tuple[complex, ...], ...]] = []
        max_delta = 0.0
        for ikpt, matrix in enumerate(amn.matrices):
            reference_matrices = [current[ikpt]]
            reference_matrices.extend(current[neighbor] for neighbor in neighbors[ikpt])
            reference = _average_amn_matrices(reference_matrices)
            aligned = _align_amn_matrix_to_reference(reference, matrix)
            updated.append(aligned)
            max_delta = max(max_delta, _amn_matrix_delta(current[ikpt], aligned))
        current = updated
        if max_delta <= tolerance:
            break

    return AmnMatrix(
        seed=amn.seed,
        path=amn.path,
        num_bands=amn.num_bands,
        num_kpts=amn.num_kpts,
        num_wann=amn.num_wann,
        matrices=tuple(current),
    )


def _unitary_smooth_amn_columns(
    amn: AmnMatrix,
    mmn: MmnSummary,
    *,
    max_iterations: int = 50,
    tolerance: float = 1.0e-12,
) -> AmnMatrix:
    if (amn.num_bands, amn.num_kpts) != (mmn.num_bands, mmn.num_kpts):
        raise ValueError(
            "amn/mmn dimensions differ: "
            f"{amn.num_bands}x{amn.num_kpts} != {mmn.num_bands}x{mmn.num_kpts}"
        )

    neighbors = _build_mmn_neighbor_graph(amn.num_kpts, mmn)
    current = [tuple(tuple(row) for row in matrix) for matrix in amn.matrices]
    for _ in range(max_iterations):
        updated: list[tuple[tuple[complex, ...], ...]] = []
        max_delta = 0.0
        for ikpt, matrix in enumerate(amn.matrices):
            reference_matrices = [current[ikpt]]
            reference_matrices.extend(current[neighbor] for neighbor in neighbors[ikpt])
            reference = _average_amn_matrices(reference_matrices)
            aligned = _unitary_align_amn_matrix_to_reference(reference, matrix)
            updated.append(aligned)
            max_delta = max(max_delta, _amn_matrix_delta(current[ikpt], aligned))
        current = updated
        if max_delta <= tolerance:
            break

    return AmnMatrix(
        seed=amn.seed,
        path=amn.path,
        num_bands=amn.num_bands,
        num_kpts=amn.num_kpts,
        num_wann=amn.num_wann,
        matrices=tuple(current),
    )


def _transport_smooth_amn_columns(
    amn: AmnMatrix,
    mmn: MmnSummary,
    *,
    mode: str,
    max_iterations: int = 50,
    tolerance: float = 1.0e-12,
) -> AmnMatrix:
    if (amn.num_bands, amn.num_kpts) != (mmn.num_bands, mmn.num_kpts):
        raise ValueError(
            "amn/mmn dimensions differ: "
            f"{amn.num_bands}x{amn.num_kpts} != {mmn.num_bands}x{mmn.num_kpts}"
        )
    if mode not in {"unitary", "polar"}:
        raise ValueError(f"unknown transport smoothing mode: {mode}")

    outgoing = _mmn_outgoing_blocks(amn.num_kpts, mmn)
    current = [tuple(tuple(row) for row in matrix) for matrix in amn.matrices]
    for _ in range(max_iterations):
        updated: list[tuple[tuple[complex, ...], ...]] = []
        max_delta = 0.0
        for ikpt, matrix in enumerate(amn.matrices):
            reference_matrices = [current[ikpt]]
            reference_matrices.extend(
                _transport_amn_matrix(block, current[block.jkpt - 1])
                for block in outgoing[ikpt]
            )
            reference = _average_amn_matrices(reference_matrices)
            if mode == "unitary":
                aligned = _unitary_align_amn_matrix_to_reference(reference, matrix)
            else:
                aligned = _polar_amn_matrix(reference)
            updated.append(aligned)
            max_delta = max(max_delta, _amn_matrix_delta(current[ikpt], aligned))
        current = updated
        if max_delta <= tolerance:
            break

    return AmnMatrix(
        seed=amn.seed,
        path=amn.path,
        num_bands=amn.num_bands,
        num_kpts=amn.num_kpts,
        num_wann=amn.num_wann,
        matrices=tuple(current),
    )


def _scdm_polar_amn_columns(amn: AmnMatrix) -> AmnMatrix:
    matrices = [_polar_amn_matrix(matrix) for matrix in amn.matrices]
    return AmnMatrix(
        seed=amn.seed,
        path=amn.path,
        num_bands=amn.num_bands,
        num_kpts=amn.num_kpts,
        num_wann=amn.num_wann,
        matrices=tuple(matrices),
    )


def _global_qr_pivot_order(
    amn: AmnMatrix, tolerance: float = 1.0e-14
) -> tuple[int, ...]:
    vectors: list[list[complex]] = []
    for iwann in range(amn.num_wann):
        vector: list[complex] = []
        for matrix in amn.matrices:
            vector.extend(matrix[iband][iwann] for iband in range(amn.num_bands))
        vectors.append(vector)

    residuals = [list(vector) for vector in vectors]
    used = [False] * amn.num_wann
    order: list[int] = []
    for _ in range(amn.num_wann):
        pivot = -1
        pivot_norm = -1.0
        for iwann, residual in enumerate(residuals):
            if used[iwann]:
                continue
            norm = sum(abs(value) ** 2 for value in residual)
            if norm > pivot_norm:
                pivot_norm = norm
                pivot = iwann
        if pivot < 0:
            raise ValueError("internal QR pivoting error: missing pivot")
        used[pivot] = True
        order.append(pivot)
        if pivot_norm <= tolerance:
            continue
        pivot_scale = 1.0 / math.sqrt(pivot_norm)
        q_column = [pivot_scale * value for value in residuals[pivot]]
        for iwann, residual in enumerate(residuals):
            if used[iwann]:
                continue
            overlap = sum(
                q_value.conjugate() * value
                for q_value, value in zip(q_column, residual)
            )
            residuals[iwann] = [
                value - overlap * q_value for q_value, value in zip(q_column, residual)
            ]
    return tuple(order)


def _reorder_amn_columns(amn: AmnMatrix, order: Sequence[int]) -> AmnMatrix:
    if len(order) != amn.num_wann or sorted(order) != list(range(amn.num_wann)):
        raise ValueError("AMN column order must be a permutation of Wannier columns")
    return AmnMatrix(
        seed=amn.seed,
        path=amn.path,
        num_bands=amn.num_bands,
        num_kpts=amn.num_kpts,
        num_wann=amn.num_wann,
        matrices=tuple(
            tuple(tuple(row[iwann] for iwann in order) for row in matrix)
            for matrix in amn.matrices
        ),
    )


def smooth_amn_file(seed: str, workdir: Path) -> None:
    amn_path = workdir / f"{seed}.amn"
    ensure_amn_file(seed, workdir)
    smoothed = _smooth_amn_columns(parse_amn_matrix(seed, workdir))
    write_amn_matrix(
        smoothed, amn_path, "Wannier90 smoothed projections generated by CP2K tooling"
    )


def neighbor_smooth_amn_file(seed: str, workdir: Path) -> None:
    amn_path = workdir / f"{seed}.amn"
    ensure_amn_file(seed, workdir)
    smoothed = _neighbor_smooth_amn_columns(
        parse_amn_matrix(seed, workdir), parse_mmn(seed, workdir)
    )
    write_amn_matrix(
        smoothed,
        amn_path,
        "Wannier90 neighbor-smoothed projections generated by CP2K tooling",
    )


def global_smooth_amn_file(seed: str, workdir: Path) -> None:
    amn_path = workdir / f"{seed}.amn"
    ensure_amn_file(seed, workdir)
    smoothed = _global_smooth_amn_columns(
        parse_amn_matrix(seed, workdir), parse_mmn(seed, workdir)
    )
    write_amn_matrix(
        smoothed,
        amn_path,
        "Wannier90 globally smoothed projections generated by CP2K tooling",
    )


def unitary_smooth_amn_file(seed: str, workdir: Path) -> None:
    amn_path = workdir / f"{seed}.amn"
    ensure_amn_file(seed, workdir)
    smoothed = _unitary_smooth_amn_columns(
        parse_amn_matrix(seed, workdir), parse_mmn(seed, workdir)
    )
    write_amn_matrix(
        smoothed,
        amn_path,
        "Wannier90 unitary-smoothed projections generated by CP2K tooling",
    )


def transport_unitary_amn_file(seed: str, workdir: Path) -> None:
    amn_path = workdir / f"{seed}.amn"
    ensure_amn_file(seed, workdir)
    smoothed = _transport_smooth_amn_columns(
        parse_amn_matrix(seed, workdir), parse_mmn(seed, workdir), mode="unitary"
    )
    write_amn_matrix(
        smoothed,
        amn_path,
        "Wannier90 transport-unitary projections generated by CP2K tooling",
    )


def transport_polar_amn_file(seed: str, workdir: Path) -> None:
    amn_path = workdir / f"{seed}.amn"
    ensure_amn_file(seed, workdir)
    smoothed = _transport_smooth_amn_columns(
        parse_amn_matrix(seed, workdir), parse_mmn(seed, workdir), mode="polar"
    )
    write_amn_matrix(
        smoothed,
        amn_path,
        "Wannier90 transport-polar projections generated by CP2K tooling",
    )


def scdm_polar_amn_file(seed: str, workdir: Path) -> None:
    amn_path = workdir / f"{seed}.amn"
    ensure_amn_file(seed, workdir)
    smoothed = _scdm_polar_amn_columns(parse_amn_matrix(seed, workdir))
    write_amn_matrix(
        smoothed, amn_path, "Wannier90 SCDM-polar projections generated by CP2K tooling"
    )


def qr_pivot_amn_file(seed: str, workdir: Path) -> None:
    amn_path = workdir / f"{seed}.amn"
    ensure_amn_file(seed, workdir)
    amn = parse_amn_matrix(seed, workdir)
    reordered = _reorder_amn_columns(amn, _global_qr_pivot_order(amn))
    write_amn_matrix(
        reordered,
        amn_path,
        "Wannier90 QR-pivoted projections generated by CP2K tooling",
    )


def phase_fix_amn_file(seed: str, workdir: Path) -> None:
    amn_path = workdir / f"{seed}.amn"
    ensure_amn_file(seed, workdir)
    amn = parse_amn_matrix(seed, workdir)
    phase_fixed = AmnMatrix(
        seed=amn.seed,
        path=amn.path,
        num_bands=amn.num_bands,
        num_kpts=amn.num_kpts,
        num_wann=amn.num_wann,
        matrices=tuple(_phase_fixed_columns(matrix) for matrix in amn.matrices),
    )
    write_amn_matrix(
        phase_fixed,
        amn_path,
        "Wannier90 phase-fixed projections generated by CP2K tooling",
    )


def parse_eig(seed: str, workdir: Path) -> EigSummary:
    path = workdir / f"{seed}.eig"
    entries: list[tuple[int, int, float]] = []
    num_bands = 0
    num_kpts = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("!"):
            continue
        fields = stripped.split()
        if len(fields) < 3:
            raise ValueError(f"{path}: malformed eig entry: {line!r}")
        iband = int(fields[0])
        ikpt = int(fields[1])
        value = float(fields[2])
        if iband <= 0 or ikpt <= 0:
            raise ValueError(f"{path}: eig index out of range: {line!r}")
        entries.append((iband, ikpt, value))
        num_bands = max(num_bands, iband)
        num_kpts = max(num_kpts, ikpt)
    if num_bands == 0 or num_kpts == 0:
        raise ValueError(f"{path}: missing eig entries")

    values = [[math.nan for _ in range(num_bands)] for _ in range(num_kpts)]
    for iband, ikpt, value in entries:
        if not math.isnan(values[ikpt - 1][iband - 1]):
            raise ValueError(f"{path}: duplicate eig entry for band {iband}, k {ikpt}")
        values[ikpt - 1][iband - 1] = value
    missing = [
        (ikpt + 1, iband + 1)
        for ikpt, row in enumerate(values)
        for iband, value in enumerate(row)
        if math.isnan(value)
    ]
    if missing:
        ikpt, iband = missing[0]
        raise ValueError(f"{path}: missing eig entry for band {iband}, k {ikpt}")
    return EigSummary(
        seed=seed,
        path=path,
        num_bands=num_bands,
        num_kpts=num_kpts,
        values=tuple(tuple(row) for row in values),
    )


def parse_mmn(seed: str, workdir: Path) -> MmnSummary:
    path = workdir / f"{seed}.mmn"
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header_index = -1
    num_bands = 0
    num_kpts = 0
    num_neighbors = 0
    for iline, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("!"):
            continue
        fields = stripped.split()
        if len(fields) >= 3:
            num_bands = int(fields[0])
            num_kpts = int(fields[1])
            num_neighbors = int(fields[2])
            header_index = iline
            break
    if header_index < 0:
        raise ValueError(f"{path}: missing mmn dimensions")
    if num_bands <= 0 or num_kpts <= 0 or num_neighbors <= 0:
        raise ValueError(f"{path}: invalid mmn dimensions")

    blocks: list[MmnBlock] = []
    iline = header_index + 1
    expected_blocks = num_kpts * num_neighbors
    for _ in range(expected_blocks):
        while iline < len(lines) and not lines[iline].strip():
            iline += 1
        if iline >= len(lines):
            raise ValueError(f"{path}: missing mmn block header")
        fields = lines[iline].split()
        iline += 1
        if len(fields) < 5:
            raise ValueError(
                f"{path}: malformed mmn block header: {lines[iline - 1]!r}"
            )
        ikpt = int(fields[0])
        jkpt = int(fields[1])
        cell = (int(fields[2]), int(fields[3]), int(fields[4]))
        if not (1 <= ikpt <= num_kpts and 1 <= jkpt <= num_kpts):
            raise ValueError(f"{path}: mmn k-point index out of range")
        matrix = [[0.0j for _ in range(num_bands)] for _ in range(num_bands)]
        for jband in range(num_bands):
            for iband in range(num_bands):
                if iline >= len(lines):
                    raise ValueError(f"{path}: truncated mmn block")
                fields = lines[iline].split()
                iline += 1
                if len(fields) < 2:
                    raise ValueError(
                        f"{path}: malformed mmn matrix entry: {lines[iline - 1]!r}"
                    )
                matrix[iband][jband] = complex(float(fields[0]), float(fields[1]))
        blocks.append(
            MmnBlock(
                ikpt=ikpt,
                jkpt=jkpt,
                cell=cell,
                values=tuple(tuple(row) for row in matrix),
            )
        )
    return MmnSummary(
        seed=seed,
        path=path,
        num_bands=num_bands,
        num_kpts=num_kpts,
        num_neighbors=num_neighbors,
        blocks=tuple(blocks),
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


def _win_keyword(line: str) -> str:
    stripped = line.strip()
    if not stripped or stripped.startswith("!") or stripped.startswith("#"):
        return ""
    key = re.split(r"\s+|=", stripped, maxsplit=1)[0].strip().lower()
    return key


def _win_optional_keyword_value(seed: str, workdir: Path, keyword: str) -> str | None:
    path = workdir / f"{seed}.win"
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if _win_keyword(line) != keyword.lower():
            continue
        stripped = line.split("!", maxsplit=1)[0].split("#", maxsplit=1)[0].strip()
        fields = [field for field in re.split(r"\s+|=", stripped) if field]
        if len(fields) >= 2:
            return fields[1]
    return None


def _win_keyword_value(seed: str, workdir: Path, keyword: str) -> str:
    value = _win_optional_keyword_value(seed, workdir, keyword)
    if value is not None:
        return value
    path = workdir / f"{seed}.win"
    raise ValueError(f"{path}: missing {keyword}")


def _win_integer_value(seed: str, workdir: Path, keyword: str) -> int:
    try:
        value = int(_win_keyword_value(seed, workdir, keyword))
    except ValueError as exc:
        raise ValueError(f"{workdir / f'{seed}.win'}: invalid {keyword}") from exc
    if value <= 0:
        raise ValueError(f"{workdir / f'{seed}.win'}: invalid {keyword}")
    return value


def identity_amn_file(seed: str, workdir: Path) -> None:
    eig = parse_eig(seed, workdir)
    num_wann = _win_integer_value(seed, workdir, "num_wann")
    if _win_optional_keyword_value(seed, workdir, "num_bands") is None:
        num_bands = eig.num_bands
    else:
        num_bands = _win_integer_value(seed, workdir, "num_bands")
    if num_bands != eig.num_bands:
        raise ValueError(
            f"{workdir / f'{seed}.eig'}: num_bands differs from {seed}.win"
        )
    if num_bands < num_wann:
        raise ValueError("identity AMN requires num_bands >= num_wann")
    matrices = tuple(
        tuple(
            tuple(1.0 + 0.0j if iband == iwann else 0.0j for iwann in range(num_wann))
            for iband in range(num_bands)
        )
        for _ in range(eig.num_kpts)
    )
    write_amn_matrix(
        AmnMatrix(
            seed=seed,
            path=workdir / f"{seed}.amn",
            num_bands=num_bands,
            num_kpts=eig.num_kpts,
            num_wann=num_wann,
            matrices=matrices,
        ),
        workdir / f"{seed}.amn",
        "Wannier90 identity projections generated by CP2K tooling",
    )


def ensure_amn_file(seed: str, workdir: Path) -> None:
    if (workdir / f"{seed}.amn").exists():
        return
    identity_amn_file(seed, workdir)


def parse_run_spec(spec: str) -> RunSpec:
    fields = spec.split(":", maxsplit=2)
    if len(fields) != 3:
        raise ValueError(
            "run specifications must use LABEL:SEED:WORKDIR, got " f"{spec!r}"
        )
    label, seed, workdir = fields
    if not label or not seed or not workdir:
        raise ValueError(f"incomplete run specification: {spec!r}")
    return RunSpec(label=label, seed=seed, workdir=Path(workdir))


def parse_strategy_case(spec: str) -> StrategyCase:
    fields = spec.split(":", maxsplit=2)
    if len(fields) != 3:
        raise ValueError(
            "case specifications must use LABEL:SEED:SOURCE_WORKDIR, got " f"{spec!r}"
        )
    label, seed, source_workdir = fields
    if not label or not seed or not source_workdir:
        raise ValueError(f"incomplete case specification: {spec!r}")
    return StrategyCase(label=label, seed=seed, source_workdir=Path(source_workdir))


def parse_cp2k_strategy_case(spec: str) -> CP2KStrategyCase:
    fields = spec.split(":", maxsplit=2)
    if len(fields) not in {2, 3}:
        raise ValueError(
            "CP2K case specifications must use LABEL:INPUT_PATH[:SEED], got "
            f"{spec!r}"
        )
    label = fields[0]
    input_path = fields[1]
    seed = fields[2] if len(fields) == 3 else None
    if not label or not input_path or seed == "":
        raise ValueError(f"incomplete CP2K case specification: {spec!r}")
    return CP2KStrategyCase(label=label, input_path=Path(input_path), seed=seed)


def set_win_num_iter(seed: str, workdir: Path, num_iter: int) -> None:
    win_path = workdir / f"{seed}.win"
    lines = win_path.read_text(encoding="utf-8", errors="replace").splitlines()
    filtered_lines = [line for line in lines if _win_keyword(line) != "num_iter"]
    filtered_lines.extend(["", f"num_iter = {num_iter}"])
    win_path.write_text("\n".join(filtered_lines) + "\n", encoding="utf-8")


def copy_workdir(source_workdir: Path, output_workdir: Path, *, force: bool) -> None:
    if output_workdir.exists():
        if not force:
            raise ValueError(f"{output_workdir}: output directory already exists")
        shutil.rmtree(output_workdir)
    shutil.copytree(source_workdir, output_workdir)


def prepare_variant(
    seed: str, source_workdir: Path, output_workdir: Path, variant: str, *, force: bool
) -> PreparedVariant:
    copy_workdir(source_workdir, output_workdir, force=force)

    win_path = output_workdir / f"{seed}.win"
    lines = win_path.read_text(encoding="utf-8", errors="replace").splitlines()
    filtered_lines = [
        line for line in lines if _win_keyword(line) not in CONTROLLED_WIN_KEYS
    ]
    removed_amn = False
    suffix: list[str] = []
    if variant == "identity":
        pass
    elif variant == "zero-iter":
        ensure_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "identity-amn":
        identity_amn_file(seed, output_workdir)
    elif variant == "identity-amn-zero-iter":
        identity_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "phase-amn":
        phase_fix_amn_file(seed, output_workdir)
    elif variant == "phase-amn-zero-iter":
        phase_fix_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "smooth-amn":
        smooth_amn_file(seed, output_workdir)
    elif variant == "smooth-amn-zero-iter":
        smooth_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "neighbor-smooth-amn":
        neighbor_smooth_amn_file(seed, output_workdir)
    elif variant == "neighbor-smooth-amn-zero-iter":
        neighbor_smooth_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "global-smooth-amn":
        global_smooth_amn_file(seed, output_workdir)
    elif variant == "global-smooth-amn-zero-iter":
        global_smooth_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "unitary-smooth-amn":
        unitary_smooth_amn_file(seed, output_workdir)
    elif variant == "unitary-smooth-amn-zero-iter":
        unitary_smooth_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "transport-unitary-amn":
        transport_unitary_amn_file(seed, output_workdir)
    elif variant == "transport-unitary-amn-zero-iter":
        transport_unitary_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "transport-polar-amn":
        transport_polar_amn_file(seed, output_workdir)
    elif variant == "transport-polar-amn-zero-iter":
        transport_polar_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "transport-unitary-guiding-centres":
        transport_unitary_amn_file(seed, output_workdir)
        suffix.append("guiding_centres = true")
    elif variant == "transport-polar-guiding-centres":
        transport_polar_amn_file(seed, output_workdir)
        suffix.append("guiding_centres = true")
    elif variant == "scdm-polar-amn":
        scdm_polar_amn_file(seed, output_workdir)
    elif variant == "scdm-polar-amn-zero-iter":
        scdm_polar_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "qr-pivot-amn":
        qr_pivot_amn_file(seed, output_workdir)
    elif variant == "qr-pivot-amn-zero-iter":
        qr_pivot_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    elif variant == "native-auto-projections":
        amn_path = output_workdir / f"{seed}.amn"
        if amn_path.exists():
            amn_path.unlink()
            removed_amn = True
        suffix.append("auto_projections = true")
    elif variant == "native-auto-projections-zero-iter":
        amn_path = output_workdir / f"{seed}.amn"
        if amn_path.exists():
            amn_path.unlink()
            removed_amn = True
        suffix.extend(["auto_projections = true", "num_iter = 0"])
    elif variant == "external-bloch":
        amn_path = output_workdir / f"{seed}.amn"
        if amn_path.exists():
            amn_path.unlink()
            removed_amn = True
        suffix.append("use_bloch_phases = true")
    elif variant == "external-bloch-zero-iter":
        amn_path = output_workdir / f"{seed}.amn"
        if amn_path.exists():
            amn_path.unlink()
            removed_amn = True
        suffix.extend(["use_bloch_phases = true", "num_iter = 0"])
    elif variant == "guiding-centres":
        suffix.append("guiding_centres = true")
    elif variant == "guiding-centres-zero-iter":
        suffix.extend(["guiding_centres = true", "num_iter = 0"])
    elif variant == "identity-zero-iter":
        ensure_amn_file(seed, output_workdir)
        suffix.append("num_iter = 0")
    else:
        raise ValueError(f"unknown Wannier90 variant: {variant}")

    if suffix:
        filtered_lines.extend(["", *suffix])
    win_path.write_text("\n".join(filtered_lines) + "\n", encoding="utf-8")
    return PreparedVariant(
        seed=seed, variant=variant, workdir=output_workdir, removed_amn=removed_amn
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


def center_delta(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    cell: W90Cell | None,
) -> tuple[float, float, float]:
    delta = (left[0] - right[0], left[1] - right[1], left[2] - right[2])
    if cell is None:
        return delta
    frac_delta = _matvec(cell.inverse, delta)
    frac_delta = (
        frac_delta[0] - round(frac_delta[0]),
        frac_delta[1] - round(frac_delta[1]),
        frac_delta[2] - round(frac_delta[2]),
    )
    return _matvec(cell.vectors, frac_delta)


def center_distance(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    cell: W90Cell | None,
) -> float:
    delta = center_delta(left, right, cell)
    return math.sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2])


def solve_assignment(cost: Sequence[Sequence[float]]) -> tuple[int, ...]:
    nrows = len(cost)
    if nrows == 0:
        return ()
    ncols = len(cost[0])
    if ncols != nrows or any(len(row) != ncols for row in cost):
        raise ValueError("assignment cost matrix must be square")

    u = [0.0] * (nrows + 1)
    v = [0.0] * (ncols + 1)
    p = [0] * (ncols + 1)
    way = [0] * (ncols + 1)
    for irow in range(1, nrows + 1):
        p[0] = irow
        j0 = 0
        minv = [float("inf")] * (ncols + 1)
        used = [False] * (ncols + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = float("inf")
            j1 = 0
            for jcol in range(1, ncols + 1):
                if used[jcol]:
                    continue
                cur = cost[i0 - 1][jcol - 1] - u[i0] - v[jcol]
                if cur < minv[jcol]:
                    minv[jcol] = cur
                    way[jcol] = j0
                if minv[jcol] < delta:
                    delta = minv[jcol]
                    j1 = jcol
            for jcol in range(0, ncols + 1):
                if used[jcol]:
                    u[p[jcol]] += delta
                    v[jcol] -= delta
                else:
                    minv[jcol] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break

    assignment = [0] * nrows
    for jcol in range(1, ncols + 1):
        if p[jcol] > 0:
            assignment[p[jcol] - 1] = jcol - 1
    return tuple(assignment)


def compare_summaries(
    left: WoutSummary,
    right: WoutSummary,
    *,
    match_centers: bool,
    periodic_cell: W90Cell | None,
) -> ComparisonSummary:
    if len(left.functions) != len(right.functions):
        raise ValueError(
            f"number of Wannier functions differs: "
            f"{len(left.functions)} != {len(right.functions)}"
        )
    if match_centers:
        cost = [
            [
                center_distance(
                    left_function.center, right_function.center, periodic_cell
                )
                + abs(left_function.spread - right_function.spread)
                for right_function in right.functions
            ]
            for left_function in left.functions
        ]
        assignment = solve_assignment(cost)
    else:
        assignment = tuple(range(len(left.functions)))

    spread_delta = 0.0
    center_delta_value = 0.0
    for ileft, iright in enumerate(assignment):
        left_function = left.functions[ileft]
        right_function = right.functions[iright]
        spread_delta = max(
            spread_delta, abs(left_function.spread - right_function.spread)
        )
        if match_centers or periodic_cell is not None:
            center_delta_value = max(
                center_delta_value,
                center_distance(
                    left_function.center, right_function.center, periodic_cell
                ),
            )
        else:
            center_delta_value = max(
                center_delta_value,
                max_abs_delta(left_function.center, right_function.center),
            )

    omega_delta = abs(left.omega_total - right.omega_total)
    return ComparisonSummary(
        omega_total_delta=omega_delta,
        spread_delta=spread_delta,
        center_delta=center_delta_value,
        assignment=assignment,
    )


def _jacobi_symmetric_eigenvalues(
    matrix: Sequence[Sequence[float]],
) -> tuple[float, ...]:
    values = [list(row) for row in matrix]
    nrows = len(values)
    if nrows == 0 or any(len(row) != nrows for row in values):
        raise ValueError("Jacobi eigenvalue input must be a non-empty square matrix")

    tolerance = 1.0e-13
    max_iterations = max(100, 30 * nrows * nrows)
    for _ in range(max_iterations):
        p = 0
        q = 1
        max_offdiag = 0.0
        for irow in range(nrows):
            for jcol in range(irow + 1, nrows):
                offdiag = abs(values[irow][jcol])
                if offdiag > max_offdiag:
                    max_offdiag = offdiag
                    p = irow
                    q = jcol
        if max_offdiag <= tolerance:
            break

        app = values[p][p]
        aqq = values[q][q]
        apq = values[p][q]
        tau = (aqq - app) / (2.0 * apq)
        if tau >= 0.0:
            tangent = 1.0 / (tau + math.sqrt(1.0 + tau * tau))
        else:
            tangent = -1.0 / (-tau + math.sqrt(1.0 + tau * tau))
        cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
        sine = tangent * cosine

        for krow in range(nrows):
            if krow in {p, q}:
                continue
            akp = values[krow][p]
            akq = values[krow][q]
            values[krow][p] = cosine * akp - sine * akq
            values[p][krow] = values[krow][p]
            values[krow][q] = sine * akp + cosine * akq
            values[q][krow] = values[krow][q]
        values[p][p] = (
            cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq
        )
        values[q][q] = (
            sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq
        )
        values[p][q] = 0.0
        values[q][p] = 0.0

    return tuple(values[irow][irow] for irow in range(nrows))


def _singular_values(matrix: Sequence[Sequence[complex]]) -> tuple[float, ...]:
    nrows = len(matrix)
    if nrows == 0:
        return ()
    ncols = len(matrix[0])
    if any(len(row) != ncols for row in matrix):
        raise ValueError("singular value input must be rectangular")

    gram = [[0.0j for _ in range(ncols)] for _ in range(ncols)]
    for icol in range(ncols):
        for jcol in range(ncols):
            value = 0.0j
            for irow in range(nrows):
                value += matrix[irow][icol].conjugate() * matrix[irow][jcol]
            gram[icol][jcol] = value

    doubled = [[0.0 for _ in range(2 * ncols)] for _ in range(2 * ncols)]
    for irow in range(ncols):
        for jcol in range(ncols):
            real_part = gram[irow][jcol].real
            imag_part = gram[irow][jcol].imag
            doubled[irow][jcol] = real_part
            doubled[irow][jcol + ncols] = -imag_part
            doubled[irow + ncols][jcol] = imag_part
            doubled[irow + ncols][jcol + ncols] = real_part

    eigenvalues = sorted(_jacobi_symmetric_eigenvalues(doubled), reverse=True)
    singular_values = []
    for icol in range(ncols):
        eig_pair = 0.5 * (eigenvalues[2 * icol] + eigenvalues[2 * icol + 1])
        singular_values.append(math.sqrt(max(eig_pair, 0.0)))
    return tuple(singular_values)


def _orthonormalized_columns(
    matrix: Sequence[Sequence[complex]], tolerance: float = 1.0e-12
) -> tuple[tuple[complex, ...], ...]:
    nrows = len(matrix)
    if nrows == 0:
        return ()
    ncols = len(matrix[0])
    if any(len(row) != ncols for row in matrix):
        raise ValueError("orthonormalization input must be rectangular")

    columns = [[matrix[irow][icol] for irow in range(nrows)] for icol in range(ncols)]
    q_columns: list[list[complex]] = []
    for column in columns:
        vector = list(column)
        for q_column in q_columns:
            overlap = sum(
                q_column[irow].conjugate() * vector[irow] for irow in range(nrows)
            )
            vector = [
                value - overlap * q_column[irow] for irow, value in enumerate(vector)
            ]
        norm = math.sqrt(sum(abs(value) ** 2 for value in vector))
        if norm <= tolerance:
            q_columns.append([0.0j for _ in range(nrows)])
        else:
            q_columns.append([value / norm for value in vector])
    return tuple(
        tuple(q_columns[icol][irow] for icol in range(ncols)) for irow in range(nrows)
    )


def _projection_subspace_singular_values(
    left_matrix: Sequence[Sequence[complex]], right_matrix: Sequence[Sequence[complex]]
) -> tuple[float, ...]:
    left_q = _orthonormalized_columns(left_matrix)
    right_q = _orthonormalized_columns(right_matrix)
    nrows = len(left_q)
    if nrows == 0 or len(right_q) != nrows:
        raise ValueError("projection subspace inputs must have the same row count")
    nleft = len(left_q[0])
    nright = len(right_q[0])
    overlap = [[0.0j for _ in range(nright)] for _ in range(nleft)]
    for ileft in range(nleft):
        for iright in range(nright):
            overlap[ileft][iright] = sum(
                left_q[irow][ileft].conjugate() * right_q[irow][iright]
                for irow in range(nrows)
            )
    return _singular_values(overlap)


def compare_amn_files(
    left_seed: str, left_workdir: Path, right_seed: str, right_workdir: Path
) -> AmnComparison:
    left_amn = parse_amn_matrix(left_seed, left_workdir)
    right_amn = parse_amn_matrix(right_seed, right_workdir)
    if (left_amn.num_bands, left_amn.num_kpts, left_amn.num_wann) != (
        right_amn.num_bands,
        right_amn.num_kpts,
        right_amn.num_wann,
    ):
        raise ValueError(
            "amn dimensions differ: "
            f"{left_amn.num_bands}x{left_amn.num_kpts}x{left_amn.num_wann} != "
            f"{right_amn.num_bands}x{right_amn.num_kpts}x{right_amn.num_wann}"
        )

    max_abs_delta = 0.0
    sum_sq_delta = 0.0
    entries = 0
    max_singular_delta = 0.0
    sum_sq_singular_delta = 0.0
    singular_entries = 0
    max_subspace_deviation = 0.0
    sum_sq_subspace_deviation = 0.0
    subspace_entries = 0
    min_subspace_svalue = float("inf")
    max_subspace_kpt_index = 0
    max_kpt_delta = 0.0
    max_kpt_index = 0
    for ikpt, (left_matrix, right_matrix) in enumerate(
        zip(left_amn.matrices, right_amn.matrices), start=1
    ):
        kpt_delta = 0.0
        for iband in range(left_amn.num_bands):
            for iwann in range(left_amn.num_wann):
                delta = abs(left_matrix[iband][iwann] - right_matrix[iband][iwann])
                max_abs_delta = max(max_abs_delta, delta)
                kpt_delta = max(kpt_delta, delta)
                sum_sq_delta += delta * delta
                entries += 1
        if kpt_delta > max_kpt_delta:
            max_kpt_delta = kpt_delta
            max_kpt_index = ikpt

        left_singular = _singular_values(left_matrix)
        right_singular = _singular_values(right_matrix)
        for left_value, right_value in zip(left_singular, right_singular):
            singular_delta = abs(left_value - right_value)
            max_singular_delta = max(max_singular_delta, singular_delta)
            sum_sq_singular_delta += singular_delta * singular_delta
            singular_entries += 1
        subspace_singular = _projection_subspace_singular_values(
            left_matrix, right_matrix
        )
        kpt_subspace_deviation = 0.0
        for singular_value in subspace_singular:
            deviation = abs(1.0 - singular_value)
            kpt_subspace_deviation = max(kpt_subspace_deviation, deviation)
            min_subspace_svalue = min(min_subspace_svalue, singular_value)
            sum_sq_subspace_deviation += deviation * deviation
            subspace_entries += 1
        if kpt_subspace_deviation > max_subspace_deviation:
            max_subspace_deviation = kpt_subspace_deviation
            max_subspace_kpt_index = ikpt

    return AmnComparison(
        max_abs_delta=max_abs_delta,
        rms_abs_delta=math.sqrt(sum_sq_delta / float(entries)),
        max_singular_delta=max_singular_delta,
        rms_singular_delta=math.sqrt(sum_sq_singular_delta / float(singular_entries)),
        max_subspace_deviation=max_subspace_deviation,
        rms_subspace_deviation=math.sqrt(
            sum_sq_subspace_deviation / float(subspace_entries)
        ),
        min_subspace_svalue=min_subspace_svalue,
        max_subspace_kpt_index=max_subspace_kpt_index,
        max_kpt_delta=max_kpt_delta,
        max_kpt_index=max_kpt_index,
    )


def compare_exports(
    left_seed: str, left_workdir: Path, right_seed: str, right_workdir: Path
) -> ExportComparison:
    left_eig = parse_eig(left_seed, left_workdir)
    right_eig = parse_eig(right_seed, right_workdir)
    if (left_eig.num_bands, left_eig.num_kpts) != (
        right_eig.num_bands,
        right_eig.num_kpts,
    ):
        raise ValueError(
            "eig dimensions differ: "
            f"{left_eig.num_bands}x{left_eig.num_kpts} != "
            f"{right_eig.num_bands}x{right_eig.num_kpts}"
        )
    eig_max_abs_delta = 0.0
    for ikpt in range(left_eig.num_kpts):
        for iband in range(left_eig.num_bands):
            eig_max_abs_delta = max(
                eig_max_abs_delta,
                abs(left_eig.values[ikpt][iband] - right_eig.values[ikpt][iband]),
            )

    left_mmn = parse_mmn(left_seed, left_workdir)
    right_mmn = parse_mmn(right_seed, right_workdir)
    if (
        left_mmn.num_bands,
        left_mmn.num_kpts,
        left_mmn.num_neighbors,
        len(left_mmn.blocks),
    ) != (
        right_mmn.num_bands,
        right_mmn.num_kpts,
        right_mmn.num_neighbors,
        len(right_mmn.blocks),
    ):
        raise ValueError(
            "mmn dimensions differ: "
            f"{left_mmn.num_bands}x{left_mmn.num_kpts}x{left_mmn.num_neighbors} != "
            f"{right_mmn.num_bands}x{right_mmn.num_kpts}x{right_mmn.num_neighbors}"
        )

    mmn_max_abs_delta = 0.0
    mmn_sum_sq_delta = 0.0
    mmn_entries = 0
    mmn_max_singular_delta = 0.0
    mmn_sum_sq_singular_delta = 0.0
    mmn_singular_entries = 0
    mmn_max_block_delta = 0.0
    mmn_max_block_index = 0
    for iblock, (left_block, right_block) in enumerate(
        zip(left_mmn.blocks, right_mmn.blocks), start=1
    ):
        if (left_block.ikpt, left_block.jkpt, left_block.cell) != (
            right_block.ikpt,
            right_block.jkpt,
            right_block.cell,
        ):
            raise ValueError(
                "mmn block metadata differ at block "
                f"{iblock}: {left_block.ikpt}/{left_block.jkpt}/{left_block.cell} "
                f"!= {right_block.ikpt}/{right_block.jkpt}/{right_block.cell}"
            )
        block_delta = 0.0
        for iband in range(left_mmn.num_bands):
            for jband in range(left_mmn.num_bands):
                delta = abs(
                    left_block.values[iband][jband] - right_block.values[iband][jband]
                )
                mmn_max_abs_delta = max(mmn_max_abs_delta, delta)
                block_delta = max(block_delta, delta)
                mmn_sum_sq_delta += delta * delta
                mmn_entries += 1
        if block_delta > mmn_max_block_delta:
            mmn_max_block_delta = block_delta
            mmn_max_block_index = iblock
        left_singular = _singular_values(left_block.values)
        right_singular = _singular_values(right_block.values)
        for left_value, right_value in zip(left_singular, right_singular):
            singular_delta = abs(left_value - right_value)
            mmn_max_singular_delta = max(mmn_max_singular_delta, singular_delta)
            mmn_sum_sq_singular_delta += singular_delta * singular_delta
            mmn_singular_entries += 1

    return ExportComparison(
        eig_max_abs_delta=eig_max_abs_delta,
        mmn_max_abs_delta=mmn_max_abs_delta,
        mmn_rms_abs_delta=math.sqrt(mmn_sum_sq_delta / float(mmn_entries)),
        mmn_max_singular_delta=mmn_max_singular_delta,
        mmn_rms_singular_delta=math.sqrt(
            mmn_sum_sq_singular_delta / float(mmn_singular_entries)
        ),
        mmn_max_block_delta=mmn_max_block_delta,
        mmn_max_block_index=mmn_max_block_index,
    )


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
    periodic_cell = None
    if args.periodic_centers:
        periodic_cell = parse_win_cell(args.left_seed, Path(args.left_workdir))
    comparison = compare_summaries(
        left, right, match_centers=bool(args.match_centers), periodic_cell=periodic_cell
    )
    print(f"omega_total_delta = {comparison.omega_total_delta:.12e}")
    print(f"spread_delta = {comparison.spread_delta:.12e}")
    print(f"center_delta = {comparison.center_delta:.12e}")
    print("center_matching = " f"{'optimal' if args.match_centers else 'direct'}")
    print(f"periodic_centers = {'yes' if args.periodic_centers else 'no'}")
    if args.match_centers:
        pairs = ", ".join(
            f"{ileft + 1}->{iright + 1}"
            for ileft, iright in enumerate(comparison.assignment)
        )
        print(f"center_assignment = {pairs}")

    failed = False
    if (
        args.omega_total_tol is not None
        and comparison.omega_total_delta > args.omega_total_tol
    ):
        failed = True
    if args.spread_tol is not None and comparison.spread_delta > args.spread_tol:
        failed = True
    if args.center_tol is not None and comparison.center_delta > args.center_tol:
        failed = True
    return 1 if failed else 0


def command_compare_export(args: argparse.Namespace) -> int:
    comparison = compare_exports(
        args.left_seed,
        Path(args.left_workdir),
        args.right_seed,
        Path(args.right_workdir),
    )
    print(f"eig_max_abs_delta = {comparison.eig_max_abs_delta:.12e}")
    print(f"mmn_max_abs_delta = {comparison.mmn_max_abs_delta:.12e}")
    print(f"mmn_rms_abs_delta = {comparison.mmn_rms_abs_delta:.12e}")
    print(f"mmn_max_singular_delta = {comparison.mmn_max_singular_delta:.12e}")
    print(f"mmn_rms_singular_delta = {comparison.mmn_rms_singular_delta:.12e}")
    print(f"mmn_max_block_delta = {comparison.mmn_max_block_delta:.12e}")
    print(f"mmn_max_block_index = {comparison.mmn_max_block_index}")
    failed = False
    if args.eig_tol is not None and comparison.eig_max_abs_delta > args.eig_tol:
        failed = True
    if args.mmn_tol is not None and comparison.mmn_max_abs_delta > args.mmn_tol:
        failed = True
    if (
        args.mmn_singular_tol is not None
        and comparison.mmn_max_singular_delta > args.mmn_singular_tol
    ):
        failed = True
    return 1 if failed else 0


def command_compare_amn(args: argparse.Namespace) -> int:
    comparison = compare_amn_files(
        args.left_seed,
        Path(args.left_workdir),
        args.right_seed,
        Path(args.right_workdir),
    )
    print(f"amn_max_abs_delta = {comparison.max_abs_delta:.12e}")
    print(f"amn_rms_abs_delta = {comparison.rms_abs_delta:.12e}")
    print(f"amn_max_singular_delta = {comparison.max_singular_delta:.12e}")
    print(f"amn_rms_singular_delta = {comparison.rms_singular_delta:.12e}")
    print(f"amn_max_subspace_deviation = {comparison.max_subspace_deviation:.12e}")
    print(f"amn_rms_subspace_deviation = {comparison.rms_subspace_deviation:.12e}")
    print(f"amn_min_subspace_svalue = {comparison.min_subspace_svalue:.12e}")
    print(f"amn_max_subspace_kpt_index = {comparison.max_subspace_kpt_index}")
    print(f"amn_max_kpt_delta = {comparison.max_kpt_delta:.12e}")
    print(f"amn_max_kpt_index = {comparison.max_kpt_index}")
    failed = False
    if args.amn_tol is not None and comparison.max_abs_delta > args.amn_tol:
        failed = True
    if (
        args.amn_singular_tol is not None
        and comparison.max_singular_delta > args.amn_singular_tol
    ):
        failed = True
    if (
        args.amn_subspace_tol is not None
        and comparison.max_subspace_deviation > args.amn_subspace_tol
    ):
        failed = True
    return 1 if failed else 0


def command_inspect_amn(args: argparse.Namespace) -> int:
    summary = parse_amn(args.seed, Path(args.workdir), args.identity_tol)
    print(f"seed = {summary.seed}")
    print(f"amn = {summary.path}")
    print(f"num_bands = {summary.num_bands}")
    print(f"num_kpts = {summary.num_kpts}")
    print(f"num_wann = {summary.num_wann}")
    print(f"max_identity_deviation = {summary.max_identity_deviation:.12e}")
    print(f"max_column_norm_delta = {summary.max_column_norm_delta:.12e}")
    print(f"max_offdiag = {summary.max_offdiag:.12e}")
    print(f"nonzero_fraction = {summary.nonzero_fraction:.12e}")
    print(f"is_identity = {'yes' if summary.is_identity else 'no'}")
    if args.require_identity and not summary.is_identity:
        return 1
    if args.max_identity_deviation is not None and (
        summary.max_identity_deviation > args.max_identity_deviation
    ):
        return 1
    return 0


def command_prepare_variant(args: argparse.Namespace) -> int:
    result = prepare_variant(
        args.seed,
        Path(args.source_workdir),
        Path(args.output_workdir),
        args.variant,
        force=bool(args.force),
    )
    print(f"seed = {result.seed}")
    print(f"variant = {result.variant}")
    print(f"workdir = {result.workdir}")
    print(f"removed_amn = {'yes' if result.removed_amn else 'no'}")
    return 0


def _zero_iter_label(label: str) -> str | None:
    if label.endswith("-zero-iter"):
        return None
    if ITERATION_LABEL_RE.match(label):
        return None
    if label == "source":
        return "zero-iter"
    return f"{label}-zero-iter"


def print_relaxation_summary(
    specs: Sequence[RunSpec],
    reference_label: str,
    *,
    allow_failed: bool,
    match_centers: bool,
    periodic_cell: W90Cell | None,
) -> None:
    spec_by_label = {spec.label: spec for spec in specs}
    if reference_label not in spec_by_label:
        raise ValueError(f"unknown relaxation reference label: {reference_label}")
    reference_spec = spec_by_label[reference_label]
    reference = parse_wout(reference_spec.seed, reference_spec.workdir)

    header = [
        "final_label",
        "start_label",
        "status",
        "start_omega_total",
        "final_omega_total",
        "omega_relax_delta",
        "start_vs_reference_delta",
        "final_vs_reference_delta",
        "spread_relax_delta",
        "center_relax_delta",
        "error",
    ]
    print("")
    print("# relaxation_summary")
    print("\t".join(header))

    for final_spec in specs:
        start_label = _zero_iter_label(final_spec.label)
        if start_label is None or start_label not in spec_by_label:
            continue
        start_spec = spec_by_label[start_label]
        try:
            start = parse_wout(start_spec.seed, start_spec.workdir)
            final = parse_wout(final_spec.seed, final_spec.workdir)
        except (OSError, ValueError) as exc:
            if not allow_failed:
                raise
            row = [
                final_spec.label,
                start_label,
                "failed",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                str(exc),
            ]
            print("\t".join(row))
            continue

        relaxation = compare_summaries(
            start, final, match_centers=match_centers, periodic_cell=periodic_cell
        )
        start_reference = compare_summaries(
            reference, start, match_centers=match_centers, periodic_cell=periodic_cell
        )
        final_reference = compare_summaries(
            reference, final, match_centers=match_centers, periodic_cell=periodic_cell
        )
        row = [
            final_spec.label,
            start_label,
            "ok",
            f"{start.omega_total:.12e}",
            f"{final.omega_total:.12e}",
            f"{final.omega_total - start.omega_total:.12e}",
            f"{start_reference.omega_total_delta:.12e}",
            f"{final_reference.omega_total_delta:.12e}",
            f"{relaxation.spread_delta:.12e}",
            f"{relaxation.center_delta:.12e}",
            "",
        ]
        print("\t".join(row))


def print_iteration_summary(
    specs: Sequence[RunSpec],
    reference_label: str,
    *,
    allow_failed: bool,
    match_centers: bool,
    periodic_cell: W90Cell | None,
) -> None:
    spec_by_label = {spec.label: spec for spec in specs}
    if reference_label not in spec_by_label:
        raise ValueError(f"unknown iteration reference label: {reference_label}")
    reference_spec = spec_by_label[reference_label]
    reference = parse_wout(reference_spec.seed, reference_spec.workdir)

    iteration_specs: list[tuple[str, int, RunSpec]] = []
    for spec in specs:
        match = ITERATION_LABEL_RE.match(spec.label)
        if match is None:
            continue
        iteration_specs.append(
            (match.group("variant"), int(match.group("iteration")), spec)
        )
    if not iteration_specs:
        return

    header = [
        "variant",
        "iteration",
        "status",
        "omega_total",
        "omega_from_start_delta",
        "omega_from_final_delta",
        "omega_vs_reference_delta",
        "spread_vs_reference_delta",
        "center_vs_reference_delta",
        "error",
    ]
    print("")
    print("# iteration_summary")
    print("\t".join(header))

    for variant, iteration, spec in sorted(iteration_specs):
        try:
            summary = parse_wout(spec.seed, spec.workdir)
            start = parse_wout(
                spec_by_label[f"{variant}-zero-iter"].seed,
                spec_by_label[f"{variant}-zero-iter"].workdir,
            )
            final = parse_wout(
                spec_by_label[variant].seed, spec_by_label[variant].workdir
            )
        except (KeyError, OSError, ValueError) as exc:
            if not allow_failed:
                raise
            row = [
                variant,
                str(iteration),
                "failed",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                str(exc),
            ]
            print("\t".join(row))
            continue

        reference_comparison = compare_summaries(
            reference, summary, match_centers=match_centers, periodic_cell=periodic_cell
        )
        row = [
            variant,
            str(iteration),
            "ok",
            f"{summary.omega_total:.12e}",
            f"{summary.omega_total - start.omega_total:.12e}",
            f"{summary.omega_total - final.omega_total:.12e}",
            f"{reference_comparison.omega_total_delta:.12e}",
            f"{reference_comparison.spread_delta:.12e}",
            f"{reference_comparison.center_delta:.12e}",
            "",
        ]
        print("\t".join(row))


def _is_strategy_candidate(label: str) -> bool:
    return (
        not label.endswith("-zero-iter")
        and "-cap-" not in label
        and ITERATION_LABEL_RE.match(label) is None
    )


def _strategy_action(history: WoutHistorySummary) -> str:
    if history.classification == "late_drift":
        return f"cap_num_iter={history.best.iteration}"
    if history.classification == "zero_iter":
        return "start_only"
    return "full_minimization"


def _has_wannier90_export(seed: str, workdir: Path) -> bool:
    return all(
        (workdir / f"{seed}.{suffix}").exists() for suffix in ("win", "eig", "mmn")
    )


def _safe_path_component(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "case"


def selected_scan_variants(
    variants: Sequence[str] | None, *, extended_strategy_scan: bool
) -> tuple[str, ...]:
    if variants:
        return tuple(variants)
    if extended_strategy_scan:
        return EXTENDED_SCAN_VARIANTS
    return DEFAULT_SCAN_VARIANTS


def read_strategy_cases(
    case_specs: Sequence[str], case_files: Sequence[str] | None
) -> list[StrategyCase]:
    specs = list(case_specs)
    for case_file in case_files or ():
        path = Path(case_file)
        for iline, line in enumerate(
            path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1
        ):
            stripped = line.split("#", maxsplit=1)[0].strip()
            if not stripped:
                continue
            try:
                parse_strategy_case(stripped)
            except ValueError as exc:
                raise ValueError(f"{path}:{iline}: {exc}") from exc
            specs.append(stripped)
    if not specs:
        raise ValueError("provide at least one case or --case-file")
    return [parse_strategy_case(spec) for spec in specs]


def read_cp2k_strategy_cases(
    case_specs: Sequence[str], case_files: Sequence[str] | None
) -> list[CP2KStrategyCase]:
    specs = list(case_specs)
    for case_file in case_files or ():
        path = Path(case_file)
        for iline, line in enumerate(
            path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1
        ):
            stripped = line.split("#", maxsplit=1)[0].strip()
            if not stripped:
                continue
            try:
                parse_cp2k_strategy_case(stripped)
            except ValueError as exc:
                raise ValueError(f"{path}:{iline}: {exc}") from exc
            specs.append(stripped)
    if not specs:
        raise ValueError("provide at least one CP2K case or --case-file")
    return [parse_cp2k_strategy_case(spec) for spec in specs]


def _cp2k_input_value(input_path: Path, keyword: str) -> str | None:
    upper_keyword = keyword.upper()
    for line in input_path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.split("!", maxsplit=1)[0].split("#", maxsplit=1)[0].strip()
        if not stripped:
            continue
        fields = stripped.split()
        if len(fields) >= 2 and fields[0].upper() == upper_keyword:
            return fields[1]
    return None


def infer_cp2k_wannier_seed(case: CP2KStrategyCase) -> str:
    if case.seed is not None:
        return case.seed
    seed = _cp2k_input_value(case.input_path, "SEED_NAME")
    if seed is not None:
        return seed
    project = _cp2k_input_value(case.input_path, "PROJECT")
    if project is not None:
        return project
    return case.input_path.stem


def _cp2k_projection_seed(
    base_seed: str,
    projection: str,
    projection_smoothing: str | None = None,
    auto_projection_score: str | None = None,
) -> str:
    projection = projection.upper()
    suffix_parts = [CP2K_PROJECTION_SEED_SUFFIXES[projection]]
    if auto_projection_score is not None:
        auto_projection_score = auto_projection_score.upper()
        suffix_parts.append(
            CP2K_AUTO_PROJECTION_SCORE_SEED_SUFFIXES[auto_projection_score]
        )
    if projection_smoothing is not None:
        projection_smoothing = projection_smoothing.upper()
        suffix_parts.append(
            CP2K_PROJECTION_SMOOTHING_SEED_SUFFIXES[projection_smoothing]
        )
    suffix = "_".join(suffix_parts)
    seed = f"{base_seed}_{suffix}"
    if len(seed) <= MAX_WANNIER90_SEED_LENGTH:
        return seed
    digest = hashlib.sha1(seed.encode("utf-8")).hexdigest()[:8]
    max_base = MAX_WANNIER90_SEED_LENGTH - len(suffix) - len(digest) - 2
    return f"{base_seed[:max_base]}_{digest}_{suffix}"


def _is_section_start(stripped: str, name: str) -> bool:
    fields = stripped.split()
    return bool(fields) and fields[0].upper() == f"&{name.upper()}"


def _is_section_end(stripped: str, name: str) -> bool:
    fields = stripped.split()
    if not fields:
        return False
    head = fields[0].upper()
    return head == "&END" and (len(fields) == 1 or fields[1].upper() == name.upper())


def _cp2k_input_with_wannier_projection(
    source_path: Path,
    *,
    seed_name: str,
    projection: str,
    projection_smoothing: str | None = None,
    auto_projection_score: str | None = None,
) -> str:
    projection = projection.upper()
    if projection not in CP2K_INITIAL_PROJECTIONS:
        raise ValueError(
            f"unknown CP2K WANNIER90 INITIAL_PROJECTIONS mode: {projection}"
        )
    if projection_smoothing is not None:
        projection_smoothing = projection_smoothing.upper()
        if projection_smoothing not in CP2K_PROJECTION_SMOOTHINGS:
            raise ValueError(
                f"unknown CP2K WANNIER90 PROJECTION_SMOOTHING mode: {projection_smoothing}"
            )
    if auto_projection_score is not None:
        auto_projection_score = auto_projection_score.upper()
        if auto_projection_score not in CP2K_AUTO_PROJECTION_SCORES:
            raise ValueError(
                f"unknown CP2K WANNIER90 AUTO_PROJECTION_SCORE mode: {auto_projection_score}"
            )
        if projection != "AUTO_SCDM_RAW":
            raise ValueError(
                "CP2K WANNIER90 AUTO_PROJECTION_SCORE variants require "
                "INITIAL_PROJECTIONS AUTO_SCDM_RAW"
            )

    output_lines: list[str] = []
    in_wannier90 = False
    found_wannier90 = False
    handled_wannier90 = False
    override_keys = {
        "SEED_NAME",
        "USE_BLOCH_PHASES",
        "INITIAL_PROJECTIONS",
        "AUTO_PROJECTION_SCORE",
    }
    if projection_smoothing is not None:
        override_keys.update({"PROJECTION_SMOOTHING", "SMOOTH_PROJECTIONS"})
    for raw_line in source_path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        stripped = raw_line.split("!", maxsplit=1)[0].split("#", maxsplit=1)[0].strip()
        if _is_section_start(stripped, "WANNIER90"):
            in_wannier90 = True
            found_wannier90 = True
            output_lines.append(raw_line)
            continue
        if in_wannier90 and _is_section_end(stripped, "WANNIER90"):
            indent = re.match(r"\s*", raw_line).group(0)
            keyword_indent = f"{indent}  "
            output_lines.append(f"{keyword_indent}SEED_NAME {seed_name}")
            if projection == "IDENTITY":
                output_lines.append(f"{keyword_indent}USE_BLOCH_PHASES T")
            output_lines.append(f"{keyword_indent}INITIAL_PROJECTIONS {projection}")
            if auto_projection_score is not None:
                output_lines.append(
                    f"{keyword_indent}AUTO_PROJECTION_SCORE {auto_projection_score}"
                )
            if projection_smoothing is not None:
                output_lines.append(
                    f"{keyword_indent}PROJECTION_SMOOTHING {projection_smoothing}"
                )
            output_lines.append(raw_line)
            in_wannier90 = False
            handled_wannier90 = True
            continue
        if in_wannier90:
            fields = stripped.split()
            if fields and fields[0].upper() in override_keys:
                continue
        output_lines.append(raw_line)

    if not found_wannier90:
        raise ValueError(f"{source_path}: missing DFT%PRINT%WANNIER90 section")
    if not handled_wannier90:
        raise ValueError(f"{source_path}: unterminated WANNIER90 section")
    return "\n".join(output_lines) + "\n"


def expand_cp2k_projection_cases(
    cases: Sequence[CP2KStrategyCase],
    input_root: Path,
    projections: Sequence[str] | None,
    projection_smoothings: Sequence[str] | None,
    auto_projection_scores: Sequence[str] | None,
    *,
    write_inputs: bool,
) -> list[CP2KStrategyCase]:
    if auto_projection_scores and not projections:
        raise ValueError(
            "--cp2k-auto-projection-score requires --cp2k-initial-projection AUTO_SCDM_RAW"
        )
    if auto_projection_scores:
        has_auto_projection = any(
            projection.upper() == "AUTO_SCDM_RAW" for projection in projections or ()
        )
        if not has_auto_projection:
            raise ValueError(
                "--cp2k-auto-projection-score requires "
                "--cp2k-initial-projection AUTO_SCDM_RAW"
            )
    if not projections and not projection_smoothings and not auto_projection_scores:
        return list(cases)
    projection_list = tuple(projections) if projections else None
    smoothings = tuple(projection_smoothings or (None,))
    auto_scores = tuple(auto_projection_scores or (None,))
    expanded: list[CP2KStrategyCase] = []
    if write_inputs:
        input_root.mkdir(parents=True, exist_ok=True)
    for case in cases:
        base_seed = infer_cp2k_wannier_seed(case)
        case_projections = projection_list or (
            _cp2k_input_value(case.input_path, "INITIAL_PROJECTIONS") or "IDENTITY",
        )
        for projection in case_projections:
            projection = projection.upper()
            case_auto_scores = auto_scores if projection == "AUTO_SCDM_RAW" else (None,)
            for auto_score in case_auto_scores:
                auto_score_label = ""
                auto_score_mode = None
                if auto_score is not None:
                    auto_score_mode = auto_score.upper()
                    if auto_score_mode not in CP2K_AUTO_PROJECTION_SCORES:
                        raise ValueError(
                            f"unknown CP2K WANNIER90 AUTO_PROJECTION_SCORE mode: {auto_score_mode}"
                        )
                    auto_score_label = f"_{auto_score_mode.lower()}"
                for smoothing in smoothings:
                    smoothing_label = ""
                    smoothing_mode = None
                    if smoothing is not None:
                        smoothing_mode = smoothing.upper()
                        if smoothing_mode not in CP2K_PROJECTION_SMOOTHINGS:
                            raise ValueError(
                                f"unknown CP2K WANNIER90 PROJECTION_SMOOTHING mode: {smoothing_mode}"
                            )
                        smoothing_label = f"_{smoothing_mode.lower()}"
                    seed = _cp2k_projection_seed(
                        base_seed, projection, smoothing_mode, auto_score_mode
                    )
                    label = f"{case.label}_{projection.lower()}{auto_score_label}{smoothing_label}"
                    input_path = input_root / f"{_safe_path_component(label)}.inp"
                    if write_inputs:
                        input_path.write_text(
                            _cp2k_input_with_wannier_projection(
                                case.input_path,
                                seed_name=seed,
                                projection=projection,
                                projection_smoothing=smoothing_mode,
                                auto_projection_score=auto_score_mode,
                            ),
                            encoding="utf-8",
                        )
                    expanded.append(
                        CP2KStrategyCase(
                            label=label,
                            input_path=input_path if write_inputs else case.input_path,
                            seed=seed,
                        )
                    )
    return expanded


def run_cp2k_export(
    executable: str,
    case: CP2KStrategyCase,
    output_workdir: Path,
    *,
    timeout: float,
    copy_input_dir_files: bool,
) -> None:
    if not case.input_path.exists():
        raise ValueError(f"{case.input_path}: CP2K input file does not exist")
    output_workdir.mkdir(parents=True, exist_ok=True)
    if copy_input_dir_files:
        for source in case.input_path.parent.iterdir():
            if source.is_file():
                shutil.copy2(source, output_workdir / source.name)
    else:
        shutil.copy2(case.input_path, output_workdir / case.input_path.name)

    output_name = f"{_safe_path_component(case.label)}.cp2k.out"
    result = subprocess.run(
        [executable, "-i", case.input_path.name, "-o", output_name],
        cwd=output_workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        check=False,
    )
    (output_workdir / f"{_safe_path_component(case.label)}.cp2k.stdout").write_text(
        result.stdout, encoding="utf-8"
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"{executable} -i {case.input_path.name} failed with exit code "
            f"{result.returncode}"
        )


def _strategy_matrix_row(values: Sequence[str]) -> dict[str, str]:
    if len(values) != len(STRATEGY_MATRIX_HEADER):
        raise ValueError("internal strategy-matrix row length mismatch")
    return dict(zip(STRATEGY_MATRIX_HEADER, values))


def _print_strategy_matrix_row(row: dict[str, str]) -> None:
    print("\t".join(row[column] for column in STRATEGY_MATRIX_HEADER))


def _float_from_row(row: dict[str, str], field: str) -> float:
    try:
        return float(row[field])
    except (KeyError, ValueError):
        return float("nan")


def _strategy_cap_verified(row: dict[str, str], cap_tol: float) -> bool:
    if row["cap_label"] == "na" or row["cap_classification"] == "late_drift":
        return False
    for field in STRATEGY_MATRIX_CAP_TOL_FIELDS:
        value = _float_from_row(row, field)
        if not math.isfinite(value) or abs(value) > cap_tol:
            return False
    return True


def strategy_matrix_assessment(row: dict[str, str], *, cap_tol: float) -> str:
    if row["status"] == "export_only":
        return "export_only"
    if row["status"] != "ok":
        return "failed"
    classification = row["selected_classification"]
    if classification == "stable":
        return "stable"
    if classification == "zero_iter":
        return "zero_iter"
    if classification == "late_drift":
        if _strategy_cap_verified(row, cap_tol):
            return "cap_verified"
        return "late_drift_unverified"
    return "unknown"


def strategy_matrix_recommendation(row: dict[str, str], *, cap_tol: float) -> str:
    assessment = strategy_matrix_assessment(row, cap_tol=cap_tol)
    if assessment == "stable":
        return "use selected full minimization"
    if assessment == "zero_iter":
        return "use start-only export"
    if assessment == "cap_verified":
        return f"use capped minimization from {row['cap_label']}"
    if assessment == "late_drift_unverified":
        return "run/inspect capped minimization before accepting this strategy"
    if assessment == "export_only":
        return "CP2K export only; no external MLWF strategy candidates"
    if assessment == "failed":
        return "inspect failed case"
    return "inspect selected strategy"


def strategy_matrix_assessments(
    rows: Sequence[dict[str, str]], *, cap_tol: float
) -> list[dict[str, str]]:
    assessments: list[dict[str, str]] = []
    for row in rows:
        assessments.append(
            {
                "case": row["case"],
                "seed": row["seed"],
                "assessment": strategy_matrix_assessment(row, cap_tol=cap_tol),
                "recommendation": strategy_matrix_recommendation(row, cap_tol=cap_tol),
            }
        )
    return assessments


def _split_cp2k_projection_case(case_label: str) -> tuple[str, str] | None:
    for projection in sorted(CP2K_INITIAL_PROJECTIONS, key=len, reverse=True):
        suffix = f"_{projection.lower()}"
        if case_label.endswith(suffix):
            return case_label[: -len(suffix)], projection
    return None


def _best_projection_row(rows: Sequence[dict[str, str]]) -> dict[str, str] | None:
    rankable = [
        row
        for row in rows
        if row["status"] == "ok"
        and math.isfinite(_float_from_row(row, "selected_best_omega_total"))
    ]
    if not rankable:
        return None
    return min(
        rankable, key=lambda row: _float_from_row(row, "selected_best_omega_total")
    )


def summarize_strategy_projection_families(
    rows: Sequence[dict[str, str]], *, cap_tol: float
) -> list[dict[str, str]]:
    grouped: dict[str, list[tuple[str, dict[str, str]]]] = {}
    for row in rows:
        split = _split_cp2k_projection_case(row["case"])
        if split is None:
            continue
        base_case, projection = split
        grouped.setdefault(base_case, []).append((projection, row))

    summary_rows: list[dict[str, str]] = []
    for base_case in sorted(grouped):
        entries = grouped[base_case]
        if len(entries) < 2:
            continue
        projection_rows = [row for _, row in entries]
        best = _best_projection_row(projection_rows)
        ortho_best = _best_projection_row(
            [row for projection, row in entries if projection.endswith("_ORTHO")]
        )
        nonortho_best = _best_projection_row(
            [row for projection, row in entries if not projection.endswith("_ORTHO")]
        )
        raw_best = _best_projection_row(
            [row for projection, row in entries if projection.endswith("_RAW")]
        )
        nonraw_best = _best_projection_row(
            [row for projection, row in entries if not projection.endswith("_RAW")]
        )
        if best is None:
            continue
        best_split = _split_cp2k_projection_case(best["case"])
        best_projection = best_split[1] if best_split is not None else "na"
        best_omega = _float_from_row(best, "selected_best_omega_total")
        if ortho_best is not None and nonortho_best is not None:
            ortho_delta = _float_from_row(
                ortho_best, "selected_best_omega_total"
            ) - _float_from_row(nonortho_best, "selected_best_omega_total")
            ortho_delta_text = f"{ortho_delta:.12e}"
        else:
            ortho_delta_text = "nan"
        if raw_best is not None and nonraw_best is not None:
            raw_delta = _float_from_row(
                raw_best, "selected_best_omega_total"
            ) - _float_from_row(nonraw_best, "selected_best_omega_total")
            raw_delta_text = f"{raw_delta:.12e}"
        else:
            raw_delta_text = "nan"
        summary_rows.append(
            {
                "base_case": base_case,
                "num_projection_cases": str(len(entries)),
                "num_ok": str(sum(row["status"] == "ok" for row in projection_rows)),
                "num_export_only": str(
                    sum(row["status"] == "export_only" for row in projection_rows)
                ),
                "num_failed": str(
                    sum(
                        row["status"] not in {"ok", "export_only"}
                        for row in projection_rows
                    )
                ),
                "best_case": best["case"],
                "best_projection": best_projection,
                "best_omega_total": f"{best_omega:.12e}",
                "best_assessment": strategy_matrix_assessment(best, cap_tol=cap_tol),
                "best_recommendation": strategy_matrix_recommendation(
                    best, cap_tol=cap_tol
                ),
                "best_ortho_minus_best_nonortho": ortho_delta_text,
                "best_raw_minus_best_nonraw": raw_delta_text,
            }
        )
    return summary_rows


def summarize_strategy_matrix(
    rows: Sequence[dict[str, str]], *, cap_tol: float
) -> dict[str, Any]:
    ok_rows = [row for row in rows if row["status"] == "ok"]
    export_only_rows = [row for row in rows if row["status"] == "export_only"]
    failed_rows = [row for row in rows if row["status"] not in {"ok", "export_only"}]
    late_drift_rows = [
        row for row in ok_rows if row["selected_classification"] == "late_drift"
    ]
    capped_rows = [row for row in ok_rows if row["cap_label"] != "na"]
    verified_cap_rows = [
        row for row in capped_rows if _strategy_cap_verified(row, cap_tol)
    ]
    selected_other_rows = [
        row
        for row in ok_rows
        if row["selected_classification"] not in {"stable", "late_drift", "zero_iter"}
    ]
    assessments = [strategy_matrix_assessment(row, cap_tol=cap_tol) for row in rows]
    return {
        "num_cases": len(rows),
        "num_ok": len(ok_rows),
        "num_export_only": len(export_only_rows),
        "num_failed": len(failed_rows),
        "num_selected_stable": sum(
            row["selected_classification"] == "stable" for row in ok_rows
        ),
        "num_selected_late_drift": len(late_drift_rows),
        "num_selected_zero_iter": sum(
            row["selected_classification"] == "zero_iter" for row in ok_rows
        ),
        "num_selected_other": len(selected_other_rows),
        "num_failed_candidates": sum(
            int(row["num_failed_candidates"])
            for row in rows
            if row["num_failed_candidates"].isdigit()
        ),
        "num_capped": len(capped_rows),
        "num_verified_caps": len(verified_cap_rows),
        "num_unverified_caps": len(capped_rows) - len(verified_cap_rows),
        "num_uncapped_late_drift": sum(
            not _strategy_cap_verified(row, cap_tol) for row in late_drift_rows
        ),
        "num_assessed_stable": assessments.count("stable"),
        "num_assessed_cap_verified": assessments.count("cap_verified"),
        "num_assessed_late_drift_unverified": assessments.count(
            "late_drift_unverified"
        ),
        "num_assessed_zero_iter": assessments.count("zero_iter"),
        "num_assessed_export_only": assessments.count("export_only"),
        "num_assessed_failed": assessments.count("failed"),
        "num_assessed_unknown": assessments.count("unknown"),
    }


def print_strategy_matrix_summary(summary: dict[str, Any]) -> None:
    print("")
    print("# strategy_matrix_summary")
    print("metric\tvalue")
    for key in sorted(summary):
        print(f"{key}\t{summary[key]}")


def print_strategy_projection_family_summary(
    projection_summary: Sequence[dict[str, str]],
) -> None:
    if not projection_summary:
        return
    columns = (
        "base_case",
        "num_projection_cases",
        "num_ok",
        "num_export_only",
        "num_failed",
        "best_projection",
        "best_omega_total",
        "best_assessment",
        "best_ortho_minus_best_nonortho",
        "best_raw_minus_best_nonraw",
    )
    print("")
    print("# strategy_projection_family_summary")
    print("\t".join(columns))
    for row in projection_summary:
        print("\t".join(row[column] for column in columns))


def _markdown_cell(value: Any) -> str:
    text = str(value)
    return text.replace("|", "\\|").replace("\n", " ")


def strategy_matrix_markdown(
    rows: Sequence[dict[str, str]],
    summary: dict[str, Any],
    projection_summary: Sequence[dict[str, str]],
    *,
    cap_tol: float,
) -> str:
    lines = [
        "# Wannier90 Strategy Matrix Report",
        "",
        "## Summary",
        "",
        "| Metric | Value |",
        "| --- | ---: |",
    ]
    for key in sorted(summary):
        lines.append(f"| `{_markdown_cell(key)}` | {_markdown_cell(summary[key])} |")

    if projection_summary:
        lines.extend(
            [
                "",
                "## Projection Families",
                "",
                (
                    "| Base case | Cases | OK | Export-only | Failed | Best projection | "
                    "Best omega | Assessment | Ortho - non-ortho | Raw - non-raw |"
                ),
                "| --- | ---: | ---: | ---: | ---: | --- | ---: | --- | ---: | ---: |",
            ]
        )
        for row in projection_summary:
            lines.append(
                "| "
                + " | ".join(
                    [
                        _markdown_cell(row["base_case"]),
                        _markdown_cell(row["num_projection_cases"]),
                        _markdown_cell(row["num_ok"]),
                        _markdown_cell(row["num_export_only"]),
                        _markdown_cell(row["num_failed"]),
                        _markdown_cell(row["best_projection"]),
                        _markdown_cell(row["best_omega_total"]),
                        _markdown_cell(row["best_assessment"]),
                        _markdown_cell(row["best_ortho_minus_best_nonortho"]),
                        _markdown_cell(row["best_raw_minus_best_nonraw"]),
                    ]
                )
                + " |"
            )

    lines.extend(
        [
            "",
            "## Cases",
            "",
            (
                "| Case | Seed | Status | Selected | Classification | Best iter | "
                "Final-best | Cap | Assessment | Recommendation |"
            ),
            "| --- | --- | --- | --- | --- | ---: | ---: | --- | --- | --- |",
        ]
    )
    for row in rows:
        assessment = strategy_matrix_assessment(row, cap_tol=cap_tol)
        recommendation = strategy_matrix_recommendation(row, cap_tol=cap_tol)
        lines.append(
            "| "
            + " | ".join(
                [
                    _markdown_cell(row["case"]),
                    _markdown_cell(row["seed"]),
                    _markdown_cell(row["status"]),
                    _markdown_cell(row["selected_label"]),
                    _markdown_cell(row["selected_classification"]),
                    _markdown_cell(row["selected_best_iteration"]),
                    _markdown_cell(row["selected_final_minus_best"]),
                    _markdown_cell(row["cap_label"]),
                    _markdown_cell(assessment),
                    _markdown_cell(recommendation),
                ]
            )
            + " |"
        )
    lines.append("")
    return "\n".join(lines)


def read_strategy_matrix_payload(
    path: Path,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload.get("rows")
    if not isinstance(rows, list):
        raise ValueError(f"{path}: missing rows list")
    parsed_rows: list[dict[str, str]] = []
    for irow, row in enumerate(rows, start=1):
        if not isinstance(row, dict):
            raise ValueError(f"{path}: row {irow} is not an object")
        missing = [key for key in STRATEGY_MATRIX_HEADER if key not in row]
        if missing:
            raise ValueError(
                f"{path}: row {irow} misses field(s): {', '.join(missing)}"
            )
        parsed_rows.append({key: str(row[key]) for key in STRATEGY_MATRIX_HEADER})
    summary = payload.get("summary")
    if not isinstance(summary, dict):
        summary = {}
    return parsed_rows, dict(summary)


def read_strategy_records(
    specs: Sequence[RunSpec], *, allow_failed: bool, drift_tol: float, uphill_tol: float
) -> tuple[list[WoutHistoryRecord], list[list[str]]]:
    candidates: list[WoutHistoryRecord] = []
    failed_rows: list[list[str]] = []
    for spec in specs:
        if not _is_strategy_candidate(spec.label):
            continue
        try:
            candidates.append(
                read_wout_history_record(
                    spec, drift_tol=drift_tol, uphill_tol=uphill_tol
                )
            )
        except (OSError, ValueError) as exc:
            if not allow_failed:
                raise
            failed_rows.append(
                [
                    spec.label,
                    "no",
                    "failed",
                    "nan",
                    "na",
                    "nan",
                    "nan",
                    "nan",
                    "na",
                    "skip",
                    str(exc),
                ]
            )
    return candidates, failed_rows


def _strategy_candidate_payload(record: WoutHistoryRecord) -> dict[str, Any]:
    history = record.history
    return {
        "label": record.spec.label,
        "status": "ok",
        "classification": history.classification,
        "start_omega_total": history.start.omega_total,
        "best_iteration": history.best.iteration,
        "best_omega_total": history.best.omega_total,
        "final_omega_total": history.final_omega_total,
        "final_minus_best": history.final_minus_best,
        "first_uphill_iteration": history.first_uphill_iteration,
        "recommended_action": _strategy_action(history),
        "workdir": str(record.spec.workdir),
    }


def _strategy_failed_candidate_payload(row: Sequence[str]) -> dict[str, Any]:
    return {
        "label": row[0] if len(row) >= 1 else "unknown",
        "status": "failed",
        "error": row[-1] if row else "unknown strategy candidate failure",
    }


def select_strategy_cap_records(
    specs: Sequence[RunSpec],
    *,
    allow_failed: bool,
    drift_tol: float,
    uphill_tol: float,
    all_late_drift: bool,
) -> list[WoutHistoryRecord]:
    candidates, _ = read_strategy_records(
        specs, allow_failed=allow_failed, drift_tol=drift_tol, uphill_tol=uphill_tol
    )
    late_drift = [
        record
        for record in candidates
        if record.history.classification == "late_drift"
        and record.history.best.iteration < record.history.final.iteration
    ]
    if all_late_drift:
        return sorted(late_drift, key=lambda record: record.history.best.omega_total)
    if not late_drift:
        return []
    return [min(late_drift, key=lambda record: record.history.best.omega_total)]


def prepare_strategy_cap_run(
    seed: str,
    source_workdir: Path,
    output_workdir: Path,
    source_label: str,
    num_iter: int,
) -> None:
    if source_label == "source":
        copy_workdir(source_workdir, output_workdir, force=True)
    elif source_label in VARIANT_CHOICES:
        prepare_variant(seed, source_workdir, output_workdir, source_label, force=True)
    else:
        raise ValueError(f"cannot apply strategy cap for unknown label: {source_label}")
    set_win_num_iter(seed, output_workdir, num_iter)


def compare_strategy_cap_run(
    cap_run: StrategyCapRun,
    *,
    drift_tol: float,
    uphill_tol: float,
    match_centers: bool,
    periodic_cell: W90Cell | None,
) -> StrategyCapComparison:
    source_history = cap_run.source_record.history
    source_states = parse_wout_states(
        cap_run.source_record.spec.seed, cap_run.source_record.spec.workdir
    )
    target_state = next(
        state
        for state in source_states
        if state.iteration == source_history.best.iteration
    )
    cap_record = read_wout_history_record(
        cap_run.cap_spec, drift_tol=drift_tol, uphill_tol=uphill_tol
    )
    target_summary = wout_state_as_summary(
        cap_run.source_record.spec.seed,
        cap_run.source_record.summary.path,
        target_state,
    )
    cap_comparison = compare_summaries(
        target_summary,
        cap_record.summary,
        match_centers=match_centers,
        periodic_cell=periodic_cell,
    )
    return StrategyCapComparison(
        source_label=cap_run.source_record.spec.label,
        cap_label=cap_run.cap_spec.label,
        target_num_iter=source_history.best.iteration,
        target_best_omega_total=source_history.best.omega_total,
        target_omega_d=target_state.omega_d,
        target_omega_od=target_state.omega_od,
        cap_final_omega_total=cap_record.history.final_omega_total,
        cap_omega_d=cap_record.summary.omega_d,
        cap_omega_od=cap_record.summary.omega_od,
        cap_minus_target=(
            cap_record.history.final_omega_total - source_history.best.omega_total
        ),
        omega_d_delta=cap_record.summary.omega_d - target_state.omega_d,
        omega_od_delta=cap_record.summary.omega_od - target_state.omega_od,
        spread_delta=cap_comparison.spread_delta,
        center_delta=cap_comparison.center_delta,
        cap_classification=cap_record.history.classification,
    )


def print_history_summary(
    specs: Sequence[RunSpec],
    *,
    allow_failed: bool,
    drift_tol: float,
    uphill_tol: float,
    fail_on_drift: bool,
) -> int:
    header = [
        "label",
        "seed",
        "status",
        "classification",
        "n_iterations",
        "start_iteration",
        "start_omega_total",
        "best_iteration",
        "best_omega_total",
        "suggested_num_iter",
        "final_iteration",
        "final_omega_total",
        "final_minus_best",
        "first_uphill_iteration",
        "max_step_increase",
        "final_rms_gradient",
        "error",
    ]
    print("")
    print("# history_summary")
    print("\t".join(header))

    status = 0
    for spec in specs:
        try:
            record = read_wout_history_record(
                spec, drift_tol=drift_tol, uphill_tol=uphill_tol
            )
        except (OSError, ValueError) as exc:
            if not allow_failed:
                raise
            row = [
                spec.label,
                spec.seed,
                "failed",
                "na",
                "0",
                "na",
                "nan",
                "na",
                "nan",
                "na",
                "na",
                "nan",
                "nan",
                "na",
                "nan",
                "nan",
                str(exc),
            ]
            print("\t".join(row))
            continue

        history = record.history
        if fail_on_drift and history.classification == "late_drift":
            status = 1
        row = [
            spec.label,
            spec.seed,
            "ok",
            history.classification,
            str(history.n_iterations),
            str(history.start.iteration),
            f"{history.start.omega_total:.12e}",
            str(history.best.iteration),
            f"{history.best.omega_total:.12e}",
            str(history.best.iteration),
            str(history.final.iteration),
            f"{history.final_omega_total:.12e}",
            f"{history.final_minus_best:.12e}",
            (
                str(history.first_uphill_iteration)
                if history.first_uphill_iteration is not None
                else "none"
            ),
            f"{history.max_step_increase:.12e}",
            f"{history.final.rms_gradient:.12e}",
            "",
        ]
        print("\t".join(row))
    return status


def print_strategy_summary(
    specs: Sequence[RunSpec], *, allow_failed: bool, drift_tol: float, uphill_tol: float
) -> None:
    header = [
        "label",
        "selected",
        "classification",
        "start_omega_total",
        "best_iteration",
        "best_omega_total",
        "final_omega_total",
        "final_minus_best",
        "first_uphill_iteration",
        "recommended_action",
        "error",
    ]
    print("")
    print("# strategy_summary")
    print("\t".join(header))

    candidates: list[WoutHistoryRecord] = []
    failed_rows: list[list[str]] = []
    for spec in specs:
        if not _is_strategy_candidate(spec.label):
            continue
        try:
            candidates.append(
                read_wout_history_record(
                    spec, drift_tol=drift_tol, uphill_tol=uphill_tol
                )
            )
        except (OSError, ValueError) as exc:
            if not allow_failed:
                raise
            failed_rows.append(
                [
                    spec.label,
                    "no",
                    "failed",
                    "nan",
                    "na",
                    "nan",
                    "nan",
                    "nan",
                    "na",
                    "skip",
                    str(exc),
                ]
            )

    selected_label = ""
    if candidates:
        selected = min(candidates, key=lambda record: record.history.best.omega_total)
        selected_label = selected.spec.label
    for record in sorted(candidates, key=lambda item: item.history.best.omega_total):
        history = record.history
        row = [
            record.spec.label,
            "yes" if record.spec.label == selected_label else "no",
            history.classification,
            f"{history.start.omega_total:.12e}",
            str(history.best.iteration),
            f"{history.best.omega_total:.12e}",
            f"{history.final_omega_total:.12e}",
            f"{history.final_minus_best:.12e}",
            (
                str(history.first_uphill_iteration)
                if history.first_uphill_iteration is not None
                else "none"
            ),
            _strategy_action(history),
            "",
        ]
        print("\t".join(row))
    for row in failed_rows:
        print("\t".join(row))


def print_strategy_cap_summary(
    cap_runs: Sequence[StrategyCapRun],
    *,
    allow_failed: bool,
    drift_tol: float,
    uphill_tol: float,
    match_centers: bool,
    periodic_cell: W90Cell | None,
) -> None:
    if not cap_runs:
        return

    header = [
        "source_label",
        "cap_label",
        "status",
        "target_num_iter",
        "target_best_omega_total",
        "target_omega_d",
        "target_omega_od",
        "cap_final_omega_total",
        "cap_omega_d",
        "cap_omega_od",
        "cap_minus_target",
        "omega_d_delta",
        "omega_od_delta",
        "spread_delta",
        "center_delta",
        "cap_classification",
        "error",
    ]
    print("")
    print("# strategy_cap_summary")
    print("\t".join(header))

    for cap_run in cap_runs:
        try:
            cap_comparison = compare_strategy_cap_run(
                cap_run,
                drift_tol=drift_tol,
                uphill_tol=uphill_tol,
                match_centers=match_centers,
                periodic_cell=periodic_cell,
            )
        except (OSError, StopIteration, ValueError) as exc:
            if not allow_failed:
                raise
            row = [
                cap_run.source_record.spec.label,
                cap_run.cap_spec.label,
                "failed",
                str(cap_run.source_record.history.best.iteration),
                f"{cap_run.source_record.history.best.omega_total:.12e}",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "na",
                str(exc),
            ]
            print("\t".join(row))
            continue

        row = [
            cap_comparison.source_label,
            cap_comparison.cap_label,
            "ok",
            str(cap_comparison.target_num_iter),
            f"{cap_comparison.target_best_omega_total:.12e}",
            f"{cap_comparison.target_omega_d:.12e}",
            f"{cap_comparison.target_omega_od:.12e}",
            f"{cap_comparison.cap_final_omega_total:.12e}",
            f"{cap_comparison.cap_omega_d:.12e}",
            f"{cap_comparison.cap_omega_od:.12e}",
            f"{cap_comparison.cap_minus_target:.12e}",
            f"{cap_comparison.omega_d_delta:.12e}",
            f"{cap_comparison.omega_od_delta:.12e}",
            f"{cap_comparison.spread_delta:.12e}",
            f"{cap_comparison.center_delta:.12e}",
            cap_comparison.cap_classification,
            "",
        ]
        print("\t".join(row))


def run_strategy_case_scan(
    case: StrategyCase,
    output_workdir: Path,
    *,
    variants: Sequence[str],
    executable: str,
    timeout: float,
    run_source: bool,
    allow_failed: bool,
    run_strategy_caps: bool,
    strategy_cap_all: bool,
    drift_tol: float,
    uphill_tol: float,
) -> tuple[list[RunSpec], list[StrategyCapRun]]:
    output_workdir.mkdir(parents=True, exist_ok=True)

    if run_source:
        try:
            run_wannier90(executable, case.seed, case.source_workdir, timeout)
        except (OSError, RuntimeError, subprocess.TimeoutExpired):
            if not allow_failed:
                raise

    run_specs = []
    source_spec = RunSpec("source", case.seed, case.source_workdir)
    if run_source:
        run_specs.append(source_spec)
    else:
        try:
            read_wout_history_record(
                source_spec, drift_tol=drift_tol, uphill_tol=uphill_tol
            )
            run_specs.append(source_spec)
        except (OSError, ValueError):
            pass
    cap_runs: list[StrategyCapRun] = []
    for variant in variants:
        variant_workdir = output_workdir / variant
        try:
            prepare_variant(
                case.seed, case.source_workdir, variant_workdir, variant, force=True
            )
            run_wannier90(executable, case.seed, variant_workdir, timeout)
        except (OSError, RuntimeError, subprocess.TimeoutExpired, ValueError):
            pass
        run_specs.append(RunSpec(variant, case.seed, variant_workdir))

    if run_strategy_caps:
        for record in select_strategy_cap_records(
            run_specs,
            allow_failed=True,
            drift_tol=drift_tol,
            uphill_tol=uphill_tol,
            all_late_drift=strategy_cap_all,
        ):
            cap_iteration = record.history.best.iteration
            cap_label = f"{record.spec.label}-cap-{cap_iteration}"
            cap_workdir = output_workdir / cap_label
            prepare_strategy_cap_run(
                case.seed,
                case.source_workdir,
                cap_workdir,
                record.spec.label,
                cap_iteration,
            )
            try:
                run_wannier90(executable, case.seed, cap_workdir, timeout)
            except (OSError, RuntimeError, subprocess.TimeoutExpired):
                if not allow_failed:
                    raise
            cap_spec = RunSpec(cap_label, case.seed, cap_workdir)
            cap_runs.append(StrategyCapRun(record, cap_spec))
            run_specs.append(cap_spec)

    return run_specs, cap_runs


def command_strategy_matrix(args: argparse.Namespace) -> int:
    cases = read_strategy_cases(args.cases, args.case_file)
    variants = selected_scan_variants(
        args.variant, extended_strategy_scan=args.extended_strategy_scan
    )
    run_strategy_caps = args.run_strategy_caps or args.robust_minimization
    strategy_cap_all = args.strategy_cap_all or args.robust_minimization
    output_workdir = Path(args.output_workdir)
    if output_workdir.exists():
        if not args.force:
            raise ValueError(f"{output_workdir}: output directory already exists")
        shutil.rmtree(output_workdir)
    output_workdir.mkdir(parents=True)

    print("# strategy_matrix")
    print("\t".join(STRATEGY_MATRIX_HEADER))

    status = 0
    rows: list[dict[str, str]] = []
    case_details: list[dict[str, Any]] = []
    for case in cases:
        case_workdir = output_workdir / _safe_path_component(case.label)
        try:
            specs, cap_runs = run_strategy_case_scan(
                case,
                case_workdir,
                variants=variants,
                executable=args.executable,
                timeout=args.timeout,
                run_source=args.run_source,
                allow_failed=args.allow_failed,
                run_strategy_caps=run_strategy_caps,
                strategy_cap_all=strategy_cap_all,
                drift_tol=args.history_drift_tol,
                uphill_tol=args.history_uphill_tol,
            )
            candidates, failed_rows = read_strategy_records(
                specs,
                allow_failed=True,
                drift_tol=args.history_drift_tol,
                uphill_tol=args.history_uphill_tol,
            )
            case_detail = {
                "case": case.label,
                "seed": case.seed,
                "status": "ok",
                "candidates": [
                    _strategy_candidate_payload(record)
                    for record in sorted(
                        candidates, key=lambda item: item.history.best.omega_total
                    )
                ],
                "failed_candidates": [
                    _strategy_failed_candidate_payload(row) for row in failed_rows
                ],
            }
            if not candidates:
                if _has_wannier90_export(case.seed, case.source_workdir):
                    case_detail["status"] = "export_only"
                    row = _strategy_matrix_row(
                        [
                            case.label,
                            case.seed,
                            "export_only",
                            "na",
                            "na",
                            "na",
                            "nan",
                            "nan",
                            "nan",
                            "export_only",
                            "0",
                            str(len(failed_rows)),
                            "na",
                            "nan",
                            "nan",
                            "nan",
                            "nan",
                            "nan",
                            "na",
                            str(case_workdir),
                            "no external MLWF strategy candidates",
                        ]
                    )
                    rows.append(row)
                    case_details.append(case_detail)
                    _print_strategy_matrix_row(row)
                    continue
                raise ValueError("no usable strategy candidates")
            selected = min(
                candidates, key=lambda record: record.history.best.omega_total
            )

            cap_comparison: StrategyCapComparison | None = None
            if cap_runs:
                periodic_cell = None
                if args.periodic_centers:
                    periodic_cell = parse_win_cell(case.seed, case.source_workdir)
                cap_comparisons = [
                    compare_strategy_cap_run(
                        cap_run,
                        drift_tol=args.history_drift_tol,
                        uphill_tol=args.history_uphill_tol,
                        match_centers=bool(args.match_centers),
                        periodic_cell=periodic_cell,
                    )
                    for cap_run in cap_runs
                ]
                cap_comparison = next(
                    (
                        comparison
                        for comparison in cap_comparisons
                        if comparison.source_label == selected.spec.label
                    ),
                    cap_comparisons[0],
                )

            if (
                args.fail_on_history_drift
                and selected.history.classification == "late_drift"
            ):
                status = 1

            row = _strategy_matrix_row(
                [
                    case.label,
                    case.seed,
                    "ok",
                    selected.spec.label,
                    selected.history.classification,
                    str(selected.history.best.iteration),
                    f"{selected.history.best.omega_total:.12e}",
                    f"{selected.history.final_omega_total:.12e}",
                    f"{selected.history.final_minus_best:.12e}",
                    _strategy_action(selected.history),
                    str(len(candidates)),
                    str(len(failed_rows)),
                    cap_comparison.cap_label if cap_comparison is not None else "na",
                    (
                        f"{cap_comparison.cap_minus_target:.12e}"
                        if cap_comparison is not None
                        else "nan"
                    ),
                    (
                        f"{cap_comparison.omega_d_delta:.12e}"
                        if cap_comparison is not None
                        else "nan"
                    ),
                    (
                        f"{cap_comparison.omega_od_delta:.12e}"
                        if cap_comparison is not None
                        else "nan"
                    ),
                    (
                        f"{cap_comparison.spread_delta:.12e}"
                        if cap_comparison is not None
                        else "nan"
                    ),
                    (
                        f"{cap_comparison.center_delta:.12e}"
                        if cap_comparison is not None
                        else "nan"
                    ),
                    (
                        cap_comparison.cap_classification
                        if cap_comparison is not None
                        else "na"
                    ),
                    str(case_workdir),
                    "",
                ]
            )
            if (
                args.fail_on_uncapped_drift
                and selected.history.classification == "late_drift"
                and not _strategy_cap_verified(row, args.cap_tol)
            ):
                status = 1
            case_details.append(case_detail)
        except (
            OSError,
            StopIteration,
            RuntimeError,
            subprocess.TimeoutExpired,
            ValueError,
        ) as exc:
            if not args.allow_failed:
                raise
            status = 1
            row = _strategy_matrix_row(
                [
                    case.label,
                    case.seed,
                    "failed",
                    "na",
                    "na",
                    "na",
                    "nan",
                    "nan",
                    "nan",
                    "na",
                    "0",
                    "0",
                    "na",
                    "nan",
                    "nan",
                    "nan",
                    "nan",
                    "nan",
                    "na",
                    str(case_workdir),
                    str(exc),
                ]
            )
            case_details.append(
                {
                    "case": case.label,
                    "seed": case.seed,
                    "status": "failed",
                    "candidates": [],
                    "failed_candidates": [],
                    "error": str(exc),
                }
            )
        rows.append(row)
        _print_strategy_matrix_row(row)
    summary = summarize_strategy_matrix(rows, cap_tol=args.cap_tol)
    projection_summary = summarize_strategy_projection_families(
        rows, cap_tol=args.cap_tol
    )
    if not args.no_summary:
        print_strategy_matrix_summary(summary)
        print_strategy_projection_family_summary(projection_summary)
    if args.json_output is not None:
        Path(args.json_output).write_text(
            json.dumps(
                {
                    "rows": rows,
                    "summary": summary,
                    "projection_summary": projection_summary,
                    "assessments": strategy_matrix_assessments(
                        rows, cap_tol=args.cap_tol
                    ),
                    "case_details": case_details,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    if args.markdown_output is not None:
        Path(args.markdown_output).write_text(
            strategy_matrix_markdown(
                rows, summary, projection_summary, cap_tol=args.cap_tol
            ),
            encoding="utf-8",
        )
    return status


def command_cp2k_strategy_matrix(args: argparse.Namespace) -> int:
    cases = read_cp2k_strategy_cases(args.cases, args.case_file)
    output_workdir = Path(args.output_workdir)
    if output_workdir.exists():
        if args.no_cp2k_run:
            matrix_output_workdir = output_workdir / "strategy_matrix"
            if matrix_output_workdir.exists():
                if not args.force:
                    raise ValueError(
                        f"{matrix_output_workdir}: output directory already exists"
                    )
                shutil.rmtree(matrix_output_workdir)
        elif not args.force:
            raise ValueError(f"{output_workdir}: output directory already exists")
        else:
            shutil.rmtree(output_workdir)
    output_workdir.mkdir(parents=True, exist_ok=args.no_cp2k_run)

    cases = expand_cp2k_projection_cases(
        cases,
        output_workdir / "cp2k_inputs",
        args.cp2k_initial_projection,
        args.cp2k_projection_smoothing,
        args.cp2k_auto_projection_score,
        write_inputs=not args.no_cp2k_run,
    )

    source_root = output_workdir / "cp2k_exports"
    strategy_cases: list[str] = []
    for case in cases:
        source_workdir = source_root / _safe_path_component(case.label)
        if not args.no_cp2k_run:
            try:
                run_cp2k_export(
                    args.cp2k_executable,
                    case,
                    source_workdir,
                    timeout=args.cp2k_timeout,
                    copy_input_dir_files=args.copy_input_dir_files,
                )
            except (OSError, RuntimeError, subprocess.TimeoutExpired, ValueError):
                if not args.allow_failed:
                    raise
        seed = infer_cp2k_wannier_seed(case)
        strategy_cases.append(f"{case.label}:{seed}:{source_workdir}")

    matrix_output_workdir = output_workdir / "strategy_matrix"
    matrix_args = argparse.Namespace(
        cases=strategy_cases,
        case_file=None,
        output_workdir=str(matrix_output_workdir),
        variant=args.variant,
        extended_strategy_scan=args.extended_strategy_scan,
        run_strategy_caps=args.run_strategy_caps or args.robust_minimization,
        strategy_cap_all=args.strategy_cap_all or args.robust_minimization,
        robust_minimization=False,
        executable=args.executable,
        timeout=args.timeout,
        force=False,
        allow_failed=args.allow_failed,
        run_source=not args.no_run_source,
        history_drift_tol=args.history_drift_tol,
        history_uphill_tol=args.history_uphill_tol,
        fail_on_history_drift=args.fail_on_history_drift,
        fail_on_uncapped_drift=args.fail_on_uncapped_drift,
        cap_tol=args.cap_tol,
        no_summary=args.no_summary,
        json_output=args.json_output,
        markdown_output=args.markdown_output,
        match_centers=args.match_centers,
        periodic_centers=args.periodic_centers,
    )
    return command_strategy_matrix(matrix_args)


def _selected_auto_projection(line: str) -> str:
    match = re.search(r"AUTO_SCDM_RAW selected\s+(\S+)", line)
    return match.group(1) if match is not None else ""


def _read_auto_selection(seed: str, workdir: Path) -> tuple[str, str]:
    selected = ""
    margin = ""
    for path in sorted(workdir.glob("*.out")):
        text = path.read_text(encoding="utf-8", errors="replace")
        for line in text.splitlines():
            if "AUTO_SCDM_RAW selected" in line:
                selected = line.strip()
            elif "AUTO_SCDM_RAW best-minus-runner-up score margin" in line:
                fields = line.split()
                if fields:
                    margin = fields[-1]
    if not selected:
        cp2k_out = workdir / f"{seed}.out"
        if cp2k_out.exists():
            text = cp2k_out.read_text(encoding="utf-8", errors="replace")
            for line in text.splitlines():
                if "AUTO_SCDM_RAW selected" in line:
                    selected = line.strip()
                elif "AUTO_SCDM_RAW best-minus-runner-up score margin" in line:
                    fields = line.split()
                    if fields:
                        margin = fields[-1]
    return selected, margin


def _read_auto_candidate_scores(workdir: Path) -> list[dict[str, str]]:
    pattern = re.compile(
        r"AUTO_SCDM_RAW candidate\s+(\S+)\s+"
        r"total/conditioning/spread-proxy(?:/Mmn-proxy(?:/locality-prior)?)? scores\s+(.+)"
    )
    rows: list[dict[str, str]] = []
    for path in sorted(workdir.glob("*.out")):
        text = path.read_text(encoding="utf-8", errors="replace")
        for line in text.splitlines():
            match = pattern.search(line)
            if match is None:
                continue
            values = match.group(2).split()
            rows.append(
                {
                    "candidate_projection": match.group(1),
                    "total_score": values[0] if len(values) > 0 else "nan",
                    "conditioning_score": values[1] if len(values) > 1 else "nan",
                    "spread_proxy_score": values[2] if len(values) > 2 else "nan",
                    "mmn_proxy_score": values[3] if len(values) > 3 else "nan",
                    "locality_prior_score": values[4] if len(values) > 4 else "nan",
                    "source_file": str(path),
                }
            )
    return rows


def _projection_matrix_summary(rows: Sequence[dict[str, str]]) -> list[dict[str, str]]:
    summaries: list[dict[str, str]] = []
    cases = sorted({row["case"] for row in rows})
    for case in cases:
        case_rows = [
            row for row in rows if row["case"] == case and row["status"] == "ok"
        ]
        direct_rows = [
            row for row in case_rows if row["initial_projection"] != "AUTO_SCDM_RAW"
        ]
        auto_rows = [
            row for row in case_rows if row["initial_projection"] == "AUTO_SCDM_RAW"
        ]
        if not direct_rows and not auto_rows:
            continue
        best_direct = (
            min(direct_rows, key=lambda row: float(row["omega_total"]))
            if direct_rows
            else None
        )
        best_auto = (
            min(auto_rows, key=lambda row: float(row["omega_total"]))
            if auto_rows
            else None
        )
        auto_regret = "nan"
        if best_direct is not None and best_auto is not None:
            auto_regret = f"{float(best_auto['omega_total']) - float(best_direct['omega_total']):.12e}"
        summaries.append(
            {
                "case": case,
                "best_direct_projection": (
                    best_direct["initial_projection"]
                    if best_direct is not None
                    else "na"
                ),
                "best_direct_omega_total": (
                    best_direct["omega_total"] if best_direct is not None else "nan"
                ),
                "best_auto_score": (
                    best_auto["auto_projection_score"]
                    if best_auto is not None
                    else "na"
                ),
                "best_auto_selected_projection": (
                    best_auto["auto_selected_projection"]
                    if best_auto is not None
                    else "na"
                ),
                "best_auto_omega_total": (
                    best_auto["omega_total"] if best_auto is not None else "nan"
                ),
                "best_auto_regret": auto_regret,
            }
        )
    return summaries


def _projection_score_mode_summary(
    rows: Sequence[dict[str, str]], score_rows: Sequence[dict[str, str]]
) -> list[dict[str, str]]:
    summaries: list[dict[str, str]] = []
    cases = sorted({row["case"] for row in rows})
    for case in cases:
        direct_rows = [
            row
            for row in rows
            if row["case"] == case
            and row["status"] == "ok"
            and row["initial_projection"] != "AUTO_SCDM_RAW"
        ]
        if not direct_rows:
            continue
        best_direct = min(direct_rows, key=lambda row: float(row["omega_total"]))
        auto_rows = [
            row
            for row in rows
            if row["case"] == case
            and row["status"] == "ok"
            and row["initial_projection"] == "AUTO_SCDM_RAW"
        ]
        for auto_row in sorted(auto_rows, key=lambda row: row["auto_projection_score"]):
            score_mode = auto_row["auto_projection_score"]
            candidate_scores = [
                row
                for row in score_rows
                if row["case"] == case and row["auto_projection_score"] == score_mode
            ]
            score_selected = ""
            best_direct_score_rank = ""
            if candidate_scores:
                ranked_scores = sorted(
                    candidate_scores,
                    key=lambda row: float(row["total_score"]),
                    reverse=True,
                )
                score_selected = ranked_scores[0]["candidate_projection"]
                for rank, score_row in enumerate(ranked_scores, start=1):
                    if (
                        score_row["candidate_projection"]
                        == best_direct["initial_projection"]
                    ):
                        best_direct_score_rank = str(rank)
                        break
            regret = float(auto_row["omega_total"]) - float(best_direct["omega_total"])
            summaries.append(
                {
                    "case": case,
                    "auto_projection_score": score_mode,
                    "auto_selected_projection": auto_row["auto_selected_projection"],
                    "score_selected_projection": score_selected,
                    "auto_omega_total": auto_row["omega_total"],
                    "best_direct_projection": best_direct["initial_projection"],
                    "best_direct_omega_total": best_direct["omega_total"],
                    "best_direct_score_rank": best_direct_score_rank,
                    "auto_regret": f"{regret:.12e}",
                    "matches_best_direct": (
                        "yes"
                        if auto_row["auto_selected_projection"]
                        == best_direct["initial_projection"]
                        else "no"
                    ),
                }
            )
    return summaries


def _projection_early_minimization_rows(
    rows: Sequence[dict[str, str]],
    output_workdir: Path,
    iterations: Sequence[int],
    *,
    executable: str,
    timeout: float,
    allow_failed: bool,
    drift_tol: float,
    uphill_tol: float,
) -> list[dict[str, str]]:
    early_rows: list[dict[str, str]] = []
    direct_rows = [
        row
        for row in rows
        if row["status"] == "ok"
        and row["initial_projection"] != "AUTO_SCDM_RAW"
        and row["workdir"]
    ]
    if not direct_rows:
        return early_rows

    early_root = output_workdir / "early_minimization"
    for iteration in sorted(set(iterations)):
        if iteration < 0:
            raise ValueError("--early-minimization-iteration must be non-negative")
        for row in direct_rows:
            early_workdir = (
                early_root
                / f"iter_{iteration}"
                / _safe_path_component(f"{row['case']}_{row['initial_projection']}")
            )
            status = "ok"
            error = ""
            omega_total = omega_i = omega_d = omega_od = "nan"
            best_iteration = final_iteration = "na"
            best_omega_total = final_omega_total = final_minus_best = "nan"
            classification = "na"
            first_uphill_iteration = "na"
            max_step_increase = final_rms_gradient = "nan"
            all_done = "no"
            try:
                copy_workdir(Path(row["workdir"]), early_workdir, force=True)
                for suffix in (
                    "wout",
                    "chk",
                    "chk.fmt",
                    "r2mn",
                    "uHu",
                    "uIu",
                    "spn",
                    "wannier90.stdout",
                ):
                    stale_output = early_workdir / f"{row['seed']}.{suffix}"
                    if stale_output.exists():
                        stale_output.unlink()
                set_win_num_iter(row["seed"], early_workdir, iteration)
                run_wannier90(executable, row["seed"], early_workdir, timeout)
                summary = parse_wout(row["seed"], early_workdir)
                history = summarize_wout_history(
                    parse_wout_iterations(row["seed"], early_workdir),
                    summary.omega_total,
                    drift_tol=drift_tol,
                    uphill_tol=uphill_tol,
                )
                omega_i = f"{summary.omega_i:.12e}"
                omega_d = f"{summary.omega_d:.12e}"
                omega_od = f"{summary.omega_od:.12e}"
                omega_total = f"{summary.omega_total:.12e}"
                all_done = "yes" if summary.all_done else "no"
                best_iteration = str(history.best.iteration)
                best_omega_total = f"{history.best.omega_total:.12e}"
                final_iteration = str(history.final.iteration)
                final_omega_total = f"{history.final_omega_total:.12e}"
                final_minus_best = f"{history.final_minus_best:.12e}"
                classification = history.classification
                first_uphill_iteration = (
                    str(history.first_uphill_iteration)
                    if history.first_uphill_iteration is not None
                    else "none"
                )
                max_step_increase = f"{history.max_step_increase:.12e}"
                final_rms_gradient = f"{history.final.rms_gradient:.12e}"
            except (
                OSError,
                RuntimeError,
                subprocess.TimeoutExpired,
                ValueError,
            ) as exc:
                status = "failed"
                error = str(exc)
                if not allow_failed:
                    raise
            early_rows.append(
                {
                    "case": row["case"],
                    "iteration_budget": str(iteration),
                    "initial_projection": row["initial_projection"],
                    "status": status,
                    "seed": row["seed"],
                    "omega_total": omega_total,
                    "omega_i": omega_i,
                    "omega_d": omega_d,
                    "omega_od": omega_od,
                    "all_done": all_done,
                    "best_iteration": best_iteration,
                    "best_omega_total": best_omega_total,
                    "final_iteration": final_iteration,
                    "final_omega_total": final_omega_total,
                    "final_minus_best": final_minus_best,
                    "classification": classification,
                    "first_uphill_iteration": first_uphill_iteration,
                    "max_step_increase": max_step_increase,
                    "final_rms_gradient": final_rms_gradient,
                    "source_workdir": row["workdir"],
                    "workdir": str(early_workdir),
                    "error": error,
                }
            )
    return early_rows


def _projection_early_minimization_summary(
    rows: Sequence[dict[str, str]], early_rows: Sequence[dict[str, str]]
) -> list[dict[str, str]]:
    summaries: list[dict[str, str]] = []
    cases = sorted({row["case"] for row in early_rows})
    for case in cases:
        direct_rows = [
            row
            for row in rows
            if row["case"] == case
            and row["status"] == "ok"
            and row["initial_projection"] != "AUTO_SCDM_RAW"
        ]
        auto_rows = [
            row
            for row in rows
            if row["case"] == case
            and row["status"] == "ok"
            and row["initial_projection"] == "AUTO_SCDM_RAW"
        ]
        best_direct = (
            min(direct_rows, key=lambda row: float(row["omega_total"]))
            if direct_rows
            else None
        )
        best_auto = (
            min(auto_rows, key=lambda row: float(row["omega_total"]))
            if auto_rows
            else None
        )
        for iteration in sorted(
            {row["iteration_budget"] for row in early_rows if row["case"] == case},
            key=int,
        ):
            case_rows = [
                row
                for row in early_rows
                if row["case"] == case
                and row["iteration_budget"] == iteration
                and row["status"] == "ok"
            ]
            failed_rows = [
                row
                for row in early_rows
                if row["case"] == case
                and row["iteration_budget"] == iteration
                and row["status"] != "ok"
            ]
            if not case_rows:
                summaries.append(
                    {
                        "case": case,
                        "iteration_budget": iteration,
                        "selected_projection": "na",
                        "selected_best_iteration": "na",
                        "selected_best_omega_total": "nan",
                        "selected_final_omega_total": "nan",
                        "selected_classification": "na",
                        "best_direct_projection": (
                            best_direct["initial_projection"]
                            if best_direct is not None
                            else "na"
                        ),
                        "best_direct_omega_total": (
                            best_direct["omega_total"]
                            if best_direct is not None
                            else "nan"
                        ),
                        "auto_selected_projection": (
                            best_auto["auto_selected_projection"]
                            if best_auto is not None
                            else "na"
                        ),
                        "auto_omega_total": (
                            best_auto["omega_total"] if best_auto is not None else "nan"
                        ),
                        "matches_best_direct": "no",
                        "matches_auto_selected": "no",
                        "num_candidates": "0",
                        "num_failed_candidates": str(len(failed_rows)),
                    }
                )
                continue
            selected = min(case_rows, key=lambda row: float(row["best_omega_total"]))
            summaries.append(
                {
                    "case": case,
                    "iteration_budget": iteration,
                    "selected_projection": selected["initial_projection"],
                    "selected_best_iteration": selected["best_iteration"],
                    "selected_best_omega_total": selected["best_omega_total"],
                    "selected_final_omega_total": selected["final_omega_total"],
                    "selected_classification": selected["classification"],
                    "best_direct_projection": (
                        best_direct["initial_projection"]
                        if best_direct is not None
                        else "na"
                    ),
                    "best_direct_omega_total": (
                        best_direct["omega_total"] if best_direct is not None else "nan"
                    ),
                    "auto_selected_projection": (
                        best_auto["auto_selected_projection"]
                        if best_auto is not None
                        else "na"
                    ),
                    "auto_omega_total": (
                        best_auto["omega_total"] if best_auto is not None else "nan"
                    ),
                    "matches_best_direct": (
                        "yes"
                        if best_direct is not None
                        and selected["initial_projection"]
                        == best_direct["initial_projection"]
                        else "no"
                    ),
                    "matches_auto_selected": (
                        "yes"
                        if best_auto is not None
                        and selected["initial_projection"]
                        == best_auto["auto_selected_projection"]
                        else "no"
                    ),
                    "num_candidates": str(len(case_rows)),
                    "num_failed_candidates": str(len(failed_rows)),
                }
            )
    return summaries


def _projection_robust_minimization_summary(
    rows: Sequence[dict[str, str]], early_rows: Sequence[dict[str, str]]
) -> list[dict[str, str]]:
    summaries: list[dict[str, str]] = []
    cases = sorted({row["case"] for row in early_rows})
    for case in cases:
        case_rows = [
            row for row in early_rows if row["case"] == case and row["status"] == "ok"
        ]
        failed_rows = [
            row for row in early_rows if row["case"] == case and row["status"] != "ok"
        ]
        auto_rows = [
            row
            for row in rows
            if row["case"] == case
            and row["status"] == "ok"
            and row["initial_projection"] == "AUTO_SCDM_RAW"
        ]
        auto_selected = auto_rows[0]["auto_selected_projection"] if auto_rows else "na"
        stable_rows = [row for row in case_rows if row["classification"] == "stable"]
        if stable_rows:
            selected = min(stable_rows, key=lambda row: float(row["final_omega_total"]))
            lower_late_rows = [
                row
                for row in case_rows
                if row["classification"] == "late_drift"
                and float(row["best_omega_total"])
                < float(selected["final_omega_total"]) - 1.0e-6
            ]
            confidence = "high"
            if lower_late_rows:
                confidence = "medium"
            stable_selected_count = sum(
                1
                for row in stable_rows
                if row["initial_projection"] == selected["initial_projection"]
            )
            if stable_selected_count < 2:
                confidence = "medium" if confidence == "high" else confidence
            selection_basis = "stable_final_omega"
        elif case_rows:
            selected = min(case_rows, key=lambda row: float(row["best_omega_total"]))
            lower_late_rows = []
            stable_selected_count = 0
            confidence = "low"
            selection_basis = "late_drift_best_omega"
        else:
            summaries.append(
                {
                    "case": case,
                    "robust_projection": "na",
                    "confidence": "low",
                    "selection_basis": "no_successful_candidates",
                    "selected_iteration_budget": "na",
                    "selected_best_iteration": "na",
                    "selected_best_omega_total": "nan",
                    "selected_final_omega_total": "nan",
                    "selected_classification": "na",
                    "auto_selected_projection": auto_selected,
                    "matches_auto_selected": "no",
                    "num_stable_candidates": "0",
                    "num_late_drift_candidates": "0",
                    "num_failed_candidates": str(len(failed_rows)),
                    "num_lower_late_drift_candidates": "0",
                    "stable_selected_count": "0",
                    "all_iteration_budgets": "",
                }
            )
            continue
        summaries.append(
            {
                "case": case,
                "robust_projection": selected["initial_projection"],
                "confidence": confidence,
                "selection_basis": selection_basis,
                "selected_iteration_budget": selected["iteration_budget"],
                "selected_best_iteration": selected["best_iteration"],
                "selected_best_omega_total": selected["best_omega_total"],
                "selected_final_omega_total": selected["final_omega_total"],
                "selected_classification": selected["classification"],
                "auto_selected_projection": auto_selected,
                "matches_auto_selected": (
                    "yes" if selected["initial_projection"] == auto_selected else "no"
                ),
                "num_stable_candidates": str(len(stable_rows)),
                "num_late_drift_candidates": str(
                    sum(1 for row in case_rows if row["classification"] == "late_drift")
                ),
                "num_failed_candidates": str(len(failed_rows)),
                "num_lower_late_drift_candidates": str(len(lower_late_rows)),
                "stable_selected_count": str(stable_selected_count),
                "all_iteration_budgets": ",".join(
                    sorted({row["iteration_budget"] for row in case_rows}, key=int)
                ),
            }
        )
    return summaries


def _write_tsv(
    path: Path, rows: Sequence[dict[str, str]], header: Sequence[str]
) -> None:
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join(row.get(column, "") for column in header))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _projection_matrix_markdown(
    rows: Sequence[dict[str, str]],
    summary: Sequence[dict[str, str]],
    score_mode_summary: Sequence[dict[str, str]],
    early_summary: Sequence[dict[str, str]],
    robust_summary: Sequence[dict[str, str]],
) -> str:
    lines = [
        "# CP2K Wannier90 Projection Matrix",
        "",
        "## Summary",
        "",
        "| case | best direct | best direct Omega | best auto score | auto selected | auto Omega | regret |",
        "| --- | --- | ---: | --- | --- | ---: | ---: |",
    ]
    for row in summary:
        lines.append(
            "| {case} | {best_direct_projection} | {best_direct_omega_total} | "
            "{best_auto_score} | {best_auto_selected_projection} | "
            "{best_auto_omega_total} | {best_auto_regret} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Runs",
            "",
            "| case | projection | score | smoothing | status | Omega | selected | margin |",
            "| --- | --- | --- | --- | --- | ---: | --- | ---: |",
        ]
    )
    for row in rows:
        lines.append(
            "| {case} | {initial_projection} | {auto_projection_score} | "
            "{projection_smoothing} | {status} | {omega_total} | "
            "{auto_selected_projection} | {auto_score_margin} |".format(**row)
        )
    if score_mode_summary:
        lines.extend(
            [
                "",
                "## Auto Score Modes",
                "",
                "| case | score | auto selected | score selected | best direct | rank | regret |",
                "| --- | --- | --- | --- | --- | ---: | ---: |",
            ]
        )
        for row in score_mode_summary:
            lines.append(
                "| {case} | {auto_projection_score} | {auto_selected_projection} | "
                "{score_selected_projection} | {best_direct_projection} | "
                "{best_direct_score_rank} | {auto_regret} |".format(**row)
            )
    if early_summary:
        lines.extend(
            [
                "",
                "## Early Minimization Selection",
                "",
                "| case | num_iter | selected | best iter | best Omega | final Omega | class | auto selected | static best | candidates | failed |",
                "| --- | ---: | --- | ---: | ---: | ---: | --- | --- | --- | ---: | ---: |",
            ]
        )
        for row in early_summary:
            lines.append(
                "| {case} | {iteration_budget} | {selected_projection} | "
                "{selected_best_iteration} | {selected_best_omega_total} | "
                "{selected_final_omega_total} | {selected_classification} | "
                "{auto_selected_projection} | {best_direct_projection} | "
                "{num_candidates} | {num_failed_candidates} |".format(**row)
            )
    if robust_summary:
        lines.extend(
            [
                "",
                "## Robust Early-Minimization Selection",
                "",
                "| case | robust projection | confidence | basis | num_iter | best iter | best Omega | final Omega | class | auto selected | stable | late drift | lower late drift |",
                "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: |",
            ]
        )
        for row in robust_summary:
            lines.append(
                "| {case} | {robust_projection} | {confidence} | {selection_basis} | "
                "{selected_iteration_budget} | {selected_best_iteration} | "
                "{selected_best_omega_total} | {selected_final_omega_total} | "
                "{selected_classification} | {auto_selected_projection} | "
                "{num_stable_candidates} | {num_late_drift_candidates} | "
                "{num_lower_late_drift_candidates} |".format(**row)
            )
    return "\n".join(lines) + "\n"


def command_cp2k_projection_matrix(args: argparse.Namespace) -> int:
    base_cases = read_cp2k_strategy_cases(args.cases, args.case_file)
    output_workdir = Path(args.output_workdir)
    if output_workdir.exists():
        if args.no_cp2k_run:
            pass
        elif not args.force:
            raise ValueError(f"{output_workdir}: output directory already exists")
        else:
            shutil.rmtree(output_workdir)
    output_workdir.mkdir(parents=True, exist_ok=args.no_cp2k_run)

    projections = tuple(
        args.cp2k_initial_projection or DEFAULT_CP2K_PROJECTION_MATRIX_PROJECTIONS
    )
    auto_scores = tuple(args.cp2k_auto_projection_score or CP2K_AUTO_PROJECTION_SCORES)
    smoothings = tuple(args.cp2k_projection_smoothing or ("MMN_TRANSPORT",))

    rows: list[dict[str, str]] = []
    score_rows: list[dict[str, str]] = []
    expanded_input_root = output_workdir / "cp2k_inputs"
    export_root = output_workdir / "cp2k_exports"
    for base_case in base_cases:
        expanded_cases = expand_cp2k_projection_cases(
            [base_case],
            expanded_input_root,
            projections,
            smoothings,
            auto_scores if "AUTO_SCDM_RAW" in projections else None,
            write_inputs=True,
        )
        for case in expanded_cases:
            source_workdir = export_root / _safe_path_component(case.label)
            seed = infer_cp2k_wannier_seed(case)
            projection = _cp2k_input_value(case.input_path, "INITIAL_PROJECTIONS") or ""
            auto_score = (
                _cp2k_input_value(case.input_path, "AUTO_PROJECTION_SCORE") or ""
            )
            smoothing = _cp2k_input_value(case.input_path, "PROJECTION_SMOOTHING") or ""
            cp2k_status = "ok"
            wannier90_status = "ok"
            error = ""
            if not args.no_cp2k_run:
                try:
                    run_cp2k_export(
                        args.cp2k_executable,
                        case,
                        source_workdir,
                        timeout=args.cp2k_timeout,
                        copy_input_dir_files=args.copy_input_dir_files,
                    )
                except (
                    OSError,
                    RuntimeError,
                    subprocess.TimeoutExpired,
                    ValueError,
                ) as exc:
                    cp2k_status = "failed"
                    wannier90_status = "not-run"
                    error = str(exc)
                    if not args.allow_failed:
                        raise
            if cp2k_status == "ok" and not args.no_wannier90_run:
                try:
                    run_wannier90(args.executable, seed, source_workdir, args.timeout)
                except (OSError, RuntimeError, subprocess.TimeoutExpired) as exc:
                    wannier90_status = "failed"
                    error = str(exc)
                    if not args.allow_failed:
                        raise
            status = (
                "ok" if cp2k_status == "ok" and wannier90_status == "ok" else "failed"
            )
            omega_i = omega_d = omega_od = omega_total = "nan"
            all_done = "no"
            if wannier90_status == "ok":
                try:
                    summary = parse_wout(seed, source_workdir)
                    omega_i = f"{summary.omega_i:.12e}"
                    omega_d = f"{summary.omega_d:.12e}"
                    omega_od = f"{summary.omega_od:.12e}"
                    omega_total = f"{summary.omega_total:.12e}"
                    all_done = "yes" if summary.all_done else "no"
                except (OSError, ValueError) as exc:
                    status = "failed"
                    error = str(exc)
                    if not args.allow_failed:
                        raise
            auto_selected_line, auto_margin = _read_auto_selection(seed, source_workdir)
            rows.append(
                {
                    "case": base_case.label,
                    "label": case.label,
                    "seed": seed,
                    "initial_projection": projection,
                    "auto_projection_score": auto_score,
                    "projection_smoothing": smoothing,
                    "status": status,
                    "cp2k_status": cp2k_status,
                    "wannier90_status": wannier90_status,
                    "omega_total": omega_total,
                    "omega_i": omega_i,
                    "omega_d": omega_d,
                    "omega_od": omega_od,
                    "all_done": all_done,
                    "auto_selected_projection": _selected_auto_projection(
                        auto_selected_line
                    ),
                    "auto_score_margin": auto_margin,
                    "auto_selected_line": auto_selected_line,
                    "workdir": str(source_workdir),
                    "error": error,
                }
            )
            for score_row in _read_auto_candidate_scores(source_workdir):
                score_rows.append(
                    {
                        "case": base_case.label,
                        "label": case.label,
                        "seed": seed,
                        "auto_projection_score": auto_score,
                        "projection_smoothing": smoothing,
                        **score_row,
                    }
                )

    header = [
        "case",
        "label",
        "seed",
        "initial_projection",
        "auto_projection_score",
        "projection_smoothing",
        "status",
        "cp2k_status",
        "wannier90_status",
        "omega_total",
        "omega_i",
        "omega_d",
        "omega_od",
        "all_done",
        "auto_selected_projection",
        "auto_score_margin",
        "auto_selected_line",
        "workdir",
        "error",
    ]
    score_header = [
        "case",
        "label",
        "seed",
        "auto_projection_score",
        "projection_smoothing",
        "candidate_projection",
        "total_score",
        "conditioning_score",
        "spread_proxy_score",
        "mmn_proxy_score",
        "locality_prior_score",
        "source_file",
    ]
    summary = _projection_matrix_summary(rows)
    summary_header = [
        "case",
        "best_direct_projection",
        "best_direct_omega_total",
        "best_auto_score",
        "best_auto_selected_projection",
        "best_auto_omega_total",
        "best_auto_regret",
    ]
    score_mode_summary = _projection_score_mode_summary(rows, score_rows)
    score_mode_header = [
        "case",
        "auto_projection_score",
        "auto_selected_projection",
        "score_selected_projection",
        "auto_omega_total",
        "best_direct_projection",
        "best_direct_omega_total",
        "best_direct_score_rank",
        "auto_regret",
        "matches_best_direct",
    ]
    early_rows: list[dict[str, str]] = []
    early_summary: list[dict[str, str]] = []
    robust_summary: list[dict[str, str]] = []
    early_header = [
        "case",
        "iteration_budget",
        "initial_projection",
        "status",
        "seed",
        "omega_total",
        "omega_i",
        "omega_d",
        "omega_od",
        "all_done",
        "best_iteration",
        "best_omega_total",
        "final_iteration",
        "final_omega_total",
        "final_minus_best",
        "classification",
        "first_uphill_iteration",
        "max_step_increase",
        "final_rms_gradient",
        "source_workdir",
        "workdir",
        "error",
    ]
    early_summary_header = [
        "case",
        "iteration_budget",
        "selected_projection",
        "selected_best_iteration",
        "selected_best_omega_total",
        "selected_final_omega_total",
        "selected_classification",
        "best_direct_projection",
        "best_direct_omega_total",
        "auto_selected_projection",
        "auto_omega_total",
        "matches_best_direct",
        "matches_auto_selected",
        "num_candidates",
        "num_failed_candidates",
    ]
    robust_summary_header = [
        "case",
        "robust_projection",
        "confidence",
        "selection_basis",
        "selected_iteration_budget",
        "selected_best_iteration",
        "selected_best_omega_total",
        "selected_final_omega_total",
        "selected_classification",
        "auto_selected_projection",
        "matches_auto_selected",
        "num_stable_candidates",
        "num_late_drift_candidates",
        "num_failed_candidates",
        "num_lower_late_drift_candidates",
        "stable_selected_count",
        "all_iteration_budgets",
    ]
    if args.early_minimization_iteration:
        early_rows = _projection_early_minimization_rows(
            rows,
            output_workdir,
            args.early_minimization_iteration,
            executable=args.executable,
            timeout=args.timeout,
            allow_failed=args.allow_failed,
            drift_tol=args.history_drift_tol,
            uphill_tol=args.history_uphill_tol,
        )
        early_summary = _projection_early_minimization_summary(rows, early_rows)
        robust_summary = _projection_robust_minimization_summary(rows, early_rows)
    _write_tsv(output_workdir / "projection_matrix.tsv", rows, header)
    _write_tsv(
        output_workdir / "projection_matrix_summary.tsv", summary, summary_header
    )
    _write_tsv(output_workdir / "projection_score_matrix.tsv", score_rows, score_header)
    _write_tsv(
        output_workdir / "projection_score_mode_summary.tsv",
        score_mode_summary,
        score_mode_header,
    )
    _write_tsv(
        output_workdir / "projection_early_minimization.tsv", early_rows, early_header
    )
    _write_tsv(
        output_workdir / "projection_early_minimization_summary.tsv",
        early_summary,
        early_summary_header,
    )
    _write_tsv(
        output_workdir / "projection_robust_minimization_summary.tsv",
        robust_summary,
        robust_summary_header,
    )
    (output_workdir / "projection_matrix.json").write_text(
        json.dumps(
            {
                "rows": rows,
                "summary": summary,
                "score_rows": score_rows,
                "score_mode_summary": score_mode_summary,
                "early_minimization_rows": early_rows,
                "early_minimization_summary": early_summary,
                "robust_minimization_summary": robust_summary,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    markdown = _projection_matrix_markdown(
        rows, summary, score_mode_summary, early_summary, robust_summary
    )
    (output_workdir / "projection_matrix.md").write_text(markdown, encoding="utf-8")
    if args.markdown_output is not None:
        Path(args.markdown_output).write_text(markdown, encoding="utf-8")
    if not args.no_summary:
        print("# cp2k_projection_matrix_summary")
        print("\t".join(summary_header))
        for row in summary:
            print("\t".join(row[column] for column in summary_header))
        if early_summary:
            print("")
            print("# cp2k_projection_early_minimization_summary")
            print("\t".join(early_summary_header))
            for row in early_summary:
                print("\t".join(row[column] for column in early_summary_header))
        if robust_summary:
            print("")
            print("# cp2k_projection_robust_minimization_summary")
            print("\t".join(robust_summary_header))
            for row in robust_summary:
                print("\t".join(row[column] for column in robust_summary_header))
    return 0


def command_strategy_report(args: argparse.Namespace) -> int:
    rows, _ = read_strategy_matrix_payload(Path(args.json_input))
    summary = summarize_strategy_matrix(rows, cap_tol=args.cap_tol)
    projection_summary = summarize_strategy_projection_families(
        rows, cap_tol=args.cap_tol
    )
    report = strategy_matrix_markdown(
        rows, summary, projection_summary, cap_tol=args.cap_tol
    )
    if args.output is None:
        print(report, end="")
    else:
        Path(args.output).write_text(report, encoding="utf-8")
    return 0


def command_scan_variants(args: argparse.Namespace) -> int:
    source_workdir = Path(args.source_workdir)
    output_workdir = Path(args.output_workdir)
    variants = selected_scan_variants(
        args.variant, extended_strategy_scan=args.extended_strategy_scan
    )
    iterations = tuple(sorted(set(args.iteration or ())))
    iteration_variants = tuple(
        args.iteration_variant or DEFAULT_ITERATION_SCAN_VARIANTS
    )

    if args.run_strategy_caps and args.no_run:
        raise ValueError("--run-strategy-caps requires running Wannier90")

    if output_workdir.exists():
        if not args.force:
            raise ValueError(f"{output_workdir}: output directory already exists")
        shutil.rmtree(output_workdir)
    output_workdir.mkdir(parents=True)

    if args.run_source:
        run_wannier90(args.executable, args.seed, source_workdir, args.timeout)

    run_specs = [f"source:{args.seed}:{source_workdir}"]
    cap_runs: list[StrategyCapRun] = []
    for variant in variants:
        variant_workdir = output_workdir / variant
        prepare_variant(args.seed, source_workdir, variant_workdir, variant, force=True)
        if not args.no_run:
            try:
                run_wannier90(args.executable, args.seed, variant_workdir, args.timeout)
            except (OSError, RuntimeError, subprocess.TimeoutExpired):
                if not args.allow_failed:
                    raise
        run_specs.append(f"{variant}:{args.seed}:{variant_workdir}")
    for variant in iteration_variants:
        if variant.endswith("-zero-iter"):
            raise ValueError(
                "iteration variants must be full variants, not '*-zero-iter' variants"
            )
        for iteration in iterations:
            if iteration < 0:
                raise ValueError("iteration budgets must be non-negative")
            label = f"{variant}-iter-{iteration}"
            variant_workdir = output_workdir / label
            prepare_variant(
                args.seed, source_workdir, variant_workdir, variant, force=True
            )
            set_win_num_iter(args.seed, variant_workdir, iteration)
            if not args.no_run:
                try:
                    run_wannier90(
                        args.executable, args.seed, variant_workdir, args.timeout
                    )
                except (OSError, RuntimeError, subprocess.TimeoutExpired):
                    if not args.allow_failed:
                        raise
            run_specs.append(f"{label}:{args.seed}:{variant_workdir}")

    if args.run_strategy_caps:
        base_specs = [parse_run_spec(spec) for spec in run_specs]
        for record in select_strategy_cap_records(
            base_specs,
            allow_failed=args.allow_failed,
            drift_tol=args.history_drift_tol,
            uphill_tol=args.history_uphill_tol,
            all_late_drift=args.strategy_cap_all,
        ):
            cap_iteration = record.history.best.iteration
            cap_label = f"{record.spec.label}-cap-{cap_iteration}"
            cap_workdir = output_workdir / cap_label
            prepare_strategy_cap_run(
                args.seed, source_workdir, cap_workdir, record.spec.label, cap_iteration
            )
            try:
                run_wannier90(args.executable, args.seed, cap_workdir, args.timeout)
            except (OSError, RuntimeError, subprocess.TimeoutExpired):
                if not args.allow_failed:
                    raise
            cap_spec = RunSpec(cap_label, args.seed, cap_workdir)
            cap_runs.append(StrategyCapRun(record, cap_spec))
            run_specs.append(f"{cap_label}:{args.seed}:{cap_workdir}")

    specs = [parse_run_spec(spec) for spec in run_specs]
    matrix_args = argparse.Namespace(
        runs=run_specs,
        reference=args.reference,
        identity_tol=args.identity_tol,
        allow_failed=args.allow_failed,
        match_centers=args.match_centers,
        periodic_centers=args.periodic_centers,
        compare_export=args.compare_export,
        compare_amn=args.compare_amn,
    )
    status = command_matrix(matrix_args)
    if not args.no_relaxation_summary:
        periodic_cell = None
        if args.periodic_centers:
            reference_spec = next(
                (spec for spec in specs if spec.label == args.reference), None
            )
            if reference_spec is None:
                raise ValueError(f"unknown reference label: {args.reference}")
            periodic_cell = parse_win_cell(reference_spec.seed, reference_spec.workdir)
        print_relaxation_summary(
            specs,
            args.reference,
            allow_failed=args.allow_failed,
            match_centers=bool(args.match_centers),
            periodic_cell=periodic_cell,
        )
    if not args.no_iteration_summary:
        periodic_cell = None
        if args.periodic_centers:
            reference_spec = next(
                (spec for spec in specs if spec.label == args.reference), None
            )
            if reference_spec is None:
                raise ValueError(f"unknown reference label: {args.reference}")
            periodic_cell = parse_win_cell(reference_spec.seed, reference_spec.workdir)
        print_iteration_summary(
            specs,
            args.reference,
            allow_failed=args.allow_failed,
            match_centers=bool(args.match_centers),
            periodic_cell=periodic_cell,
        )
    if not args.no_history_summary:
        status = max(
            status,
            print_history_summary(
                specs,
                allow_failed=args.allow_failed,
                drift_tol=args.history_drift_tol,
                uphill_tol=args.history_uphill_tol,
                fail_on_drift=args.fail_on_history_drift,
            ),
        )
    if not args.no_strategy_summary:
        print_strategy_summary(
            specs,
            allow_failed=args.allow_failed,
            drift_tol=args.history_drift_tol,
            uphill_tol=args.history_uphill_tol,
        )
    if cap_runs:
        periodic_cell = None
        if args.periodic_centers:
            reference_spec = next(
                (spec for spec in specs if spec.label == args.reference), None
            )
            if reference_spec is None:
                raise ValueError(f"unknown reference label: {args.reference}")
            periodic_cell = parse_win_cell(reference_spec.seed, reference_spec.workdir)
        print_strategy_cap_summary(
            cap_runs,
            allow_failed=args.allow_failed,
            drift_tol=args.history_drift_tol,
            uphill_tol=args.history_uphill_tol,
            match_centers=bool(args.match_centers),
            periodic_cell=periodic_cell,
        )
    return status


def command_matrix(args: argparse.Namespace) -> int:
    specs = [parse_run_spec(spec) for spec in args.runs]
    reference_label = args.reference or specs[0].label
    reference_spec = next(
        (spec for spec in specs if spec.label == reference_label), None
    )
    if reference_spec is None:
        raise ValueError(f"unknown reference label: {reference_label}")
    reference = parse_wout(reference_spec.seed, reference_spec.workdir)
    periodic_cell = None
    if args.periodic_centers:
        periodic_cell = parse_win_cell(reference_spec.seed, reference_spec.workdir)

    header = [
        "label",
        "seed",
        "status",
        "num_wann",
        "omega_total",
        "omega_i",
        "omega_d",
        "omega_od",
        "all_done",
        "amn_identity",
        "amn_max_identity_deviation",
        "amn_max_offdiag",
        "amn_nonzero_fraction",
        "omega_total_delta",
        "spread_delta",
        "center_delta",
        "error",
    ]
    if args.compare_export:
        header.extend(
            [
                "export_status",
                "eig_max_abs_delta",
                "mmn_max_abs_delta",
                "mmn_rms_abs_delta",
                "mmn_max_singular_delta",
                "mmn_rms_singular_delta",
            ]
        )
    if args.compare_amn:
        header.extend(
            [
                "amn_compare_status",
                "amn_max_abs_delta",
                "amn_rms_abs_delta",
                "amn_max_singular_delta",
                "amn_max_subspace_deviation",
                "amn_min_subspace_svalue",
            ]
        )
    print("\t".join(header))
    for spec in specs:
        try:
            summary = parse_wout(spec.seed, spec.workdir)
        except (OSError, ValueError) as exc:
            if not args.allow_failed:
                raise
            row = [
                spec.label,
                spec.seed,
                "failed",
                "na",
                "nan",
                "nan",
                "nan",
                "nan",
                "no",
                "missing",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                "nan",
                str(exc),
            ]
            if args.compare_export:
                row.extend(["not-run", "nan", "nan", "nan", "nan", "nan"])
            if args.compare_amn:
                row.extend(["not-run", "nan", "nan", "nan", "nan", "nan"])
            print("\t".join(row))
            continue
        amn_identity = "missing"
        amn_max_identity_deviation = "nan"
        amn_max_offdiag = "nan"
        amn_nonzero_fraction = "nan"
        if (spec.workdir / f"{spec.seed}.amn").exists():
            amn_summary = parse_amn(spec.seed, spec.workdir, args.identity_tol)
            amn_identity = "yes" if amn_summary.is_identity else "no"
            amn_max_identity_deviation = f"{amn_summary.max_identity_deviation:.12e}"
            amn_max_offdiag = f"{amn_summary.max_offdiag:.12e}"
            amn_nonzero_fraction = f"{amn_summary.nonzero_fraction:.12e}"
        comparison = compare_summaries(
            reference,
            summary,
            match_centers=bool(args.match_centers),
            periodic_cell=periodic_cell,
        )
        row = [
            spec.label,
            spec.seed,
            "ok",
            str(len(summary.functions)),
            f"{summary.omega_total:.12e}",
            f"{summary.omega_i:.12e}",
            f"{summary.omega_d:.12e}",
            f"{summary.omega_od:.12e}",
            "yes" if summary.all_done else "no",
            amn_identity,
            amn_max_identity_deviation,
            amn_max_offdiag,
            amn_nonzero_fraction,
            f"{comparison.omega_total_delta:.12e}",
            f"{comparison.spread_delta:.12e}",
            f"{comparison.center_delta:.12e}",
            "",
        ]
        if args.compare_export:
            try:
                export_comparison = compare_exports(
                    reference_spec.seed, reference_spec.workdir, spec.seed, spec.workdir
                )
                row.extend(
                    [
                        "ok",
                        f"{export_comparison.eig_max_abs_delta:.12e}",
                        f"{export_comparison.mmn_max_abs_delta:.12e}",
                        f"{export_comparison.mmn_rms_abs_delta:.12e}",
                        f"{export_comparison.mmn_max_singular_delta:.12e}",
                        f"{export_comparison.mmn_rms_singular_delta:.12e}",
                    ]
                )
            except (OSError, ValueError) as exc:
                if not args.allow_failed:
                    raise
                row.extend(["failed", "nan", "nan", "nan", "nan", "nan"])
                row[-6] = f"failed:{exc}"
        if args.compare_amn:
            try:
                amn_comparison = compare_amn_files(
                    reference_spec.seed, reference_spec.workdir, spec.seed, spec.workdir
                )
                row.extend(
                    [
                        "ok",
                        f"{amn_comparison.max_abs_delta:.12e}",
                        f"{amn_comparison.rms_abs_delta:.12e}",
                        f"{amn_comparison.max_singular_delta:.12e}",
                        f"{amn_comparison.max_subspace_deviation:.12e}",
                        f"{amn_comparison.min_subspace_svalue:.12e}",
                    ]
                )
            except (OSError, ValueError) as exc:
                if not args.allow_failed:
                    raise
                row.extend(["failed", "nan", "nan", "nan", "nan", "nan"])
                row[-6] = f"failed:{exc}"
        print("\t".join(row))
    return 0


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
    compare_parser.add_argument(
        "--match-centers",
        "--match-centres",
        action="store_true",
        help="compare Wannier centres after optimal one-to-one matching",
    )
    compare_parser.add_argument(
        "--periodic-centers",
        "--periodic-centres",
        action="store_true",
        help="measure centre distances modulo the left seed unit cell",
    )
    compare_parser.set_defaults(func=command_compare)

    export_parser = subparsers.add_parser(
        "compare-export", help="compare CP2K-generated .eig and .mmn files"
    )
    export_parser.add_argument("left_seed")
    export_parser.add_argument("right_seed")
    export_parser.add_argument("--left-workdir", default=".")
    export_parser.add_argument("--right-workdir", default=".")
    export_parser.add_argument("--eig-tol", type=float, default=None)
    export_parser.add_argument("--mmn-tol", type=float, default=None)
    export_parser.add_argument("--mmn-singular-tol", type=float, default=None)
    export_parser.set_defaults(func=command_compare_export)

    compare_amn_parser = subparsers.add_parser(
        "compare-amn", help="compare Wannier90 .amn projection matrices"
    )
    compare_amn_parser.add_argument("left_seed")
    compare_amn_parser.add_argument("right_seed")
    compare_amn_parser.add_argument("--left-workdir", default=".")
    compare_amn_parser.add_argument("--right-workdir", default=".")
    compare_amn_parser.add_argument("--amn-tol", type=float, default=None)
    compare_amn_parser.add_argument("--amn-singular-tol", type=float, default=None)
    compare_amn_parser.add_argument("--amn-subspace-tol", type=float, default=None)
    compare_amn_parser.set_defaults(func=command_compare_amn)

    amn_parser = subparsers.add_parser(
        "inspect-amn", help="inspect a Wannier90 projection matrix"
    )
    amn_parser.add_argument("seed", help="Wannier90 seed name")
    amn_parser.add_argument(
        "--workdir", default=".", help="directory containing the seed files"
    )
    amn_parser.add_argument("--identity-tol", type=float, default=1.0e-12)
    amn_parser.add_argument("--max-identity-deviation", type=float, default=None)
    amn_parser.add_argument("--require-identity", action="store_true")
    amn_parser.set_defaults(func=command_inspect_amn)

    variant_parser = subparsers.add_parser(
        "prepare-variant", help="prepare a Wannier90 research-input variant"
    )
    variant_parser.add_argument("seed", help="Wannier90 seed name")
    variant_parser.add_argument("variant", choices=VARIANT_CHOICES)
    variant_parser.add_argument("--source-workdir", default=".")
    variant_parser.add_argument("--output-workdir", required=True)
    variant_parser.add_argument("--force", action="store_true")
    variant_parser.set_defaults(func=command_prepare_variant)

    scan_parser = subparsers.add_parser(
        "scan-variants",
        help="prepare, run, and summarize Wannier90 projection variants",
    )
    scan_parser.add_argument("seed", help="Wannier90 seed name")
    scan_parser.add_argument("--source-workdir", default=".")
    scan_parser.add_argument("--output-workdir", required=True)
    scan_parser.add_argument(
        "--variant",
        action="append",
        choices=VARIANT_CHOICES,
        help=(
            "variant to scan; may be repeated. Defaults to transport/SCDM variants, "
            "or the broader projection set with --extended-strategy-scan."
        ),
    )
    scan_parser.add_argument(
        "--extended-strategy-scan",
        action="store_true",
        help="scan the broader set of AMN projection variants",
    )
    scan_parser.add_argument(
        "--reference",
        default="source",
        help="matrix reference label, usually source or one of the selected variants",
    )
    scan_parser.add_argument("--identity-tol", type=float, default=1.0e-12)
    scan_parser.add_argument("--executable", default="wannier90.x")
    scan_parser.add_argument("--timeout", type=float, default=300.0)
    scan_parser.add_argument(
        "--iteration",
        action="append",
        type=int,
        help="additional Wannier90 num_iter budget to run for iteration-path diagnostics",
    )
    scan_parser.add_argument(
        "--iteration-variant",
        action="append",
        choices=VARIANT_CHOICES,
        help="full variant to use for --iteration scans; defaults to SCDM/transport variants",
    )
    scan_parser.add_argument(
        "--run-strategy-caps",
        action="store_true",
        help="run capped num_iter checks for late-drift projection strategies",
    )
    scan_parser.add_argument(
        "--strategy-cap-all",
        action="store_true",
        help="with --run-strategy-caps, cap all late-drift strategies instead of the best one",
    )
    scan_parser.add_argument("--force", action="store_true")
    scan_parser.add_argument("--allow-failed", action="store_true")
    scan_parser.add_argument(
        "--no-run",
        action="store_true",
        help="prepare variants but do not run wannier90.x before summarizing",
    )
    scan_parser.add_argument(
        "--run-source",
        action="store_true",
        help="run wannier90.x in the source directory before scanning variants",
    )
    scan_parser.add_argument(
        "--no-relaxation-summary",
        action="store_true",
        help="suppress the zero-iteration to relaxed-MLWF pairing table",
    )
    scan_parser.add_argument(
        "--no-iteration-summary",
        action="store_true",
        help="suppress the fixed-num_iter iteration-path table",
    )
    scan_parser.add_argument(
        "--no-history-summary",
        action="store_true",
        help="suppress the per-run Wannier90 iteration-history table",
    )
    scan_parser.add_argument(
        "--no-strategy-summary",
        action="store_true",
        help="suppress projection-strategy ranking from Wannier90 histories",
    )
    scan_parser.add_argument(
        "--history-drift-tol",
        type=float,
        default=1.0e-6,
        help="Omega_total tolerance for classifying late Wannier90 minimization drift",
    )
    scan_parser.add_argument(
        "--history-uphill-tol",
        type=float,
        default=1.0e-6,
        help="single-step Omega_total increase tolerance for history diagnostics",
    )
    scan_parser.add_argument(
        "--fail-on-history-drift",
        action="store_true",
        help="return a nonzero status if a run is classified as late_drift",
    )
    scan_parser.add_argument("--match-centers", "--match-centres", action="store_true")
    scan_parser.add_argument(
        "--periodic-centers", "--periodic-centres", action="store_true"
    )
    scan_parser.add_argument("--compare-export", action="store_true")
    scan_parser.add_argument("--compare-amn", action="store_true")
    scan_parser.set_defaults(func=command_scan_variants)

    strategy_matrix_parser = subparsers.add_parser(
        "strategy-matrix",
        help="run projection-strategy scans for several Wannier90 cases",
    )
    strategy_matrix_parser.add_argument(
        "cases", nargs="*", help="case specification as LABEL:SEED:SOURCE_WORKDIR"
    )
    strategy_matrix_parser.add_argument(
        "--case-file",
        action="append",
        help=(
            "file with one LABEL:SEED:SOURCE_WORKDIR case per line; blank lines "
            "and comments starting with # are ignored"
        ),
    )
    strategy_matrix_parser.add_argument("--output-workdir", required=True)
    strategy_matrix_parser.add_argument(
        "--variant",
        action="append",
        choices=VARIANT_CHOICES,
        help=(
            "variant to scan; may be repeated. Defaults to transport/SCDM variants, "
            "or the broader projection set with --extended-strategy-scan."
        ),
    )
    strategy_matrix_parser.add_argument(
        "--extended-strategy-scan",
        action="store_true",
        help="scan the broader set of AMN projection variants",
    )
    strategy_matrix_parser.add_argument(
        "--run-strategy-caps",
        action="store_true",
        help="run capped num_iter checks for late-drift projection strategies",
    )
    strategy_matrix_parser.add_argument(
        "--strategy-cap-all",
        action="store_true",
        help="with --run-strategy-caps, cap all late-drift strategies instead of the best one",
    )
    strategy_matrix_parser.add_argument("--executable", default="wannier90.x")
    strategy_matrix_parser.add_argument("--timeout", type=float, default=300.0)
    strategy_matrix_parser.add_argument("--force", action="store_true")
    strategy_matrix_parser.add_argument("--allow-failed", action="store_true")
    strategy_matrix_parser.add_argument(
        "--run-source",
        action="store_true",
        help="run wannier90.x in each source directory before scanning variants",
    )
    strategy_matrix_parser.add_argument(
        "--history-drift-tol",
        type=float,
        default=1.0e-6,
        help="Omega_total tolerance for classifying late Wannier90 minimization drift",
    )
    strategy_matrix_parser.add_argument(
        "--history-uphill-tol",
        type=float,
        default=1.0e-6,
        help="single-step Omega_total increase tolerance for history diagnostics",
    )
    strategy_matrix_parser.add_argument(
        "--fail-on-history-drift",
        action="store_true",
        help="return a nonzero status if a selected run is classified as late_drift",
    )
    strategy_matrix_parser.add_argument(
        "--fail-on-uncapped-drift",
        action="store_true",
        help=(
            "return a nonzero status when a late_drift strategy has no verified "
            "capped-num_iter reproduction"
        ),
    )
    strategy_matrix_parser.add_argument(
        "--robust-minimization",
        action="store_true",
        help=(
            "automatically run and verify capped num_iter reproductions for all "
            "late-drift strategies"
        ),
    )
    strategy_matrix_parser.add_argument(
        "--cap-tol",
        type=float,
        default=1.0e-6,
        help="absolute tolerance for verifying capped late-drift reproductions",
    )
    strategy_matrix_parser.add_argument(
        "--no-summary",
        action="store_true",
        help="suppress the aggregate strategy_matrix_summary table",
    )
    strategy_matrix_parser.add_argument(
        "--json-output", help="write strategy-matrix rows and summary to this JSON file"
    )
    strategy_matrix_parser.add_argument(
        "--markdown-output",
        help="write a compact Markdown strategy report to this file",
    )
    strategy_matrix_parser.add_argument(
        "--match-centers", "--match-centres", action="store_true"
    )
    strategy_matrix_parser.add_argument(
        "--periodic-centers", "--periodic-centres", action="store_true"
    )
    strategy_matrix_parser.set_defaults(func=command_strategy_matrix)

    cp2k_strategy_matrix_parser = subparsers.add_parser(
        "cp2k-strategy-matrix",
        help="run CP2K Wannier90 exports before a projection-strategy matrix",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "cases", nargs="*", help="CP2K case specification as LABEL:INPUT_PATH[:SEED]"
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--case-file",
        action="append",
        help=(
            "file with one LABEL:INPUT_PATH[:SEED] CP2K case per line; blank "
            "lines and comments starting with # are ignored"
        ),
    )
    cp2k_strategy_matrix_parser.add_argument("--output-workdir", required=True)
    cp2k_strategy_matrix_parser.add_argument(
        "--cp2k-executable",
        default="cp2k.psmp",
        help="CP2K executable used to generate Wannier90 export files",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--cp2k-timeout",
        type=float,
        default=600.0,
        help="CP2K execution timeout per case in seconds",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--copy-input-dir-files",
        action="store_true",
        help="copy sibling files from each CP2K input directory into the run directory",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--no-cp2k-run",
        action="store_true",
        help="reuse existing cp2k_exports directories below --output-workdir",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--cp2k-initial-projection",
        action="append",
        choices=CP2K_INITIAL_PROJECTIONS,
        help=(
            "expand each CP2K case into a WANNIER90 export with the selected "
            "INITIAL_PROJECTIONS mode; IDENTITY also enables USE_BLOCH_PHASES"
        ),
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--cp2k-projection-smoothing",
        action="append",
        choices=CP2K_PROJECTION_SMOOTHINGS,
        help=(
            "expand each CP2K case into a WANNIER90 export with the selected "
            "PROJECTION_SMOOTHING mode; may be repeated. When used without "
            "--cp2k-initial-projection, the input projection mode is preserved."
        ),
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--cp2k-auto-projection-score",
        action="append",
        choices=CP2K_AUTO_PROJECTION_SCORES,
        help=(
            "expand AUTO_SCDM_RAW CP2K cases over the selected AUTO_PROJECTION_SCORE "
            "modes; requires --cp2k-initial-projection AUTO_SCDM_RAW"
        ),
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--variant",
        action="append",
        choices=VARIANT_CHOICES,
        help=(
            "variant to scan; may be repeated. Defaults to transport/SCDM variants, "
            "or the broader projection set with --extended-strategy-scan."
        ),
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--extended-strategy-scan",
        action="store_true",
        help="scan the broader set of AMN projection variants",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--run-strategy-caps",
        action="store_true",
        help="run capped num_iter checks for late-drift projection strategies",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--strategy-cap-all",
        action="store_true",
        help="with --run-strategy-caps, cap all late-drift strategies instead of the best one",
    )
    cp2k_strategy_matrix_parser.add_argument("--executable", default="wannier90.x")
    cp2k_strategy_matrix_parser.add_argument("--timeout", type=float, default=300.0)
    cp2k_strategy_matrix_parser.add_argument("--force", action="store_true")
    cp2k_strategy_matrix_parser.add_argument("--allow-failed", action="store_true")
    cp2k_strategy_matrix_parser.add_argument(
        "--no-run-source",
        action="store_true",
        help="do not run wannier90.x in CP2K export directories before scanning variants",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--history-drift-tol",
        type=float,
        default=1.0e-6,
        help="Omega_total tolerance for classifying late Wannier90 minimization drift",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--history-uphill-tol",
        type=float,
        default=1.0e-6,
        help="single-step Omega_total increase tolerance for history diagnostics",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--fail-on-history-drift",
        action="store_true",
        help="return a nonzero status if a selected run is classified as late_drift",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--fail-on-uncapped-drift",
        action="store_true",
        help=(
            "return a nonzero status when a late_drift strategy has no verified "
            "capped-num_iter reproduction"
        ),
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--robust-minimization",
        action="store_true",
        help=(
            "automatically run and verify capped num_iter reproductions for all "
            "late-drift strategies"
        ),
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--cap-tol",
        type=float,
        default=1.0e-6,
        help="absolute tolerance for verifying capped late-drift reproductions",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--no-summary",
        action="store_true",
        help="suppress the aggregate strategy_matrix_summary table",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--json-output", help="write strategy-matrix rows and summary to this JSON file"
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--markdown-output",
        help="write a compact Markdown strategy report to this file",
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--match-centers", "--match-centres", action="store_true"
    )
    cp2k_strategy_matrix_parser.add_argument(
        "--periodic-centers", "--periodic-centres", action="store_true"
    )
    cp2k_strategy_matrix_parser.set_defaults(func=command_cp2k_strategy_matrix)

    cp2k_projection_matrix_parser = subparsers.add_parser(
        "cp2k-projection-matrix",
        help="run CP2K Wannier90 projection candidates and summarize their MLWF spreads",
    )
    cp2k_projection_matrix_parser.add_argument(
        "cases", nargs="*", help="CP2K case specification as LABEL:INPUT_PATH[:SEED]"
    )
    cp2k_projection_matrix_parser.add_argument(
        "--case-file",
        action="append",
        help=(
            "file with one LABEL:INPUT_PATH[:SEED] CP2K case per line; blank "
            "lines and comments starting with # are ignored"
        ),
    )
    cp2k_projection_matrix_parser.add_argument("--output-workdir", required=True)
    cp2k_projection_matrix_parser.add_argument(
        "--cp2k-executable",
        default="cp2k.psmp",
        help="CP2K executable used to generate Wannier90 export files",
    )
    cp2k_projection_matrix_parser.add_argument(
        "--cp2k-timeout",
        type=float,
        default=600.0,
        help="CP2K execution timeout per case in seconds",
    )
    cp2k_projection_matrix_parser.add_argument(
        "--copy-input-dir-files",
        action="store_true",
        help="copy sibling files from each CP2K input directory into the run directory",
    )
    cp2k_projection_matrix_parser.add_argument(
        "--no-cp2k-run",
        action="store_true",
        help="summarize an existing projection-matrix output directory",
    )
    cp2k_projection_matrix_parser.add_argument(
        "--no-wannier90-run",
        action="store_true",
        help="do not run wannier90.x before parsing .wout files",
    )
    cp2k_projection_matrix_parser.add_argument(
        "--cp2k-initial-projection",
        action="append",
        choices=CP2K_INITIAL_PROJECTIONS,
        help=(
            "projection mode to export; defaults to AO_SCDM_RAW, AO_PAIR_SCDM_RAW, "
            "AO_PAIR_SIGNED_SCDM_RAW, CENTER_SCDM_RAW, and AUTO_SCDM_RAW"
        ),
    )
    cp2k_projection_matrix_parser.add_argument(
        "--cp2k-projection-smoothing",
        action="append",
        choices=CP2K_PROJECTION_SMOOTHINGS,
        help="projection smoothing mode to use; defaults to MMN_TRANSPORT",
    )
    cp2k_projection_matrix_parser.add_argument(
        "--cp2k-auto-projection-score",
        action="append",
        choices=CP2K_AUTO_PROJECTION_SCORES,
        help=(
            "AUTO_SCDM_RAW score to export; defaults to all AUTO_PROJECTION_SCORE modes"
        ),
    )
    cp2k_projection_matrix_parser.add_argument(
        "--early-minimization-iteration",
        action="append",
        type=int,
        help=(
            "copy each direct CP2K projection export, rerun Wannier90 with this "
            "num_iter budget, and rank projections by the best early Omega_total"
        ),
    )
    cp2k_projection_matrix_parser.add_argument(
        "--history-drift-tol",
        type=float,
        default=1.0e-6,
        help="Omega_total tolerance for classifying late Wannier90 minimization drift",
    )
    cp2k_projection_matrix_parser.add_argument(
        "--history-uphill-tol",
        type=float,
        default=1.0e-6,
        help="single-step Omega_total increase tolerance for history diagnostics",
    )
    cp2k_projection_matrix_parser.add_argument("--executable", default="wannier90.x")
    cp2k_projection_matrix_parser.add_argument("--timeout", type=float, default=300.0)
    cp2k_projection_matrix_parser.add_argument("--force", action="store_true")
    cp2k_projection_matrix_parser.add_argument("--allow-failed", action="store_true")
    cp2k_projection_matrix_parser.add_argument(
        "--no-summary", action="store_true", help="do not print the summary table"
    )
    cp2k_projection_matrix_parser.add_argument(
        "--markdown-output", help="write the Markdown summary to this file too"
    )
    cp2k_projection_matrix_parser.set_defaults(func=command_cp2k_projection_matrix)

    strategy_report_parser = subparsers.add_parser(
        "strategy-report",
        help="render a Markdown report from a strategy-matrix JSON file",
    )
    strategy_report_parser.add_argument("json_input")
    strategy_report_parser.add_argument("--output", help="write report to this file")
    strategy_report_parser.add_argument(
        "--cap-tol",
        type=float,
        default=1.0e-6,
        help="absolute tolerance for verifying capped late-drift reproductions",
    )
    strategy_report_parser.set_defaults(func=command_strategy_report)

    matrix_parser = subparsers.add_parser(
        "matrix", help="summarize and compare several Wannier90 .wout runs"
    )
    matrix_parser.add_argument(
        "runs", nargs="+", help="run specification as LABEL:SEED:WORKDIR"
    )
    matrix_parser.add_argument("--reference", default=None)
    matrix_parser.add_argument("--identity-tol", type=float, default=1.0e-12)
    matrix_parser.add_argument("--allow-failed", action="store_true")
    matrix_parser.add_argument(
        "--match-centers", "--match-centres", action="store_true"
    )
    matrix_parser.add_argument(
        "--periodic-centers", "--periodic-centres", action="store_true"
    )
    matrix_parser.add_argument(
        "--compare-export",
        action="store_true",
        help="include .eig/.mmn deltas against the reference run",
    )
    matrix_parser.add_argument(
        "--compare-amn",
        action="store_true",
        help="include .amn matrix and subspace deltas against the reference run",
    )
    matrix_parser.set_defaults(func=command_matrix)
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
