#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
#   CP2K: A general program to perform molecular dynamics simulations
#   Copyright 2000-2026 CP2K developers group <https://cp2k.org>
#
#   SPDX-License-Identifier: GPL-2.0-or-later
# --------------------------------------------------------------------------------------------------

"""Plot and export CP2K band structure (.bs) files.

Supports standard Kohn-Sham band structures and unfolded/projected band
structures produced by &DFT%PRINT%BAND_STRUCTURE (including &UNFOLDING,
spin-polarized calculations, multiple KPOINT_SET blocks, and ATOMS groups).
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

# Match src/common/physcon.F (CODATA constants used by CP2K).
ANGSTROM_PER_BOHR: float = 0.52917720859
EV_PER_HARTREE: float = (
    2.0 * 10973731.568527 * 6.62606896e-34 * 299792458.0 / 1.602176487e-19
)
COORD_MATCH_TOL: float = 1.0e-5
DEFAULT_ENERGY_GRID: int = 600
DEFAULT_DPI: int = 150
DEFAULT_SIGMA_EV: float = 0.03
DEFAULT_DOT_STRIDE: int = 4
SPECTRUM_CUTOFF_SIGMA: float = 4.5
SPECTRUM_PERCENTILE_MAX: float = 99.5

_RE_PRIM_VEC = re.compile(
    r"^#\s*UNFOLD_BAND\|\s*primitive cell vector\s+(\d+)\s+\[Bohr\]\s*=\s*"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"
)
_RE_ATOM_GROUPS = re.compile(r"^#\s*UNFOLD_BAND\|\s*atom groups:\s*(\d+)")
_RE_SET_HEADER = re.compile(
    r"^#\s*Set\s+(\d+):\s*(\d+)\s+special points,\s*(\d+)\s+k-points,\s*(\d+)\s+bands"
)
_RE_SPECIAL_POINT = re.compile(
    r"^#\s*Special point\s+(\d+)\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)(?:\s+(.*))?$"
)
_RE_POINT_HEADER = re.compile(
    r"^#\s*Point\s+(\d+)\s+Spin\s+(\d+):\s*"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"
)
_RE_FERMI = re.compile(r"Fermi energy:\s+([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)")
_RE_CELL_VEC = re.compile(
    r"CELL\|\s*Vector\s+([abc])\s+\[angstrom\]:\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+"
    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"
)


@dataclass(frozen=True)
class SpecialPoint:
    """High-symmetry special point in a KPOINT_SET."""

    coords: tuple[float, float, float]
    label: str


@dataclass(frozen=True)
class KPointData:
    """Eigenvalues, occupations, and optional unfolding weights at one k-point."""

    index: int
    coords: tuple[float, float, float]
    k_weight: float
    energies: FloatArray
    occupations: FloatArray
    unfold_weight: Optional[FloatArray]
    group_weights: Optional[FloatArray]


@dataclass(frozen=True)
class BandSet:
    """One KPOINT_SET block in a CP2K .bs file."""

    set_index: int
    special_points: list[SpecialPoint]
    points_by_spin: dict[int, list[KPointData]]

    @property
    def spins(self) -> list[int]:
        return sorted(self.points_by_spin.keys())


@dataclass(frozen=True)
class BandStructure:
    """Parsed contents of a CP2K .bs file."""

    sets: list[BandSet]
    is_unfolded: bool
    n_groups: int
    prim_cell_angstrom: Optional[FloatArray]

    @property
    def spins(self) -> list[int]:
        found: set[int] = set()
        for bset in self.sets:
            found.update(bset.points_by_spin.keys())
        return sorted(found)


@dataclass(frozen=True)
class Cp2kOutInfo:
    """Metadata extracted from a CP2K output (.out) file."""

    fermi_ev: Optional[float]
    cell_angstrom: Optional[FloatArray]


@dataclass(frozen=True)
class KPathGeometry:
    """Cumulative k-path coordinates and high-symmetry tick locations."""

    set_x: list[FloatArray]
    tick_positions: list[float]
    tick_labels: list[str]
    break_positions: list[float]


def format_kpoint_label(raw: str) -> str:
    """Format a CP2K SPECIAL_POINT label for display."""
    cleaned = raw.strip()
    if not cleaned or cleaned.lower() == "not specified":
        return ""
    if cleaned.upper() in {"GAMMA", "G"}:
        return "\u0393"
    return cleaned


def _coords_close(a: tuple[float, float, float], b: tuple[float, float, float]) -> bool:
    return (
        abs(a[0] - b[0]) < COORD_MATCH_TOL
        and abs(a[1] - b[1]) < COORD_MATCH_TOL
        and abs(a[2] - b[2]) < COORD_MATCH_TOL
    )


def parse_cp2k_out(path: Path) -> Cp2kOutInfo:
    """Extract the last Fermi energy [eV] and lattice vectors [Angstrom] from CP2K output."""
    fermi_ha: Optional[float] = None
    cell_rows: dict[str, tuple[float, float, float]] = {}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            m_ef = _RE_FERMI.search(line)
            if m_ef is not None:
                fermi_ha = float(m_ef.group(1))
                continue
            m_cell = _RE_CELL_VEC.search(line)
            if m_cell is not None:
                cell_rows[m_cell.group(1)] = (
                    float(m_cell.group(2)),
                    float(m_cell.group(3)),
                    float(m_cell.group(4)),
                )

    cell_mat: Optional[FloatArray] = None
    if {"a", "b", "c"}.issubset(cell_rows):
        cell_mat = np.array(
            [cell_rows["a"], cell_rows["b"], cell_rows["c"]], dtype=np.float64
        )
    fermi_ev = fermi_ha * EV_PER_HARTREE if fermi_ha is not None else None
    return Cp2kOutInfo(fermi_ev=fermi_ev, cell_angstrom=cell_mat)


def parse_bs_file(path: Path) -> BandStructure:
    """Parse a CP2K .bs file (standard or unfolded).

    CP2K appends to existing .bs files across reruns; if multiple Set 1 blocks
    are detected, only the last run is kept and a warning is written to stderr.
    """
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    # Detect appended reruns (repeated "# Set 1:").
    set1_lines = [
        i
        for i, line in enumerate(lines)
        if (m := _RE_SET_HEADER.match(line.strip())) is not None
        and int(m.group(1)) == 1
    ]
    if len(set1_lines) > 1:
        print(
            f"WARNING: multiple appended runs detected in '{path}'; using the last run.",
            file=sys.stderr,
        )
        start = set1_lines[-1]
        while start > 0 and lines[start - 1].lstrip().startswith("# UNFOLD_BAND|"):
            start -= 1
        lines = lines[start:]

    prim_rows: dict[int, tuple[float, float, float]] = {}
    n_groups = 0
    is_unfolded = False
    sets: list[BandSet] = []

    cur_set_idx: Optional[int] = None
    cur_specials: list[SpecialPoint] = []
    cur_points: dict[int, list[KPointData]] = {}

    def _flush_set() -> None:
        if cur_set_idx is not None:
            sets.append(
                BandSet(
                    set_index=cur_set_idx,
                    special_points=list(cur_specials),
                    points_by_spin={k: list(v) for k, v in cur_points.items()},
                )
            )

    i = 0
    n_lines = len(lines)
    while i < n_lines:
        stripped = lines[i].strip()
        if not stripped:
            i += 1
            continue

        m_pv = _RE_PRIM_VEC.match(stripped)
        if m_pv is not None:
            is_unfolded = True
            prim_rows[int(m_pv.group(1))] = (
                float(m_pv.group(2)) * ANGSTROM_PER_BOHR,
                float(m_pv.group(3)) * ANGSTROM_PER_BOHR,
                float(m_pv.group(4)) * ANGSTROM_PER_BOHR,
            )
            i += 1
            continue

        m_ag = _RE_ATOM_GROUPS.match(stripped)
        if m_ag is not None:
            is_unfolded = True
            n_groups = int(m_ag.group(1))
            i += 1
            continue

        m_set = _RE_SET_HEADER.match(stripped)
        if m_set is not None:
            _flush_set()
            cur_set_idx = int(m_set.group(1))
            cur_specials = []
            cur_points = {}
            i += 1
            continue

        m_sp = _RE_SPECIAL_POINT.match(stripped)
        if m_sp is not None:
            cur_specials.append(
                SpecialPoint(
                    coords=(
                        float(m_sp.group(2)),
                        float(m_sp.group(3)),
                        float(m_sp.group(4)),
                    ),
                    label=format_kpoint_label(m_sp.group(5) or ""),
                )
            )
            i += 1
            continue

        m_pt = _RE_POINT_HEADER.match(stripped)
        if m_pt is not None:
            pt_idx = int(m_pt.group(1))
            spin = int(m_pt.group(2))
            coords = (
                float(m_pt.group(3)),
                float(m_pt.group(4)),
                float(m_pt.group(5)),
            )
            k_weight = float(m_pt.group(6))
            i += 1
            if i < n_lines and lines[i].lstrip().startswith("#"):
                i += 1

            rows: list[list[float]] = []
            while i < n_lines:
                row_s = lines[i].strip()
                if not row_s or row_s.startswith("#"):
                    break
                rows.append([float(tok) for tok in row_s.split()])
                i += 1

            arr = np.asarray(rows, dtype=np.float64)
            if arr.ndim != 2 or arr.shape[1] < 3:
                raise ValueError(
                    f"Malformed band table at point {pt_idx} spin {spin} in '{path}'."
                )

            energies = arr[:, 1]
            occupations = arr[:, 2]
            unfold_w: Optional[FloatArray] = None
            group_w: Optional[FloatArray] = None
            if arr.shape[1] >= 4:
                is_unfolded = True
                unfold_w = arr[:, 3]
                if arr.shape[1] > 4:
                    group_w = arr[:, 4:]
                    if n_groups == 0:
                        n_groups = group_w.shape[1]

            cur_points.setdefault(spin, []).append(
                KPointData(
                    index=pt_idx,
                    coords=coords,
                    k_weight=k_weight,
                    energies=energies,
                    occupations=occupations,
                    unfold_weight=unfold_w,
                    group_weights=group_w,
                )
            )
            continue

        i += 1

    _flush_set()
    if not sets:
        raise ValueError(f"No KPOINT_SET data found in '{path}'.")

    prim_cell: Optional[FloatArray] = None
    if {1, 2, 3}.issubset(prim_rows):
        prim_cell = np.array(
            [prim_rows[1], prim_rows[2], prim_rows[3]], dtype=np.float64
        )

    return BandStructure(
        sets=sets,
        is_unfolded=is_unfolded,
        n_groups=n_groups,
        prim_cell_angstrom=prim_cell,
    )


def compute_kpath_geometry(
    bs: BandStructure, cell_angstrom: Optional[FloatArray]
) -> KPathGeometry:
    """Compute cumulative k-path distance and high-symmetry ticks across sets."""
    lat = bs.prim_cell_angstrom if bs.prim_cell_angstrom is not None else cell_angstrom
    recip: Optional[FloatArray] = None
    if lat is not None and abs(float(np.linalg.det(lat))) > 1.0e-12:
        recip = 2.0 * math.pi * np.linalg.inv(lat).T

    set_x: list[FloatArray] = []
    tick_positions: list[float] = []
    tick_labels: list[str] = []
    break_positions: list[float] = []

    offset = 0.0
    prev_end_coords: Optional[tuple[float, float, float]] = None

    for set_i, bset in enumerate(bs.sets):
        ref_pts = bset.points_by_spin[bset.spins[0]]
        coords_arr = np.array([pt.coords for pt in ref_pts], dtype=np.float64)
        cart = coords_arr @ recip if recip is not None else coords_arr
        diffs = np.diff(cart, axis=0)
        seg_len = np.linalg.norm(diffs, axis=1)
        x_arr = np.concatenate(([0.0], np.cumsum(seg_len))) + offset
        set_x.append(x_arr)

        # Match each SPECIAL_POINT monotonically along this set.
        search_from = 0
        for sp_i, sp in enumerate(bset.special_points):
            matched_idx: Optional[int] = None
            for k_idx in range(search_from, len(ref_pts)):
                if _coords_close(ref_pts[k_idx].coords, sp.coords):
                    matched_idx = k_idx
                    search_from = k_idx + 1
                    break
            if matched_idx is None:
                continue

            pos = float(x_arr[matched_idx])
            label = (
                sp.label
                or f"({sp.coords[0]:.2g},{sp.coords[1]:.2g},{sp.coords[2]:.2g})"
            )

            if (
                set_i > 0
                and sp_i == 0
                and tick_positions
                and abs(pos - tick_positions[-1]) < 1.0e-9
            ):
                if prev_end_coords is not None and not _coords_close(
                    prev_end_coords, sp.coords
                ):
                    if tick_labels[-1] != label:
                        tick_labels[-1] = f"{tick_labels[-1]}|{label}"
                    break_positions.append(pos)
            else:
                tick_positions.append(pos)
                tick_labels.append(label)

        offset = float(x_arr[-1])
        prev_end_coords = ref_pts[-1].coords

    return KPathGeometry(
        set_x=set_x,
        tick_positions=tick_positions,
        tick_labels=tick_labels,
        break_positions=break_positions,
    )


def select_weights(pt: KPointData, group: Optional[int]) -> FloatArray:
    """Return total unfolding weight (group=None) or 1-based ATOMS group weight."""
    if group is None:
        if pt.unfold_weight is None:
            raise ValueError("Band structure file does not contain unfolding weights.")
        return pt.unfold_weight
    if pt.group_weights is None:
        raise ValueError(
            "Band structure file does not contain ATOMS group weights "
            "(requested --group requires multiple ATOMS groups in &UNFOLDING)."
        )
    if group < 1 or group > pt.group_weights.shape[1]:
        raise ValueError(
            f"Requested --group {group}, but file has {pt.group_weights.shape[1]} groups."
        )
    return pt.group_weights[:, group - 1]


def compute_spectral_map(
    points: Sequence[KPointData],
    energy_zero: float,
    e_grid: FloatArray,
    sigma_ev: float,
    group: Optional[int],
) -> FloatArray:
    """Evaluate Gaussian-broadened spectral function A(k, E) on e_grid."""
    n_e = len(e_grid)
    n_k = len(points)
    spec = np.zeros((n_e, n_k), dtype=np.float64)
    inv_2s2 = 0.5 / (sigma_ev * sigma_ev)
    cutoff = SPECTRUM_CUTOFF_SIGMA * sigma_ev
    e_lo = float(e_grid[0]) - cutoff
    e_hi = float(e_grid[-1]) + cutoff

    for ik, pt in enumerate(points):
        e_rel = pt.energies - energy_zero
        w = select_weights(pt, group)
        mask = (w > 0.0) & (e_rel >= e_lo) & (e_rel <= e_hi)
        if not np.any(mask):
            continue
        de = e_grid[:, None] - e_rel[mask][None, :]
        spec[:, ik] = np.exp(-de * de * inv_2s2) @ w[mask]
    return spec


def export_csv(
    bs: BandStructure,
    geom: KPathGeometry,
    energy_zero: float,
    out_path: Path,
) -> None:
    """Write a tidy CSV table of eigenvalues and optional unfolding/group weights."""
    header = [
        "set",
        "k_index",
        "k_dist",
        "kx",
        "ky",
        "kz",
        "spin",
        "band",
        "energy_ev",
        "occupation",
    ]
    if bs.is_unfolded:
        header.append("weight")
        for g in range(1, bs.n_groups + 1):
            header.append(f"w_grp{g}")

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for bset, x_arr in zip(bs.sets, geom.set_x):
            for sp in bset.spins:
                for ik, pt in enumerate(bset.points_by_spin[sp]):
                    x_val = float(x_arr[ik])
                    e_rel = pt.energies - energy_zero
                    for ib in range(len(e_rel)):
                        row: list[Any] = [
                            bset.set_index,
                            pt.index,
                            f"{x_val:.8f}",
                            f"{pt.coords[0]:.8f}",
                            f"{pt.coords[1]:.8f}",
                            f"{pt.coords[2]:.8f}",
                            sp,
                            ib + 1,
                            f"{float(e_rel[ib]):.8f}",
                            f"{float(pt.occupations[ib]):.6f}",
                        ]
                        if bs.is_unfolded and pt.unfold_weight is not None:
                            row.append(f"{float(pt.unfold_weight[ib]):.8f}")
                            if pt.group_weights is not None:
                                for ig in range(bs.n_groups):
                                    row.append(f"{float(pt.group_weights[ib, ig]):.8f}")
                        writer.writerow(row)


def _decorate_axis(
    ax: Any,
    geom: KPathGeometry,
    lo: float,
    hi: float,
    has_zero_ref: bool,
) -> None:
    """Apply high-symmetry vertical guides, Fermi line, ticks, and limits."""
    x_min = float(geom.set_x[0][0])
    x_max = float(geom.set_x[-1][-1])
    for xp in geom.tick_positions[1:-1]:
        ax.axvline(xp, color="0.55", linewidth=0.7, zorder=1)
    for bp in geom.break_positions:
        ax.axvline(bp, color="0.20", linewidth=1.2, zorder=2)
    if has_zero_ref:
        ax.axhline(0.0, color="0.35", linewidth=0.8, linestyle="--", zorder=2)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(lo, hi)
    ax.set_xticks(geom.tick_positions)
    ax.set_xticklabels(geom.tick_labels, fontsize=11)


def _resolve_spins(available: Sequence[int], spin_arg: str) -> list[int]:
    if spin_arg == "both":
        return list(available)
    sp = int(spin_arg)
    if sp not in available:
        raise ValueError(
            f"Requested --spin {sp}, but file has spins {list(available)}."
        )
    return [sp]


def _resolve_group(group_arg: str) -> Optional[int]:
    if group_arg.strip().lower() == "total":
        return None
    val = int(group_arg)
    if val < 1:
        raise ValueError("--group must be 'total' or a positive 1-based group index.")
    return val


def plot_bands(
    bs: BandStructure,
    geom: KPathGeometry,
    energy_zero: float,
    has_zero_ref: bool,
    out_path: Path,
    emin: Optional[float],
    emax: Optional[float],
    spin_arg: str,
    group_arg: str,
    mode: str,
    sigma_ev: float,
    dot_stride: int,
    cmap: str,
    vmax: Optional[float],
    title: Optional[str],
    cbar_label: Optional[str],
    dpi: int,
) -> None:
    """Render standard or unfolded/projected band structure to an image file."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize

    spins = _resolve_spins(bs.spins, spin_arg)
    all_e = np.concatenate(
        [
            pt.energies - energy_zero
            for bset in bs.sets
            for sp in spins
            for pt in bset.points_by_spin[sp]
        ]
    )
    lo = emin if emin is not None else float(np.min(all_e)) - 0.5
    hi = emax if emax is not None else float(np.max(all_e)) + 0.5
    ylabel = r"$E - E_{\mathrm{F}}$ [eV]" if has_zero_ref else r"$E$ [eV]"

    if not bs.is_unfolded:
        fig, ax = plt.subplots(figsize=(6.4, 5.0), layout="constrained")
        spin_styles = {
            1: {"color": "#1f4e79", "linestyle": "-", "label": "Spin 1"},
            2: {"color": "#c0392b", "linestyle": "--", "label": "Spin 2"},
        }
        for sp in spins:
            style = spin_styles.get(
                sp, {"color": "C0", "linestyle": "-", "label": f"Spin {sp}"}
            )
            labeled = False
            for bset, x_arr in zip(bs.sets, geom.set_x):
                e_mat = (
                    np.array([pt.energies for pt in bset.points_by_spin[sp]])
                    - energy_zero
                )
                for ib in range(e_mat.shape[1]):
                    ax.plot(
                        x_arr,
                        e_mat[:, ib],
                        color=style["color"],
                        linestyle=style["linestyle"],
                        linewidth=1.2,
                        label=(
                            style["label"] if (len(spins) > 1 and not labeled) else None
                        ),
                    )
                    labeled = True
        _decorate_axis(ax, geom, lo, hi, has_zero_ref)
        ax.set_ylabel(ylabel, fontsize=11)
        if len(spins) > 1:
            ax.legend(frameon=False, fontsize=9.5, loc="upper right")
        if title:
            ax.set_title(title, fontsize=11)
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return

    group = _resolve_group(group_arg)
    n_panels = len(spins)
    fig, axes = plt.subplots(
        1,
        n_panels,
        figsize=(6.4 * n_panels, 5.0),
        sharey=True,
        squeeze=False,
        layout="constrained",
    )
    axes_list = list(axes[0])
    default_cbar = f"Weight (group {group})" if group is not None else "Spectral weight"
    cb_text = cbar_label if cbar_label is not None else default_cbar

    if mode == "spectrum":
        e_grid = np.linspace(lo, hi, DEFAULT_ENERGY_GRID, dtype=np.float64)
        panel_maps: list[list[tuple[FloatArray, FloatArray]]] = []
        all_vals: list[FloatArray] = []
        for sp in spins:
            maps_sp: list[tuple[FloatArray, FloatArray]] = []
            for bset, x_arr in zip(bs.sets, geom.set_x):
                spec = compute_spectral_map(
                    bset.points_by_spin[sp], energy_zero, e_grid, sigma_ev, group
                )
                maps_sp.append((x_arr, spec))
                all_vals.append(spec.ravel())
            panel_maps.append(maps_sp)

        p_vmax = (
            vmax
            if vmax is not None
            else max(
                float(np.percentile(np.concatenate(all_vals), SPECTRUM_PERCENTILE_MAX)),
                1.0e-6,
            )
        )
        mappable: Any = None
        for ax, sp, maps_sp in zip(axes_list, spins, panel_maps):
            for x_arr, spec in maps_sp:
                mappable = ax.pcolormesh(
                    x_arr,
                    e_grid,
                    spec,
                    cmap=cmap,
                    vmin=0.0,
                    vmax=p_vmax,
                    shading="gouraud" if len(x_arr) >= 2 else "nearest",
                    rasterized=True,
                )
            _decorate_axis(ax, geom, lo, hi, has_zero_ref)
            if len(spins) > 1:
                ax.set_title(f"Spin {sp}", fontsize=10.5)
    elif mode == "lines":
        p_vmax = vmax if vmax is not None else 1.0
        cmap_obj = plt.get_cmap(cmap)
        norm = Normalize(vmin=0.0, vmax=p_vmax)
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
        mappable.set_array([])

        for ax, sp in zip(axes_list, spins):
            for bset, x_arr in zip(bs.sets, geom.set_x):
                pts = bset.points_by_spin[sp]
                e_mat = np.array([pt.energies for pt in pts]) - energy_zero
                w_mat = np.array([select_weights(pt, group) for pt in pts])
                band_order = np.argsort(w_mat.max(axis=0))

                for ib in band_order:
                    eb = e_mat[:, ib]
                    if eb.max() < lo or eb.min() > hi:
                        continue
                    wb = w_mat[:, ib]
                    pts_xy = np.column_stack([x_arr, eb]).reshape(-1, 1, 2)
                    segs: Any = np.concatenate([pts_xy[:-1], pts_xy[1:]], axis=1)
                    lc = LineCollection(
                        segs, cmap=cmap_obj, norm=norm, linewidths=1.1, alpha=0.88
                    )
                    lc.set_array(0.5 * (wb[:-1] + wb[1:]))
                    ax.add_collection(lc)

                if dot_stride > 0 and len(x_arr) > 0:
                    idx_pts = np.unique(
                        np.concatenate(
                            [np.arange(0, len(x_arr), dot_stride), [len(x_arr) - 1]]
                        )
                    )
                    xs = np.repeat(x_arr[idx_pts], e_mat.shape[1])
                    ys = e_mat[idx_pts, :].ravel()
                    ws = w_mat[idx_pts, :].ravel()
                    mask = (ys >= lo) & (ys <= hi)
                    if np.any(mask):
                        order = np.argsort(ws[mask])
                        ax.scatter(
                            xs[mask][order],
                            ys[mask][order],
                            c=ws[mask][order],
                            cmap=cmap_obj,
                            norm=norm,
                            s=11.0,
                            linewidths=0.0,
                            zorder=4,
                        )
            _decorate_axis(ax, geom, lo, hi, has_zero_ref)
            if len(spins) > 1:
                ax.set_title(f"Spin {sp}", fontsize=10.5)
    else:
        raise ValueError(f"Unsupported --mode '{mode}'.")

    axes_list[0].set_ylabel(ylabel, fontsize=11)
    if title:
        if len(axes_list) == 1:
            axes_list[0].set_title(title, fontsize=10.5)
        else:
            fig.suptitle(title, fontsize=11)
    if mappable is not None:
        cbar = fig.colorbar(mappable, ax=axes_list, pad=0.02, shrink=0.98)
        cbar.set_label(cb_text, fontsize=10)
        cbar.ax.tick_params(labelsize=8.5)

    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Plot and export CP2K band structure (.bs) files."
    )
    parser.add_argument("bs_file", type=Path, help="Path to CP2K .bs file.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output image path (e.g. bands.png).",
    )
    parser.add_argument(
        "--cp2k-out",
        type=Path,
        default=None,
        help="CP2K .out file to read Fermi energy and lattice vectors.",
    )
    parser.add_argument(
        "--efermi",
        type=float,
        default=None,
        help="Explicit Fermi energy [eV] (overrides --cp2k-out).",
    )
    parser.add_argument("--emin", type=float, default=None, help="Minimum energy [eV].")
    parser.add_argument("--emax", type=float, default=None, help="Maximum energy [eV].")
    parser.add_argument(
        "--spin",
        choices=("both", "1", "2"),
        default="both",
        help="Spin channel(s) to plot (default: both).",
    )
    parser.add_argument(
        "--group",
        default="total",
        help="Weight to plot for unfolded .bs: 'total' or 1-based ATOMS group index.",
    )
    parser.add_argument(
        "--mode",
        choices=("spectrum", "lines"),
        default="spectrum",
        help="Plot mode for unfolded/projected .bs: Gaussian 'spectrum' or colored 'lines'.",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=DEFAULT_SIGMA_EV,
        help=f"Gaussian broadening [eV] for --mode spectrum (default: {DEFAULT_SIGMA_EV}).",
    )
    parser.add_argument(
        "--dot-stride",
        type=int,
        default=DEFAULT_DOT_STRIDE,
        help=f"Overlay one dot every N k-points in --mode lines; 0 disables dots (default: {DEFAULT_DOT_STRIDE}).",
    )
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap name (default: viridis).",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Upper color-scale limit (default: 99.5th percentile in spectrum, 1.0 in lines).",
    )
    parser.add_argument("--title", default=None, help="Plot title.")
    parser.add_argument("--cbar-label", default=None, help="Colorbar label.")
    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help=f"Output DPI (default: {DEFAULT_DPI}).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Export parsed band data to a tidy CSV file.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entry point."""
    args = build_arg_parser().parse_args(argv)
    bs = parse_bs_file(args.bs_file)
    out_info = (
        parse_cp2k_out(args.cp2k_out)
        if args.cp2k_out is not None
        else Cp2kOutInfo(fermi_ev=None, cell_angstrom=None)
    )
    geom = compute_kpath_geometry(bs, out_info.cell_angstrom)

    if args.efermi is not None:
        energy_zero, has_zero_ref = args.efermi, True
    elif out_info.fermi_ev is not None:
        energy_zero, has_zero_ref = out_info.fermi_ev, True
    else:
        energy_zero, has_zero_ref = 0.0, False

    if args.csv is not None:
        export_csv(bs, geom, energy_zero, args.csv)

    out_img = args.output
    if out_img is None and args.csv is None:
        out_img = args.bs_file.with_suffix(".png")

    if out_img is not None:
        plot_bands(
            bs=bs,
            geom=geom,
            energy_zero=energy_zero,
            has_zero_ref=has_zero_ref,
            out_path=out_img,
            emin=args.emin,
            emax=args.emax,
            spin_arg=args.spin,
            group_arg=args.group,
            mode=args.mode,
            sigma_ev=args.sigma,
            dot_stride=args.dot_stride,
            cmap=args.cmap,
            vmax=args.vmax,
            title=args.title,
            cbar_label=args.cbar_label,
            dpi=args.dpi,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
