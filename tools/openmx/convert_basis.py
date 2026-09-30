#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Approximate selected OpenMX numerical PAOs by CP2K Gaussian contractions.

This is a radial-function fit, NOT a lossless format conversion or a density-
matrix/restart converter. Requires NumPy and SciPy. See README.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
import re
import sys

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.linalg import lstsq
from scipy.optimize import minimize
from scipy.special import gamma

LETTERS = "spdfg"
ELEMENTS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn "
    "Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
    "In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu "
    "Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
# Indices of OpenMX real harmonics in CP2K's m=-l,...,+l order.
# Verified against Get_Orbitals.c and orbital_transformation_matrices.F.
ANGULAR_PERMUTATIONS = {
    0: [0],
    1: [1, 2, 0],
    2: [2, 4, 0, 3, 1],
    3: [6, 4, 2, 0, 1, 3, 5],
    4: [8, 6, 4, 2, 0, 1, 3, 5, 7],
}


@dataclass(frozen=True)
class BasisSpec:
    stem: str
    counts: tuple[int, ...]

    @classmethod
    def parse(cls, text: str) -> "BasisSpec":
        match = re.fullmatch(
            r"([A-Za-z][A-Za-z0-9_.]*)-((?:[spdfg][1-9][0-9]*)+)", text
        )
        if match is None or ".." in text:
            raise ValueError(
                f"Invalid PAO selection {text!r}; expected e.g. C6.0-s2p2d1"
            )
        counts = [0] * len(LETTERS)
        last = -1
        for letter, count in re.findall(r"([spdfg])([0-9]+)", match[2]):
            angular = LETTERS.index(letter)
            if angular <= last:
                raise ValueError("Shells must be unique and ordered s,p,d,f,g")
            counts[angular] = int(count)
            last = angular
        return cls(match[1], tuple(counts[: last + 1]))

    @property
    def label(self) -> str:
        return (
            self.stem
            + "-"
            + "".join(f"{LETTERS[l]}{n}" for l, n in enumerate(self.counts) if n)
        )

    @property
    def name(self) -> str:
        return "OMX-FIT-" + self.label

    @property
    def nao(self) -> int:
        return sum(n * (2 * l + 1) for l, n in enumerate(self.counts))

    def permutation(self) -> list[int]:
        result = []
        offset = 0
        for l, count in enumerate(self.counts):
            for _ in range(count):
                result.extend(offset + m for m in ANGULAR_PERMUTATIONS[l])
                offset += 2 * l + 1
        return result


@dataclass
class PaoFile:
    path: Path
    element: str
    cutoff: float
    valence: float
    radius: np.ndarray
    orbitals: dict[int, np.ndarray]
    sha256: str

    @classmethod
    def read(cls, path: Path) -> "PaoFile":
        raw = path.read_bytes()
        text = raw.decode("utf-8")

        def scalar(key: str, integer: bool = False) -> float | int:
            matches = re.findall(
                r"^\s*" + re.escape(key) + r"\s+([^\s#]+)", text, re.M | re.I
            )
            if len(matches) != 1:
                raise ValueError(f"{path}: expected exactly one {key}")
            value = float(matches[0].replace("D", "E").replace("d", "e"))
            if not math.isfinite(value) or (integer and value != int(value)):
                raise ValueError(f"{path}: invalid {key}")
            return int(value) if integer else value

        atomic_number = scalar("AtomSpecies", True)
        lmax, multiplicity = scalar("PAO.Lmax", True), scalar("PAO.Mul", True)
        rows = scalar("grid.num.output", True)
        cutoff = scalar("radial.cutoff.pao")
        valence = scalar("valence.electron")
        if not 1 <= atomic_number < len(ELEMENTS) or not 0 <= lmax < len(LETTERS):
            raise ValueError(
                f"{path}: unsupported element or angular momentum (supported: s..g)"
            )
        if multiplicity < 1 or rows < 6 or cutoff <= 0 or valence < 0:
            raise ValueError(f"{path}: invalid PAO dimensions or cutoff")
        orbitals = {}
        radius = None
        for l in range(lmax + 1):
            tag = f"pseudo.atomic.orbitals.L={l}"
            blocks = re.findall(
                r"^\s*<"
                + re.escape(tag)
                + r"\s*\n(.*?)^\s*"
                + re.escape(tag)
                + r">\s*$",
                text,
                re.M | re.S,
            )
            if len(blocks) != 1:
                raise ValueError(f"{path}: missing or duplicated block {tag}")
            table = np.array(
                [
                    [
                        float(x.replace("D", "E").replace("d", "e"))
                        for x in line.split("#", 1)[0].split()
                    ]
                    for line in blocks[0].splitlines()
                    if line.split("#", 1)[0].strip()
                ]
            )
            if table.shape != (rows, multiplicity + 2) or not np.isfinite(table).all():
                raise ValueError(
                    f"{path}: invalid dimensions/nonfinite values in {tag}"
                )
            r = table[:, 1]
            if np.any(r <= 0) or np.any(np.diff(r) <= 0):
                raise ValueError(f"{path}: radial grid must be positive and increasing")
            if not np.allclose(np.log(r), table[:, 0], atol=1e-8, rtol=0):
                raise ValueError(f"{path}: inconsistent log(r), r columns")
            if radius is not None and not np.array_equal(radius, r):
                raise ValueError(
                    f"{path}: angular channels have different radial grids"
                )
            radius = r
            orbitals[l] = table[:, 2:]
        if radius[-1] < cutoff:
            raise ValueError(
                f"{path}: grid does not reach the nominal confinement radius"
            )
        return cls(
            path,
            ELEMENTS[atomic_number],
            cutoff,
            valence,
            radius,
            orbitals,
            hashlib.sha256(raw).hexdigest(),
        )

    def radial(
        self, l: int, count: int, r: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Cubic-Hermite interpolation; regular origin and zero outside PAO grid.

        Interior slopes are the same three-point derivatives as Get_Orbitals.c.
        At the tiny inner boundary a regular r**l extrapolation is used. This
        approximation (also relevant to derivatives) is reported explicitly.
        """
        if l not in self.orbitals or count > self.orbitals[l].shape[1]:
            raise ValueError(
                f"{self.path}: requested {LETTERS[l]}{count} is unavailable"
            )
        values = self.orbitals[l][:, :count]
        slopes = np.gradient(values, self.radius, axis=0, edge_order=2)
        spline = CubicHermiteSpline(
            self.radius, values, slopes, axis=0, extrapolate=False
        )
        clipped = np.clip(r, self.radius[0], self.radius[-1])
        value, derivative = spline(clipped), spline(clipped, 1)
        low, high = r < self.radius[0], r > self.radius[-1]
        value[low] = (r[low, None] / self.radius[0]) ** l * values[0]
        derivative[low] = (
            l * r[low, None] ** (l - 1) / self.radius[0] ** l * values[0] if l else 0.0
        )
        value[high], derivative[high] = 0.0, 0.0
        return value, derivative


def quadrature(radius: np.ndarray, order: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre integration on each radial interval (weights are dr)."""
    x, w = np.polynomial.legendre.leggauss(order)
    widths = np.diff(radius)
    r = (radius[:-1, None] + (x + 1) * widths[:, None] / 2).ravel()
    weights = (widths[:, None] * w / 2).ravel()
    return r, weights


def primitives(
    r: np.ndarray, exponents: np.ndarray, l: int
) -> tuple[np.ndarray, np.ndarray]:
    """Normalized spherical radial Gaussians, integral r**2*g**2 dr = 1."""
    normalization = np.sqrt(2 * (2 * exponents) ** (l + 1.5) / gamma(l + 1.5))
    g = normalization * r[:, None] ** l * np.exp(-r[:, None] ** 2 * exponents)
    dg = g * (l / r[:, None] - 2 * r[:, None] * exponents)
    return g, dg


def gaussian_overlap(exponents: np.ndarray, l: int) -> np.ndarray:
    a, b = exponents[:, None], exponents[None, :]
    return (2 * np.sqrt(a * b) / (a + b)) ** (l + 1.5)


@dataclass
class Fit:
    l: int
    exponents: np.ndarray
    coefficients: np.ndarray
    diagnostics: dict


def fit_channel(
    pao: PaoFile,
    l: int,
    count: int,
    nprimitive: int = 24,
    derivative_weight: float = 0.2,
    rcond: float = 1e-6,
) -> Fit:
    """Variable projection with optimized even-tempered exponent endpoints.

    Shared exponents, independent contractions: no shell mixing or change in AO
    count. The objective includes radial L2 and gradient errors, including tails.
    """
    amin = 0.4 / pao.cutoff**2
    tail_end = max(pao.radius[-1] * 3, math.sqrt(45 / amin))
    knots = np.concatenate(
        ([0.0], pao.radius, np.geomspace(pao.radius[-1], tail_end, 160)[1:])
    )
    r, w = quadrature(knots, 2)
    target, dtarget = pao.radial(l, count, r)
    weight = np.sqrt(w) * r
    norm = np.sqrt(np.sum((weight[:, None] * target) ** 2, axis=0))
    if np.any(norm < 1e-12):
        raise ValueError("Cannot fit a zero radial orbital")
    target, dtarget = target / norm, dtarget / norm
    angular_weight = np.sqrt(w * l * (l + 1))
    rhs = np.vstack(
        (
            weight[:, None] * target,
            derivative_weight * weight[:, None] * dtarget,
            derivative_weight * angular_weight[:, None] * target,
        )
    )

    def solve(
        parameters: np.ndarray,
    ) -> tuple[float, np.ndarray, np.ndarray, int, np.ndarray]:
        exponents = np.geomspace(
            math.exp(parameters[0]), math.exp(parameters[1]), nprimitive
        )
        g, dg = primitives(r, exponents, l)
        design = np.vstack(
            (
                weight[:, None] * g,
                derivative_weight * weight[:, None] * dg,
                derivative_weight * angular_weight[:, None] * g,
            )
        )
        coefficients, _, rank, singular = lstsq(
            design, rhs, cond=rcond, lapack_driver="gelsd"
        )
        error = design @ coefficients - rhs
        return float(np.sum(error**2)), exponents, coefficients, rank, singular

    initial = np.log([1.0 / pao.cutoff**2, 100.0])
    bounds = [
        (math.log(amin), math.log(10 / pao.cutoff**2)),
        (math.log(20.0), math.log(4000.0)),
    ]
    optimization = minimize(
        lambda p: solve(p)[0],
        initial,
        method="Nelder-Mead",
        bounds=bounds,
        options={"maxiter": 90, "xatol": 1e-4, "fatol": 1e-12},
    )
    choices = [solve(initial), solve(optimization.x)]
    _, alpha, coeff, rank, singular = min(choices, key=lambda item: item[0])

    # Independent, higher-order quadrature for exported (CP2K-normalized) orbitals.
    rv, wv = quadrature(knots, 6)
    original, original_d = pao.radial(l, count, rv)
    target_norm = np.sqrt(np.sum(wv[:, None] * rv[:, None] ** 2 * original**2, axis=0))
    original, original_d = original / target_norm, original_d / target_norm
    gv, dgv = primitives(rv, alpha, l)
    overlap = gaussian_overlap(alpha, l)
    analytic_norm2 = np.einsum("ij,ik,kj->j", coeff, overlap, coeff)
    if np.any(analytic_norm2 <= 0) or not np.isfinite(analytic_norm2).all():
        raise ValueError("Ill-conditioned Gaussian contraction: non-positive norm")
    coeff /= np.sqrt(analytic_norm2)
    fitted, fitted_d = gv @ coeff, dgv @ coeff
    norm_check = np.sum(wv[:, None] * rv[:, None] ** 2 * fitted**2, axis=0)
    error = np.sqrt(
        np.sum(wv[:, None] * rv[:, None] ** 2 * (fitted - original) ** 2, axis=0)
    )
    kinetic_target = 0.5 * np.sum(
        wv[:, None] * (rv[:, None] ** 2 * original_d**2 + l * (l + 1) * original**2),
        axis=0,
    )
    kinetic_fit = 0.5 * np.sum(
        wv[:, None] * (rv[:, None] ** 2 * fitted_d**2 + l * (l + 1) * fitted**2), axis=0
    )
    beyond = rv > pao.cutoff
    leakage = np.sum(
        wv[beyond, None] * rv[beyond, None] ** 2 * fitted[beyond] ** 2, axis=0
    )
    source_overlap = original.T @ ((wv * rv**2)[:, None] * original)
    fit_overlap = coeff.T @ overlap @ coeff
    diagnostics = {
        "l": l,
        "n_radial": count,
        "n_primitive": nprimitive,
        "source_radial_norms": target_norm.tolist(),
        "normalized_radial_l2_errors": error.tolist(),
        "kinetic_source_hartree": kinetic_target.tolist(),
        "kinetic_fit_hartree": kinetic_fit.tolist(),
        "kinetic_absolute_errors_hartree": np.abs(
            kinetic_fit - kinetic_target
        ).tolist(),
        "norm_outside_nominal_cutoff": leakage.tolist(),
        "same_center_overlap_max_error": float(
            np.max(np.abs(fit_overlap - source_overlap))
        ),
        "analytic_vs_quadrature_norm_max_error": float(np.max(np.abs(norm_check - 1))),
        "max_abs_coefficient": float(np.max(np.abs(coeff))),
        "least_squares_rank": int(rank),
        "smallest_relative_singular_value": float(singular[-1] / singular[0]),
        "exponent_search_converged": bool(optimization.success),
        "svd_relative_cutoff": rcond,
    }
    return Fit(l, alpha[::-1], coeff[::-1], diagnostics)


def render_basis(spec: BasisSpec, pao: PaoFile, fits: list[Fit]) -> str:
    lines = [
        f"# Source: {pao.path.name}; SHA256 {pao.sha256}",
        f"# Nominal PAO cutoff {pao.cutoff:g} bohr; OpenMX valence {pao.valence:g}",
        f"# Approximate radial fit; NOT a pseudopotential conversion; {spec.nao} AOs.",
        f"# Maximum normalized radial L2 error: {max(max(f.diagnostics['normalized_radial_l2_errors']) for f in fits):.8e}",
        f"{pao.element} {spec.name}",
        str(len(fits)),
    ]
    for fit in fits:
        lines.append(
            f"{fit.l + 1} {fit.l} {fit.l} {len(fit.exponents)} {fit.coefficients.shape[1]}"
        )
        for exponent, coefficients in zip(fit.exponents, fit.coefficients):
            lines.append(" ".join(f"{v:.16e}" for v in [exponent, *coefficients]))
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "selections", nargs="+", help="Exact OpenMX selections, e.g. C6.0-s2p2d1"
    )
    parser.add_argument("--pao-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("BASIS_OMX"))
    parser.add_argument(
        "--report", type=Path, help="JSON diagnostics; default OUTPUT.json"
    )
    parser.add_argument("--nprimitive", type=int, default=24)
    parser.add_argument(
        "--derivative-weight",
        type=float,
        default=0.2,
        help="Gradient-error weight in bohr",
    )
    parser.add_argument("--max-l2-error", type=float, default=0.01)
    parser.add_argument(
        "--allow-inaccurate",
        action="store_true",
        help="Explicitly allow fits above error threshold",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    report_path = args.report or args.output.with_name(args.output.name + ".json")
    if args.output.resolve() == report_path.resolve():
        parser.error("Basis and report paths must differ")
    if not args.overwrite and (args.output.exists() or report_path.exists()):
        parser.error("Output exists; use --overwrite explicitly")
    if (
        not 4 <= args.nprimitive <= 80
        or not math.isfinite(args.derivative_weight)
        or args.derivative_weight < 0
    ):
        parser.error("Need 4..80 primitives and a finite nonnegative derivative weight")
    if not math.isfinite(args.max_l2_error) or args.max_l2_error <= 0:
        parser.error("--max-l2-error must be finite and positive")
    try:
        specs = [BasisSpec.parse(s) for s in args.selections]
        if len({s.label for s in specs}) != len(specs):
            raise ValueError("Duplicate basis selection")
        source_paths = {(args.pao_dir / (s.stem + ".pao")).resolve() for s in specs}
        if (
            args.output.resolve() in source_paths
            or report_path.resolve() in source_paths
        ):
            raise ValueError("An output path would overwrite a PAO source file")
        reports, entries = [], []
        for spec in specs:
            pao = PaoFile.read(args.pao_dir / (spec.stem + ".pao"))
            filename_element = re.match(r"[A-Z][a-z]?", spec.stem)
            if filename_element is None or filename_element[0] != pao.element:
                raise ValueError(f"{spec.label}: filename and AtomSpecies disagree")
            fits = [
                fit_channel(pao, l, n, args.nprimitive, args.derivative_weight)
                for l, n in enumerate(spec.counts)
                if n
            ]
            maximum = max(
                max(f.diagnostics["normalized_radial_l2_errors"]) for f in fits
            )
            norm_error = max(
                f.diagnostics["analytic_vs_quadrature_norm_max_error"] for f in fits
            )
            if norm_error > 1e-6:
                raise ValueError(
                    f"{spec.label}: cancellation in Gaussian normalization ({norm_error:.3g})"
                )
            if maximum > args.max_l2_error and not args.allow_inaccurate:
                raise ValueError(
                    f"{spec.label}: L2 error {maximum:.3g} exceeds {args.max_l2_error:g}; "
                    "increase primitives or explicitly use --allow-inaccurate"
                )
            entries.append(render_basis(spec, pao, fits))
            reports.append(
                {
                    "selection": spec.label,
                    "cp2k_name": spec.name,
                    "element": pao.element,
                    "source_file": pao.path.name,
                    "source_sha256": pao.sha256,
                    "source_valence": pao.valence,
                    "nominal_cutoff_bohr": pao.cutoff,
                    "nao": spec.nao,
                    "cp2k_from_openmx_indices_zero_based": spec.permutation(),
                    "radial_fits": [f.diagnostics for f in fits],
                }
            )
            print(
                f"{spec.label}: {spec.nao} AOs; maximum radial L2 error {maximum:.6g}",
                flush=True,
            )
        header = (
            "# BASIS_OMX: approximate Gaussian fits to selected OpenMX PAOs.\n"
            "# Generated with tools/openmx/convert_basis.py. Default CP2K orbital normalization.\n"
            "# OpenMX database 2019, T. Ozaki and H. Kawai, distributed under GNU GPL.\n"
            "# https://www.openmx-square.org/vps_pao2019/\n"
            "# EXPERIMENTAL: NOT lossless and NOT validated for production energies/forces.\n"
            "# Not the full database. Exact selections and errors: companion JSON.\n"
            "# AO order differs; a raw OpenMX density matrix is NOT a CP2K restart.\n\n"
        )
        report = {
            "schema": "openmx-cp2k-radial-fit-v1",
            "exact_conversion": False,
            "cp2k_normalization": "default normalized primitives and contractions (norm_type=2)",
            "radial_interpolation": "cubic Hermite, three-point interior slopes, regular r**l origin",
            "derivative_weight_bohr": args.derivative_weight,
            "max_l2_error_requested": args.max_l2_error,
            "basis_sha256": hashlib.sha256(
                (header + "\n".join(entries)).encode()
            ).hexdigest(),
            "basis_sets": reports,
        }
        args.output.write_text(header + "\n".join(entries))
        report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Error: {exc}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
