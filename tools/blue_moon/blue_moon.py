#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
#   CP2K: A general program to perform molecular dynamics simulations
#   Copyright 2000-2026 CP2K developers group <https://cp2k.org>
#
#   SPDX-License-Identifier: GPL-2.0-or-later
# --------------------------------------------------------------------------------------------------

"""Postprocess ONE fixed collective constraint; see README.md for limitations."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence, TextIO, Union, cast

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]
Frame = tuple[int, FloatArray]
Sample = tuple[int, float, float, float, float, float, float]

# Match src/common/physcon.F, not a different vintage of CODATA constants.
ANGSTROM_PER_BOHR = 0.52917720859
KB_HARTREE_PER_K = 1.3806504e-23 / (2 * 10973731.568527 * 6.62606896e-34 * 299792458)


@dataclass
class Derivative:
    """Scalar, Cartesian gradient and Hessian (second-order forward differentiation)."""

    value: float
    gradient: FloatArray
    hessian: FloatArray

    def __add__(self, other: Union[Derivative, float]) -> Derivative:
        if not isinstance(other, Derivative):
            return Derivative(self.value + other, self.gradient, self.hessian)
        return Derivative(
            self.value + other.value,
            self.gradient + other.gradient,
            self.hessian + other.hessian,
        )

    __radd__ = __add__

    def __neg__(self) -> Derivative:
        return self * -1

    def __sub__(self, other: Union[Derivative, float]) -> Derivative:
        return self + (-other)

    def __rsub__(self, other: Union[Derivative, float]) -> Derivative:
        return -self + other

    def __mul__(self, other: Union[Derivative, float]) -> Derivative:
        if not isinstance(other, Derivative):
            return Derivative(
                self.value * other, self.gradient * other, self.hessian * other
            )
        cross = np.outer(self.gradient, other.gradient)
        return Derivative(
            self.value * other.value,
            self.gradient * other.value + self.value * other.gradient,
            self.hessian * other.value + self.value * other.hessian + cross + cross.T,
        )

    __rmul__ = __mul__

    def compose(self, value: float, first: float, second: float) -> Derivative:
        return Derivative(
            value,
            first * self.gradient,
            first * self.hessian + second * np.outer(self.gradient, self.gradient),
        )

    def __truediv__(self, other: Union[Derivative, float]) -> Derivative:
        if not isinstance(other, Derivative):
            return self * (1 / other)
        x = other.value
        return self * other.compose(1 / x, -1 / x**2, 2 / x**3)

    def sqrt(self) -> Derivative:
        if self.value <= 1e-24:
            raise ValueError("Degenerate CV: zero distance or collinear plane/torsion")
        root = math.sqrt(self.value)
        return self.compose(root, 0.5 / root, -0.25 / root**3)


def atan2(y: Derivative, x: Derivative) -> Derivative:
    """Differentiate atan2 without dividing by x (which may be zero)."""
    r2 = x.value**2 + y.value**2
    if r2 <= 1e-24:
        raise ValueError("Degenerate angular CV")
    fx, fy = -y.value / r2, x.value / r2
    fxx, fyy = 2 * x.value * y.value / r2**2, -2 * x.value * y.value / r2**2
    fxy = (y.value**2 - x.value**2) / r2**2
    cross = np.outer(x.gradient, y.gradient)
    return Derivative(
        math.atan2(y.value, x.value),
        fx * x.gradient + fy * y.gradient,
        fx * x.hessian
        + fy * y.hessian
        + fxx * np.outer(x.gradient, x.gradient)
        + fyy * np.outer(y.gradient, y.gradient)
        + fxy * (cross + cross.T),
    )


def dot(a: Sequence[Derivative], b: Sequence[Derivative]) -> Derivative:
    result = a[0] * b[0]
    for x, y in zip(a[1:], b[1:]):
        result = result + x * y
    return result


def cross(a: Sequence[Derivative], b: Sequence[Derivative]) -> list[Derivative]:
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def subtract(a: Sequence[Derivative], b: Sequence[Derivative]) -> list[Derivative]:
    return [x - y for x, y in zip(a, b)]


def check_keys(
    mapping: object, required: Iterable[str], optional: Iterable[str] = ()
) -> dict[str, object]:
    if not isinstance(mapping, dict) or not all(isinstance(k, str) for k in mapping):
        raise ValueError("Expected a JSON object with string keys")
    # The container and its keys have been checked; values remain untrusted.
    mapping = cast(dict[str, object], mapping)
    required_keys = set(required)
    missing = required_keys - mapping.keys()
    extra = mapping.keys() - required_keys - set(optional)
    if missing or extra:
        raise ValueError(
            f"Invalid keys: missing {sorted(missing)}, unknown {sorted(extra)}"
        )
    return mapping


def finite_number(value: object) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
    ):
        raise ValueError(f"Expected a finite number, got {value!r}")
    return float(value)


@dataclass(frozen=True)
class PrimitiveCV:
    kind: str
    atoms: tuple[int, ...]
    pbc: bool
    reference: float


@dataclass(frozen=True)
class LinearCombinationCV:
    terms: tuple[tuple[float, CV], ...]


CV = Union[PrimitiveCV, LinearCombinationCV]


def parse_cv(value: object) -> CV:
    """Validate raw JSON before creating the typed coordinate tree."""
    cv = check_keys(value, ("type",), ("atoms", "pbc", "reference", "terms"))
    kind = cv["type"]
    if kind == "linear_combination":
        check_keys(cv, ("type", "terms"))
        raw_terms = cv["terms"]
        if not isinstance(raw_terms, list) or not raw_terms:
            raise ValueError("A linear combination needs nonempty terms")
        terms: list[tuple[float, CV]] = []
        for raw_term in raw_terms:
            term = check_keys(raw_term, ("coefficient", "cv"))
            terms.append((finite_number(term["coefficient"]), parse_cv(term["cv"])))
        return LinearCombinationCV(tuple(terms))
    counts = {
        "distance": 2,
        "angle": 3,
        "torsion": 4,
        "point_plane": 4,
        "point_bond_center": 3,
    }
    if not isinstance(kind, str) or kind not in counts:
        raise ValueError(f"Unsupported CV type: {kind}")
    check_keys(
        cv,
        ("type", "atoms"),
        ("pbc", "reference") if kind == "torsion" else ("pbc",),
    )
    raw_atoms = cv["atoms"]
    if not isinstance(raw_atoms, list) or len(raw_atoms) != counts[kind]:
        raise ValueError(f"Wrong atom count for {kind}")
    atoms: list[int] = []
    for atom in raw_atoms:
        if type(atom) is not int or atom < 1 or atom in atoms:
            raise ValueError(
                "Primitive CV atoms must be distinct, positive, one-based integers"
            )
        atoms.append(atom)
    pbc = cv.get("pbc", True)
    if not isinstance(pbc, bool):
        raise ValueError("pbc must be true or false")
    reference = 0.0
    if kind == "torsion":
        # CP2K follows a continuous torsion branch. The branch must be specified,
        # especially when torsions occur inside linear combinations.
        if "reference" not in cv:
            raise ValueError("A torsion needs a reference angle in radians")
        reference = finite_number(cv["reference"])
    return PrimitiveCV(kind, tuple(atoms), pbc, reference)


def cv_atoms(cv: CV) -> set[int]:
    """Shared atoms are differentiated only once, also across nested terms."""
    if isinstance(cv, PrimitiveCV):
        return set(cv.atoms)
    return {atom for _, term in cv.terms for atom in cv_atoms(term)}


class Coordinate:
    def __init__(self, definition: object, cell_bohr: ArrayLike) -> None:
        self.definition = parse_cv(definition)
        self.atoms = sorted(cv_atoms(self.definition))
        self.cell = np.asarray(cell_bohr, dtype=float)
        if self.cell.shape != (3, 3) or not np.isfinite(self.cell).all():
            raise ValueError(
                "cell_angstrom must contain three finite cell vectors (rows)"
            )
        if np.linalg.cond(self.cell) > 1e12:
            raise ValueError("Singular or ill-conditioned cell")
        self.inverse = np.linalg.inv(self.cell)

    def minimum_image(
        self, vector: Sequence[Derivative], enabled: bool
    ) -> list[Derivative]:
        if not enabled:
            return list(vector)
        fractional = np.array([x.value for x in vector]) @ self.inverse
        nearest = np.copysign(np.floor(np.abs(fractional) + 0.5), fractional)
        if np.any(np.abs(np.abs(fractional - nearest) - 0.5) < 1e-8):
            raise ValueError("CV lies on a minimum-image branch boundary")
        shift = nearest @ self.cell
        return [x - s for x, s in zip(vector, shift)]

    def evaluate(self, positions_bohr: ArrayLike) -> Derivative:
        positions = np.asarray(positions_bohr, dtype=float)
        if (
            positions.ndim != 2
            or positions.shape[1] != 3
            or not np.isfinite(positions).all()
        ):
            raise ValueError("Expected finite Cartesian positions, shape (natoms, 3)")
        if len(positions) < max(self.atoms):
            raise ValueError("CV atom index exceeds trajectory atom count")
        size = 3 * len(self.atoms)
        basis = np.eye(size)
        points = {
            atom: [
                Derivative(
                    positions[atom - 1, k], basis[3 * i + k], np.zeros((size, size))
                )
                for k in range(3)
            ]
            for i, atom in enumerate(self.atoms)
        }
        result = self._evaluate(self.definition, points)
        if not all(
            np.isfinite(x).all()
            for x in (result.value, result.gradient, result.hessian)
        ):
            raise ValueError("Non-finite CV or derivatives")
        return result

    def _evaluate(self, cv: CV, points: dict[int, list[Derivative]]) -> Derivative:
        if isinstance(cv, LinearCombinationCV):
            initial_term, *rest = cv.terms
            result = initial_term[0] * self._evaluate(initial_term[1], points)
            for coefficient, term in rest:
                result = result + coefficient * self._evaluate(term, points)
            return result
        kind = cv.kind
        p = [points[i] for i in cv.atoms]

        def mic(v: Sequence[Derivative]) -> list[Derivative]:
            return self.minimum_image(v, cv.pbc)

        if kind == "distance":
            v = mic(subtract(p[0], p[1]))
            return dot(v, v).sqrt()
        if kind == "point_bond_center":
            # CP2K POINT TYPE GEO_CENTER uses the arithmetic mean of coordinates;
            # it does not unwrap the constituent atoms. They must form a whole bond.
            center = [(a + b) / 2 for a, b in zip(p[1], p[2])]
            v = mic(subtract(p[0], center))
            return dot(v, v).sqrt()
        if kind == "point_plane":
            # Ordering: the three ATOMS_PLANE, then ATOM_POINT. Match CP2K's
            # signed normal and its separate wrapping of the centroid displacement.
            a, b = mic(subtract(p[0], p[1])), mic(subtract(p[2], p[1]))
            centroid = [(x + y + z) / 3 for x, y, z in zip(*p[:3])]
            v = mic(subtract(p[3], centroid))
            normal = cross(a, b)
            return dot(v, normal) / dot(normal, normal).sqrt()
        if kind == "angle":
            a, b = mic(subtract(p[0], p[1])), mic(subtract(p[2], p[1]))
            normal = cross(a, b)
            return atan2(dot(normal, normal).sqrt(), dot(a, b))
        if kind == "torsion":
            a, b, c = [mic(subtract(p[i + 1], p[i])) for i in range(3)]
            # CP2K's gradient also wraps the two 1-3 displacements. If they
            # disagree with the sums of the adjacent bonds, its gradient is
            # not the derivative of this local torsion branch. Do not mix them.
            for i, first, second in ((0, a, b), (1, b, c)):
                diagonal = mic(subtract(p[i + 2], p[i]))
                residual = [
                    d.value - x.value - y.value
                    for d, x, y in zip(diagonal, first, second)
                ]
                if not np.allclose(residual, 0, atol=1e-9, rtol=0):
                    raise ValueError(
                        "Inconsistent CP2K torsion images for adjacent and 1-3 bonds"
                    )
            t, u = cross(a, b), cross(b, c)
            dot(t, t).sqrt()
            dot(u, u).sqrt()
            phi = atan2(dot(b, cross(t, u)) / dot(b, b).sqrt(), dot(t, u))
            difference = math.remainder(phi.value - cv.reference, 2 * math.pi)
            if abs(abs(difference) - math.pi) < 1e-8:
                raise ValueError("Torsion lies on the specified branch boundary")
            return phi + (cv.reference + difference - phi.value)
        raise ValueError(f"Unsupported CV type: {kind}")


def metric(derivative: Derivative, masses: ArrayLike) -> tuple[float, float]:
    """Z = g^T M^-1 g; G = (M^-1 g)^T H (M^-1 g) / Z^2.

    A common mass-unit conversion cancels from G and normalized reweighting.
    Z is reported using masses in amu, derivatives in bohr/radians.
    """
    masses = np.asarray(masses, dtype=float)
    if (
        masses.shape != (len(derivative.gradient) // 3,)
        or not np.isfinite(masses).all()
        or np.any(masses <= 0)
    ):
        raise ValueError("Need one finite positive mass per participating atom")
    velocity = derivative.gradient / np.repeat(masses, 3)
    z = float(derivative.gradient @ velocity)
    if not math.isfinite(z) or z <= 0:
        raise ValueError("Singular mass metric: the CV gradient is zero")
    g = float(velocity @ derivative.hessian @ velocity / z**2)
    if not math.isfinite(g):
        raise ValueError("Non-finite metric correction")
    return z, g


def read_xyz(stream: TextIO) -> Iterator[Frame]:
    """Read CP2K XYZ coordinates (Angstrom) with strictly increasing MD steps."""
    previous = -1
    labels = None
    for line in stream:
        if not line.strip():
            continue
        try:
            count = int(line.strip())
        except ValueError:
            raise ValueError("Invalid XYZ atom count") from None
        if count <= 0:
            raise ValueError("XYZ atom count must be positive")
        title = stream.readline()
        match = re.search(r"\bi\s*=\s*(\d+)\s*,", title)
        if not match:
            raise ValueError("XYZ title needs the CP2K MD step: i = N,")
        step = int(match[1])
        if step <= previous:
            raise ValueError(
                "Repeated/decreasing XYZ steps: split restarted runs before analysis"
            )
        previous = step
        frame_labels, rows = [], []
        for _ in range(count):
            fields = stream.readline().split()
            if len(fields) != 4:
                raise ValueError(f"Incomplete/invalid XYZ coordinates at step {step}")
            frame_labels.append(fields[0])
            rows.append([float(x) for x in fields[1:]])
        if labels is not None and labels != frame_labels:
            raise ValueError("Atom count/order/labels changed in XYZ trajectory")
        labels = frame_labels
        positions = np.asarray(rows, dtype=np.float64) / ANGSTROM_PER_BOHR
        if not np.isfinite(positions).all():
            raise ValueError("Non-finite XYZ coordinates")
        yield step, positions


def read_multipliers(stream: TextIO) -> Iterator[float]:
    """Strict velocity-Verlet SHAKE/RATTLE pairs, with exactly one constraint.

    Accept the standard fixed-width output and reject wrapped multi-constraint
    records, missing pairs, overflow stars, and NaN rather than picking a column.
    """
    pending = None
    for line in stream:
        if not line.strip():
            continue
        match = re.fullmatch(
            r"\s*(Shake|Rattle)\s+Lagrangian Multipliers:\s*(.*?)\s*", line
        )
        if not match or len(match[2].split()) != 1:
            raise ValueError("Expected exactly one multiplier per SHAKE/RATTLE record")
        value = float(match[2].replace("D", "E").replace("d", "e"))
        if not math.isfinite(value):
            raise ValueError("Non-finite Lagrange multiplier")
        if match[1] == "Shake":
            if pending is not None:
                raise ValueError("Missing RATTLE record after SHAKE")
            pending = value
        else:
            if pending is None:
                raise ValueError(
                    "RATTLE without SHAKE: check initialization/integrator"
                )
            yield pending
            pending = None
    if pending is not None:
        raise ValueError("Incomplete final SHAKE/RATTLE pair")


def analyze(
    config: object,
    frames: Iterable[Frame],
    multipliers: Iterable[float],
    first_step: int,
    stride: int = 1,
) -> Iterator[Sample]:
    settings = check_keys(
        config,
        (
            "cv",
            "cell_angstrom",
            "masses_amu",
            "temperature_kelvin",
            "target_au",
            "target_tolerance_au",
        ),
    )
    coordinate = Coordinate(
        settings["cv"],
        np.asarray(settings["cell_angstrom"], dtype=float) / ANGSTROM_PER_BOHR,
    )
    mass_values = check_keys(settings["masses_amu"], [str(i) for i in coordinate.atoms])
    masses = [finite_number(mass_values[str(i)]) for i in coordinate.atoms]
    temperature = finite_number(settings["temperature_kelvin"])
    target = finite_number(settings["target_au"])
    tolerance = finite_number(settings["target_tolerance_au"])
    if temperature <= 0 or tolerance <= 0 or first_step < 0 or stride < 1:
        raise ValueError(
            "Temperature, tolerance and stride must be positive; first step nonnegative"
        )
    selected = (frame for frame in frames if frame[0] >= first_step)
    for n, pair in enumerate(itertools.zip_longest(selected, multipliers)):
        frame, multiplier = pair
        if frame is None or multiplier is None:
            raise ValueError(
                "Trajectory/multiplier count mismatch; check initialization, print cadence and restarts"
            )
        step, positions = frame
        expected = first_step + n * stride
        if step != expected:
            raise ValueError(
                f"Expected XYZ step {expected}, got {step}; no automatic pairing across gaps"
            )
        derivative = coordinate.evaluate(positions)
        if abs(derivative.value - target) > tolerance:
            raise ValueError(
                f"CV at step {step} ({derivative.value}) differs from fixed target {target}; check definition, units and sampling"
            )
        z, g = metric(derivative, masses)
        corrected = -multiplier + KB_HARTREE_PER_K * temperature * g
        if not math.isfinite(corrected):
            raise ValueError("Non-finite corrected force")
        yield step, derivative.value, multiplier, z, g, z**-0.5, corrected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "config",
        type=Path,
        help="JSON CV, actual masses, fixed cell, temperature and target",
    )
    parser.add_argument("trajectory", type=Path, help="CP2K XYZ trajectory in Angstrom")
    parser.add_argument(
        "multipliers", type=Path, help="CP2K LagrangeMultLog, ONE constraint only"
    )
    parser.add_argument(
        "--first-step",
        type=int,
        required=True,
        help="MD step of the FIRST SHAKE record (not inferred from XYZ)",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Identical MD print stride in BOTH input files (default: 1)",
    )
    parser.add_argument(
        "--discard",
        type=int,
        default=0,
        help="Discard this many paired samples as equilibration",
    )
    args = parser.parse_args()
    if args.discard < 0:
        parser.error("--discard must be nonnegative")
    try:
        with args.config.open() as stream:
            config: object = json.load(stream)
        with args.trajectory.open() as xyz, args.multipliers.open() as lagrange:
            samples = analyze(
                config,
                read_xyz(xyz),
                read_multipliers(lagrange),
                args.first_step,
                args.stride,
            )
            writer = csv.writer(sys.stdout)
            writer.writerow(
                (
                    "step",
                    "cv_au",
                    "lambda_au",
                    "Z_amu",
                    "G_au",
                    "weight",
                    "corrected_force_au",
                )
            )
            numerator, denominator, raw = [], [], []
            for n, sample in enumerate(samples):
                if n < args.discard:
                    continue
                writer.writerow(sample)
                numerator.append(sample[5] * sample[6])
                denominator.append(sample[5])
                raw.append(-sample[2])
            if not numerator:
                raise ValueError("No production samples remain")
            summary = {
                "samples": len(numerator),
                "free_energy_gradient_au": math.fsum(numerator)
                / math.fsum(denominator),
                "uncorrected_minus_mean_lambda_au": math.fsum(raw) / len(raw),
            }
            print(json.dumps(summary, indent=2), file=sys.stderr)
    except (
        OSError,
        ValueError,
        TypeError,
        KeyError,
        OverflowError,
        ZeroDivisionError,
    ) as error:
        parser.exit(
            1,
            f"blue_moon: {error}\nNo valid estimate produced; discard any partial CSV.\n",
        )


if __name__ == "__main__":
    main()
