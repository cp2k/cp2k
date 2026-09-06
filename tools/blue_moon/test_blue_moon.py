#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
#   CP2K: A general program to perform molecular dynamics simulations
#   Copyright 2000-2026 CP2K developers group <https://cp2k.org>
#
#   SPDX-License-Identifier: GPL-2.0-or-later
# --------------------------------------------------------------------------------------------------

import io
import json
import math
from pathlib import Path
import subprocess
import sys
import unittest

import numpy as np

from blue_moon import (
    ANGSTROM_PER_BOHR,
    KB_HARTREE_PER_K,
    Coordinate,
    analyze,
    metric,
    read_multipliers,
    read_xyz,
)


def primitive(kind, atoms, **kwargs):
    return {"type": kind, "atoms": atoms, **kwargs}


def combination(*terms):
    return {
        "type": "linear_combination",
        "terms": [{"coefficient": c, "cv": cv} for c, cv in terms],
    }


DISTANCE = primitive("distance", [1, 2])
DIFFERENCE = combination((1, DISTANCE), (-1, primitive("distance", [3, 2])))
CELL = np.eye(3) * 30
POSITIONS = np.array(
    [[2.0, 0.3, 0.1], [0.0, 0.0, 0.0], [1.0, 3.0, -0.4], [0.1, 0.2, 2.4]]
)


class DerivativeTests(unittest.TestCase):
    def test_distance(self):
        result = Coordinate(DISTANCE, CELL).evaluate(POSITIONS)
        r = np.linalg.norm(POSITIONS[0] - POSITIONS[1])
        unit = (POSITIONS[0] - POSITIONS[1]) / r
        h = (np.eye(3) - np.outer(unit, unit)) / r
        self.assertAlmostEqual(result.value, r)
        np.testing.assert_allclose(result.gradient, np.r_[unit, -unit])
        np.testing.assert_allclose(
            result.hessian, np.block([[h, -h], [-h, h]]), atol=1e-15
        )
        z, g = metric(result, [12, 1])
        self.assertAlmostEqual(z, 1 / 12 + 1)
        self.assertAlmostEqual(g, 0)

    def test_distance_difference_shared_atom(self):
        x = np.array([[2.0, 0, 0], [0, 0, 0], [1.0, 3, 0]])
        coordinate = Coordinate(DIFFERENCE, CELL)
        result = coordinate.evaluate(x)
        masses = np.array([12, 1, 16])
        z, g = metric(result, masses)
        r1, r2 = np.linalg.norm(x[0]), np.linalg.norm(x[2])
        cosine = np.dot(x[0], x[2]) / r1 / r2
        z_reference = 1 / masses[0] + 1 / masses[2] + 2 / masses[1] * (1 - cosine)
        g_reference = (
            (1 - cosine**2) / (masses[1] ** 2 * z_reference**2) * (1 / r1 - 1 / r2)
        )
        self.assertAlmostEqual(z, z_reference)
        self.assertAlmostEqual(g, g_reference)
        self.assertAlmostEqual(g, 0.07221504489510591)

        # An independent directional derivative of the analytic metric (no Hessian).
        def analytic_z(y):
            a, b = y[0] - y[1], y[2] - y[1]
            c = np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b)
            return 1 / 12 + 1 / 16 + 2 * (1 - c)

        direction = result.gradient.reshape(3, 3) / masses[:, None]
        for step in (1e-3, 1e-4, 1e-5):
            numerical = (
                analytic_z(x + step * direction) - analytic_z(x - step * direction)
            ) / (4 * step * z**2)
            self.assertAlmostEqual(g, numerical, delta=1e-8)
        # Equal lengths and collinear geometries are the special G=0 cases.
        for end in ([0, 2, 0], [-3, 0, 0]):
            x[2] = end
            self.assertAlmostEqual(metric(coordinate.evaluate(x), masses)[1], 0)

    def test_point_center(self):
        coordinate = Coordinate(primitive("point_bond_center", [1, 2, 3]), CELL)
        result = coordinate.evaluate(POSITIONS)
        z, g = metric(result, [12, 1, 16])
        self.assertAlmostEqual(
            result.value,
            np.linalg.norm(POSITIONS[0] - (POSITIONS[1] + POSITIONS[2]) / 2),
        )
        self.assertAlmostEqual(z, 1 / 12 + 1 / 4 + 1 / 64)
        self.assertAlmostEqual(g, 0)

    def test_all_primitive_derivatives(self):
        cvs = [
            DISTANCE,
            DIFFERENCE,
            primitive("angle", [1, 2, 3]),
            primitive("torsion", [1, 2, 3, 4], reference=0),
            primitive("point_plane", [1, 2, 3, 4]),
            primitive("point_bond_center", [1, 2, 3]),
            combination(
                (0.7, DISTANCE), (-1.4, primitive("point_plane", [1, 2, 3, 4]))
            ),
            combination(
                (1, DISTANCE),
                (-2, primitive("distance", [3, 2])),
                (3, primitive("distance", [3, 4])),
            ),
        ]
        for cv in cvs:
            with self.subTest(cv=cv):
                coordinate = Coordinate(cv, CELL)
                result = coordinate.evaluate(POSITIONS)
                numerical_gradient = np.zeros_like(result.gradient)
                numerical_hessian = np.zeros_like(result.hessian)
                h = 2e-5
                for i, atom in enumerate(coordinate.atoms):
                    for k in range(3):
                        dx = np.zeros_like(POSITIONS)
                        dx[atom - 1, k] = h
                        plus, minus = coordinate.evaluate(
                            POSITIONS + dx
                        ), coordinate.evaluate(POSITIONS - dx)
                        numerical_gradient[3 * i + k] = (plus.value - minus.value) / (
                            2 * h
                        )
                        numerical_hessian[:, 3 * i + k] = (
                            plus.gradient - minus.gradient
                        ) / (2 * h)
                np.testing.assert_allclose(
                    result.gradient, numerical_gradient, atol=2e-9
                )
                np.testing.assert_allclose(result.hessian, numerical_hessian, atol=2e-9)
                np.testing.assert_allclose(
                    result.gradient.reshape(-1, 3).sum(axis=0), 0, atol=2e-15
                )
                np.testing.assert_allclose(result.hessian, result.hessian.T, atol=1e-15)

    def test_angle_and_signed_plane_values(self):
        angle = (
            Coordinate(primitive("angle", [1, 2, 3]), CELL).evaluate(POSITIONS).value
        )
        a, b = POSITIONS[0] - POSITIONS[1], POSITIONS[2] - POSITIONS[1]
        self.assertAlmostEqual(
            angle, math.acos(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))
        )
        normal = np.cross(a, b)
        plane = (
            Coordinate(primitive("point_plane", [1, 2, 3, 4]), CELL)
            .evaluate(POSITIONS)
            .value
        )
        self.assertAlmostEqual(
            plane,
            np.dot(POSITIONS[3] - POSITIONS[:3].mean(axis=0), normal)
            / np.linalg.norm(normal),
        )

    def test_mass_and_coordinate_scaling(self):
        result = Coordinate(DIFFERENCE, CELL).evaluate(POSITIONS)
        z, g = metric(result, [12, 1, 16])
        z2, g2 = metric(result, [24, 2, 32])
        self.assertAlmostEqual(z2, z / 2)
        self.assertAlmostEqual(g2, g)
        for factor in (2, -3):
            scaled = Coordinate(combination((factor, DIFFERENCE)), CELL).evaluate(
                POSITIONS
            )
            zs, gs = metric(scaled, [12, 1, 16])
            self.assertAlmostEqual(zs, factor**2 * z)
            self.assertAlmostEqual(gs, g / factor)

    def test_triclinic_images(self):
        cell = np.array([[20, 0, 0], [3, 24, 0], [-1, 2, 22]])
        shifted = POSITIONS.copy()
        shifted[0] += cell[1] - 2 * cell[2]
        coordinate = Coordinate(DIFFERENCE, cell)
        a, b = coordinate.evaluate(POSITIONS), coordinate.evaluate(shifted)
        self.assertAlmostEqual(a.value, b.value)
        np.testing.assert_allclose(a.gradient, b.gradient, atol=1e-14)
        np.testing.assert_allclose(a.hessian, b.hessian, atol=1e-14)
        nonperiodic = Coordinate(primitive("distance", [1, 2], pbc=False), cell)
        self.assertAlmostEqual(
            nonperiodic.evaluate(shifted).value, np.linalg.norm(shifted[0] - shifted[1])
        )

    def test_torsion_branch(self):
        phi = Coordinate(
            primitive("torsion", [1, 2, 3, 4], reference=0), CELL
        ).evaluate(POSITIONS)
        wrapped = Coordinate(
            primitive("torsion", [1, 2, 3, 4], reference=2 * math.pi), CELL
        ).evaluate(POSITIONS)
        self.assertAlmostEqual(wrapped.value - phi.value, 2 * math.pi)
        np.testing.assert_allclose(wrapped.gradient, phi.gradient)
        reversed_phi = Coordinate(
            primitive("torsion", [1, 2, 3, 4], reference=-phi.value), CELL
        )
        # A non-zero reference selects a branch, not a shift of the coordinate.
        self.assertAlmostEqual(
            math.remainder(
                reversed_phi.evaluate(POSITIONS).value - phi.value, 2 * math.pi
            ),
            0,
        )

    def test_cp2k_coordinate_values(self):
        # CP2K 2026.1 METADYN/COLVAR output for examples/coordinates.inp.
        # Values are printed with five decimal places, in internal units.
        cvs = [
            DISTANCE,
            primitive("angle", [1, 2, 3]),
            primitive("torsion", [1, 2, 3, 4], reference=0),
            primitive("point_plane", [1, 2, 3, 4]),
            primitive("point_bond_center", [1, 2, 3]),
            DIFFERENCE,
        ]
        expected = [3.82640, 1.11171, -1.43934, 4.51310, 3.67405, -2.19705]
        for cv, value in zip(cvs, expected):
            result = Coordinate(cv, CELL / ANGSTROM_PER_BOHR).evaluate(
                (POSITIONS + 10) / ANGSTROM_PER_BOHR
            )
            self.assertAlmostEqual(result.value, value, delta=5e-6)

    def test_rotational_invariance(self):
        rotation = np.array([[0.6, -0.8, 0], [0.8, 0.6, 0], [0, 0, 1]])
        coordinate = Coordinate(DIFFERENCE, CELL)
        original = coordinate.evaluate(POSITIONS)
        rotated = coordinate.evaluate(POSITIONS @ rotation)
        np.testing.assert_allclose(
            metric(original, [12, 1, 16]), metric(rotated, [12, 1, 16]), atol=1e-14
        )

    def test_inconsistent_torsion_images(self):
        coordinate = Coordinate(primitive("torsion", [1, 2, 3, 4], reference=0), CELL)
        # Adjacent bonds are shorter than half a cell, their 1-3 sum is not.
        positions = np.array([[0, 0, 0], [9, 1, 0], [18, 0, 0], [19, 1, 2]])
        with self.assertRaisesRegex(ValueError, "Inconsistent CP2K torsion images"):
            coordinate.evaluate(positions)

    def test_invalid_coordinates(self):
        invalid = [
            primitive("distance", [1, 1]),
            primitive("distance", [0, 2]),
            primitive("distance", [True, 2]),
            primitive("unknown", [1, 2]),
            {**DISTANCE, "axis": "X"},
            primitive("torsion", [1, 2, 3, 4]),
        ]
        for cv in invalid:
            with self.subTest(cv=cv), self.assertRaises(ValueError):
                Coordinate(cv, CELL)
        for masses in ([1], [1, 0], [1, float("nan")]):
            with self.assertRaises(ValueError):
                metric(Coordinate(DISTANCE, CELL).evaluate(POSITIONS), masses)
        with self.assertRaisesRegex(ValueError, "branch"):
            Coordinate(DISTANCE, CELL).evaluate([[15, 0, 0], [0, 0, 0]])
        with self.assertRaisesRegex(ValueError, "Degenerate"):
            Coordinate(primitive("angle", [1, 2, 3]), CELL).evaluate(
                [[1, 0, 0], [0, 0, 0], [2, 0, 0]]
            )
        with self.assertRaisesRegex(ValueError, "Singular mass"):
            metric(
                Coordinate(combination((1, DISTANCE), (-1, DISTANCE)), CELL).evaluate(
                    POSITIONS
                ),
                [12, 1],
            )


def config():
    return {
        "cv": DISTANCE,
        "cell_angstrom": (CELL * ANGSTROM_PER_BOHR).tolist(),
        "masses_amu": {"1": 12.0, "2": 1.0},
        "temperature_kelvin": 300.0,
        "target_au": 2.0,
        "target_tolerance_au": 1e-6,
    }


class InputTests(unittest.TestCase):
    def test_multiplier_pairs(self):
        text = "Shake  Lagrangian Multipliers: -0.1\nRattle Lagrangian Multipliers: 999\nShake  Lagrangian Multipliers: 0.2D+0\nRattle Lagrangian Multipliers: -888\n"
        self.assertEqual(list(read_multipliers(io.StringIO(text))), [-0.1, 0.2])
        for bad in (
            "Shake Lagrangian Multipliers: 1 2\n",
            "Shake Lagrangian Multipliers: 1\n",
            "Rattle Lagrangian Multipliers: 1\n",
            " 1 2\n",
            "Shake Lagrangian Multipliers: ********\n",
            "Shake Lagrangian Multipliers: NaN\n",
            "Shake Lagrangian Multipliers: 1\nShake Lagrangian Multipliers: 2\n",
        ):
            with self.subTest(text=bad), self.assertRaises(ValueError):
                list(read_multipliers(io.StringIO(bad)))

    def test_xyz(self):
        xyz = "2\n i = 0, time = 0.0\nC 0 0 0\nH 1 0 0\n2\n i = 1, time = 0.5\nC 0 0 0\nH 1 1 0\n"
        frames = list(read_xyz(io.StringIO(xyz)))
        self.assertEqual([x[0] for x in frames], [0, 1])
        self.assertAlmostEqual(frames[0][1][1, 0], 1 / ANGSTROM_PER_BOHR)
        for bad in (
            xyz + xyz,
            xyz.replace("i = 0,", "unknown"),
            xyz.replace("H 1 1 0", "O 1 1 0"),
            xyz.replace("H 1 1 0", ""),
            xyz.replace("H 1 0 0", "H nan 0 0"),
        ):
            with self.assertRaises(ValueError):
                list(read_xyz(io.StringIO(bad)))

    def test_alignment_and_units(self):
        positions = np.array([[2.0, 0, 0], [0, 0, 0]])
        frames = [(0, positions), (1, positions), (2, positions)]
        result = list(analyze(config(), iter(frames), iter([0.1, 0.2]), 1))
        self.assertEqual([x[0] for x in result], [1, 2])
        self.assertAlmostEqual(result[0][-1], -0.1)
        for steps, multipliers in (
            (frames, [0.1]),
            (frames, [0.1, 0.2, 0.3]),
            ([(2, positions)], [0.1]),
        ):
            with self.assertRaises(ValueError):
                list(analyze(config(), iter(steps), iter(multipliers), 1))
        changed = config()
        changed["target_au"] = 3
        with self.assertRaisesRegex(ValueError, "fixed target"):
            list(analyze(changed, iter(frames), iter([0.1, 0.2]), 1))
        with self.assertRaisesRegex(ValueError, "Invalid keys"):
            list(analyze({**config(), "unknown": 1}, iter(frames), iter([0.1, 0.2]), 1))

    def test_weighted_estimator(self):
        # Two geometries on the SAME nonzero distance-difference window.
        frames = [
            (1, np.array([[2, 0, 0], [0, 0, 0], [0, 3, 0]])),
            (2, np.array([[2, 0, 0], [0, 0, 0], [1.8, 2.4, 0]])),
        ]
        cfg = {
            **config(),
            "cv": DIFFERENCE,
            "masses_amu": {"1": 12, "2": 1, "3": 16},
            "target_au": -1,
        }
        samples = list(analyze(cfg, iter(frames), iter([0.1, -0.2]), 1))
        zs = np.array([1 / 12 + 1 / 16 + 2, 1 / 12 + 1 / 16 + 0.8])
        gs = np.array([1, 1 - 0.6**2]) * (1 / 2 - 1 / 3) / zs**2
        expected = np.average(
            -np.array([0.1, -0.2]) + KB_HARTREE_PER_K * 300 * gs, weights=zs**-0.5
        )
        actual = sum(s[-2] * s[-1] for s in samples) / sum(s[-2] for s in samples)
        self.assertAlmostEqual(actual, expected)
        self.assertNotAlmostEqual(actual, 0.05)

    def test_cli_example(self):
        directory = Path(__file__).parent
        command = [
            sys.executable,
            str(directory / "blue_moon.py"),
            str(directory / "examples/distance.json"),
            str(directory / "examples/trajectory.xyz"),
            str(directory / "examples/constraint.LagrangeMultLog"),
            "--first-step",
            "1",
        ]
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        summary = json.loads(result.stderr)
        self.assertEqual(summary["samples"], 5)
        self.assertAlmostEqual(
            summary["free_energy_gradient_au"],
            summary["uncorrected_minus_mean_lambda_au"],
        )
        self.assertEqual(len(result.stdout.splitlines()), 6)
        failed = subprocess.run(
            command + ["--discard", "5"], capture_output=True, text=True
        )
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("No production samples", failed.stderr)


if __name__ == "__main__":
    unittest.main()
