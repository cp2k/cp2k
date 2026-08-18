#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
#   CP2K: A general program to perform molecular dynamics simulations
#   Copyright 2000-2026 CP2K developers group <https://cp2k.org>
#
#   SPDX-License-Identifier: GPL-2.0-or-later
# --------------------------------------------------------------------------------------------------

import math
import pathlib
import subprocess
import sys
import tempfile
import unittest

EXECUTABLE = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else None


def write_xyz(path, frames):
    with path.open("w", encoding="utf8") as handle:
        for header, atoms in frames:
            handle.write(f"{len(atoms)}\n")
            if isinstance(header, int):
                handle.write(f"i = {header}, time = 0.0, E = 0.0\n")
            else:
                handle.write(header + "\n")
            for label, coordinates in atoms:
                handle.write(
                    f"{label:2s} {coordinates[0]:16.9f} {coordinates[1]:16.9f} "
                    f"{coordinates[2]:16.9f}\n"
                )


def read_cube_values(path):
    with path.open(encoding="utf8") as handle:
        handle.readline()
        handle.readline()
        natoms = int(handle.readline().split()[0])
        shape = [int(handle.readline().split()[0]) for _ in range(3)]
        for _ in range(abs(natoms)):
            handle.readline()
        values = [float(value) for line in handle for value in line.split()]
    assert len(values) == math.prod(shape)
    return shape, values


def rotate_z(coordinates):
    x, y, z = coordinates
    return [-y, x, z]


class SpatialDistributionTest(unittest.TestCase):
    def run_tool(self, directory, input_text):
        input_path = directory / "sdf.in"
        input_path.write_text(input_text, encoding="utf8")
        result = subprocess.run(
            [str(EXECUTABLE), str(input_path)],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, msg=result.stdout + result.stderr)
        return result

    def test_water_frame_rotation_and_cell_step_matching(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = pathlib.Path(temporary_directory)
            center = [5.0, 5.0, 5.0]
            origins = [[3.0, 5.0, 5.0], [7.0, 5.0, 5.0]]
            absolute = [
                ("O", origins[0]),
                ("O", origins[1]),
                ("H", [origins[0][0] - 0.75, origins[0][1], origins[0][2] + 0.60]),
                ("H", [origins[0][0] + 0.75, origins[0][1], origins[0][2] + 0.60]),
                ("H", [origins[1][0] - 0.75, origins[1][1], origins[1][2] + 0.60]),
                ("H", [origins[1][0] + 0.75, origins[1][1], origins[1][2] + 0.60]),
                (
                    "Ne",
                    [origins[0][0] + 1.20, origins[0][1] - 0.70, origins[0][2] + 2.10],
                ),
            ]

            frames = []
            for step, transform in [(0, lambda vector: vector), (10, rotate_z)]:
                atoms = []
                for label, coordinates in absolute:
                    relative = [coordinates[i] - center[i] for i in range(3)]
                    rotated = transform(relative)
                    atoms.append((label, [center[i] + rotated[i] for i in range(3)]))
                frames.append((step, atoms))
            write_xyz(directory / "water.xyz", frames)
            (directory / "water.cell").write_text(
                "# Step Time Ax Ay Az Bx By Bz Cx Cy Cz Volume\n"
                "0 0.0 10 0 0 0 10 0 0 0 10 1000\n"
                "10 1.0 10 0 0 0 10 0 0 0 10 1000\n",
                encoding="utf8",
            )

            self.run_tool(
                directory,
                """TRAJECTORY water.xyz
CELL_FILE water.cell
REFERENCE_PATTERN 1 3 4 1 2 2 2
FRAME_MODE BISECTOR
TARGET neon Ne
BOUNDS -4 4 -4 4 -4 4
GRID 8 8 8
OUTPUT_PREFIX water
""",
            )
            shape, values = read_cube_values(directory / "water-neon.cube")
            self.assertEqual(shape, [8, 8, 8])
            nonzero = [value for value in values if abs(value) > 1.0e-12]
            self.assertEqual(len(nonzero), 2)
            for value in nonzero:
                self.assertAlmostEqual(value, 500.0, places=8)

            reference_lines = (
                (directory / "water-reference.xyz")
                .read_text(encoding="utf8")
                .splitlines()
            )
            arm1 = [float(value) for value in reference_lines[3].split()[1:]]
            arm2 = [float(value) for value in reference_lines[4].split()[1:]]
            for actual, expected in zip(arm1, [-0.75, 0.0, 0.60]):
                self.assertAlmostEqual(actual, expected, places=8)
            for actual, expected in zip(arm2, [0.75, 0.0, 0.60]):
                self.assertAlmostEqual(actual, expected, places=8)

    def test_general_triad_with_triclinic_cell(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = pathlib.Path(temporary_directory)
            output_directory = directory / "results with spaces"
            output_directory.mkdir()
            origin = [11.7, 7.2, 5.4]
            atoms = [
                ("C", origin),
                ("N", [origin[0], origin[1], origin[2] + 1.0]),
                ("O", [origin[0] + 1.0, origin[1], origin[2]]),
                ("He", [origin[0] - 1.2 - 8.0, origin[1] + 0.4, origin[2] - 0.6]),
            ]
            extxyz_header = (
                'Lattice="8 0 0 3 7 0 2 1 6" ' "Properties=species:S:1:pos:R:3 Step=7"
            )
            write_xyz(directory / "general trajectory.xyz", [(extxyz_header, atoms)])

            self.run_tool(
                directory,
                """TRAJECTORY "general trajectory.xyz"
REFERENCE 1 2 3
FRAME_MODE TRIAD
TARGET remaining *
BOUNDS -2 2 -2 2 -2 2
GRID 4 4 4
OUTPUT_PREFIX "results with spaces/general"
""",
            )
            shape, values = read_cube_values(
                output_directory / "general-remaining.cube"
            )
            self.assertEqual(shape, [4, 4, 4])
            nonzero = [value for value in values if abs(value) > 1.0e-12]
            self.assertEqual(len(nonzero), 1)
            self.assertAlmostEqual(nonzero[0], 84.0, places=8)


if __name__ == "__main__":
    if EXECUTABLE is None:
        raise SystemExit("usage: test_spatial_distribution.py EXECUTABLE")
    unittest.main(argv=[sys.argv[0]])
