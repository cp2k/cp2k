#!/usr/bin/env python3

import math
import subprocess
import tempfile
import unittest
from pathlib import Path

TOOL_DIR = Path(__file__).resolve().parents[1]
TRAJANA = TOOL_DIR / "trajana.x"
TRAJCONVERT = TOOL_DIR / "trajconvert.x"
CELL = "10 0 0 0 10 0 0 0 10"


def run(*arguments):
    completed = subprocess.run(arguments, cwd=TOOL_DIR, text=True, capture_output=True)
    if completed.returncode != 0:
        raise AssertionError(
            f"Command failed ({completed.returncode}): {' '.join(arguments)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )


def data_lines(path):
    return [
        [float(field) for field in line.split()]
        for line in Path(path).read_text(encoding="utf8").splitlines()
        if line and not line.startswith("#")
    ]


def xyz_frames(path):
    lines = Path(path).read_text(encoding="utf8").splitlines()
    frames = []
    position = 0
    while position < len(lines):
        natoms = int(lines[position])
        position += 2
        atoms = []
        for _ in range(natoms):
            fields = lines[position].split()
            atoms.append((fields[0], tuple(float(value) for value in fields[1:4])))
            position += 1
        frames.append(atoms)
    return frames


class TrajectoryToolTests(unittest.TestCase):
    def test_wrap_fractional_and_unwrap(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "crossing.xyz"
            source.write_text(
                "1\nframe 1\nX 9.5 0 0\n" "1\nframe 2\nX 0.5 0 0\n",
                encoding="utf8",
            )
            unwrapped = root / "unwrapped.xyz"
            run(
                str(TRAJCONVERT),
                "unwrap",
                "--input",
                str(source),
                "--output",
                str(unwrapped),
                "--cell",
                CELL,
            )
            frames = xyz_frames(unwrapped)
            self.assertAlmostEqual(frames[0][0][1][0], 9.5)
            self.assertAlmostEqual(frames[1][0][1][0], 10.5)

            wrapped_source = root / "outside.xyz"
            wrapped_source.write_text("1\noutside\nX -1 11 5\n", encoding="utf8")
            wrapped = root / "wrapped.xyz"
            fractional = root / "fractional.xyz"
            run(
                str(TRAJCONVERT),
                "wrap",
                "--input",
                str(wrapped_source),
                "--output",
                str(wrapped),
                "--cell",
                CELL,
            )
            run(
                str(TRAJCONVERT),
                "fractional",
                "--input",
                str(wrapped_source),
                "--output",
                str(fractional),
                "--cell",
                CELL,
            )
            self.assertEqual(
                tuple(round(x, 12) for x in xyz_frames(wrapped)[0][0][1]),
                (9.0, 1.0, 5.0),
            )
            self.assertEqual(
                tuple(round(x, 12) for x in xyz_frames(fractional)[0][0][1]),
                (0.9, 0.1, 0.5),
            )

            changing = root / "changing.xyz"
            changing.write_text(
                "1\nframe 1\nX 11 0 0\n1\nframe 2\nX 11 0 0\n",
                encoding="utf8",
            )
            cell_file = root / "changing.cell"
            cell_file.write_text(
                "# step time ax ay az bx by bz cx cy cz volume\n"
                "0 0 10 0 0 0 10 0 0 0 10 1000\n"
                "1 1 12 0 0 0 12 0 0 0 12 1728\n",
                encoding="utf8",
            )
            changing_output = root / "changing-wrapped.xyz"
            run(
                str(TRAJCONVERT),
                "wrap",
                "--input",
                str(changing),
                "--output",
                str(changing_output),
                "--cell-file",
                str(cell_file),
            )
            changing_frames = xyz_frames(changing_output)
            self.assertAlmostEqual(changing_frames[0][0][1][0], 1.0)
            self.assertAlmostEqual(changing_frames[1][0][1][0], 11.0)

    def test_group_center_reconstructs_across_boundary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "water.xyz"
            source.write_text(
                "3\nwater\nO 9.8 0 0\nH 0.2 0 0\nH 9.6 0 0\n", encoding="utf8"
            )
            groups = root / "groups.dat"
            groups.write_text("W 1 2 3\n", encoding="utf8")
            masses = root / "masses.dat"
            masses.write_text("O 15.9994\nH 1.0079\n", encoding="utf8")
            output = root / "center.xyz"
            run(
                str(TRAJCONVERT),
                "center",
                "--input",
                str(source),
                "--output",
                str(output),
                "--groups",
                str(groups),
                "--mode",
                "mass",
                "--mass-file",
                str(masses),
                "--cell",
                CELL,
                "--wrap-output",
            )
            center = xyz_frames(output)[0][0][1][0]
            expected = (15.9994 * 9.8 + 1.0079 * 10.2 + 1.0079 * 9.6) / 18.0152
            self.assertAlmostEqual(center, expected, places=12)

    def test_geometry_and_rdf(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "pair.xyz"
            source.write_text(
                "4\ngeometry\nA 1 0 0\nB 0 0 0\nC 0 1 0\nD 0 1 1\n",
                encoding="utf8",
            )
            actions = root / "actions.in"
            actions.write_text(
                "DISTANCE 1 2\nANGLE 1 2 3\nTORSION 1 2 3 4\n", encoding="utf8"
            )
            geometry = root / "geometry.dat"
            run(
                str(TRAJANA),
                "geometry",
                "--input",
                str(source),
                "--output",
                str(geometry),
                "--actions",
                str(actions),
            )
            geometry_values = data_lines(geometry)[0]
            self.assertAlmostEqual(geometry_values[1], 1.0)
            self.assertAlmostEqual(geometry_values[2], 90.0)
            self.assertAlmostEqual(geometry_values[3], 90.0)

            skewed_source = root / "skewed.xyz"
            skewed_source.write_text(
                "2\nskewed cell\nA 0 0 0\nB 9.31 0.49 0\n", encoding="utf8"
            )
            skewed_actions = root / "skewed.in"
            skewed_actions.write_text("DISTANCE 1 2\n", encoding="utf8")
            skewed_output = root / "skewed.dat"
            run(
                str(TRAJANA),
                "geometry",
                "--input",
                str(skewed_source),
                "--output",
                str(skewed_output),
                "--actions",
                str(skewed_actions),
                "--cell",
                "10 0 0 9 1 0 0 0 10",
            )
            self.assertAlmostEqual(
                data_lines(skewed_output)[0][1], math.sqrt(0.31**2 + 0.51**2)
            )

            rdf = root / "rdf.dat"
            run(
                str(TRAJANA),
                "rdf",
                "--input",
                str(source),
                "--output",
                str(rdf),
                "--cell",
                CELL,
                "--select-a",
                "name:A",
                "--select-b",
                "name:B",
                "--bins",
                "5",
                "--rmax",
                "5",
            )
            values = data_lines(rdf)
            expected_g = 1000.0 / ((4.0 * math.pi / 3.0) * (2.0**3 - 1.0**3))
            self.assertAlmostEqual(values[1][0], 1.5)
            self.assertAlmostEqual(values[1][1], expected_g, places=12)
            self.assertAlmostEqual(values[1][2], 1.0)

    def test_vacf_uses_all_time_origins(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "velocity.xyz"
            velocities = [1.0, 2.0, 4.0, 8.0]
            source.write_text(
                "".join(
                    f"1\nframe {frame}\nH {velocity} 0 0\n"
                    for frame, velocity in enumerate(velocities)
                ),
                encoding="utf8",
            )
            output = root / "vacf.dat"
            spectrum = root / "spectrum.dat"
            run(
                str(TRAJANA),
                "vacf",
                "--input",
                str(source),
                "--output",
                str(output),
                "--spectrum",
                str(spectrum),
                "--dt-fs",
                "0.5",
            )
            values = data_lines(output)
            self.assertEqual(len(values), 4)
            expected = [21.25, 14.0, 10.0, 8.0]
            for row, reference in zip(values, expected):
                self.assertAlmostEqual(row[1], reference, places=12)
                self.assertAlmostEqual(row[2], reference / expected[0], places=12)
            spectrum_values = data_lines(spectrum)
            self.assertGreater(spectrum_values[0][1], 0.0)
            self.assertGreater(sum(row[1] for row in spectrum_values[1:]), 0.0)

    def test_hbond_and_dynamic_structure_factor(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            water = root / "water.xyz"
            water.write_text(
                "6\nwater dimer\n"
                "O 0 0 0\nH 1 0 0\nH 0 1 0\n"
                "O 2.8 0 0\nH 2.8 1 0\nH 2.8 0 1\n",
                encoding="utf8",
            )
            groups = root / "groups.dat"
            groups.write_text("W1 1 2 3\nW2 4 5 6\n", encoding="utf8")
            output = root / "hbond.dat"
            summary = root / "hbond-summary.dat"
            run(
                str(TRAJANA),
                "hbond",
                "--input",
                str(water),
                "--output",
                str(output),
                "--summary",
                str(summary),
                "--groups",
                str(groups),
                "--cell",
                CELL,
            )
            hbond = data_lines(output)[0]
            self.assertEqual(int(hbond[1]), 1)
            self.assertAlmostEqual(hbond[2], 0.5)
            self.assertAlmostEqual(hbond[3], 0.75)
            populations = {
                (int(row[0]), int(row[1])): row[2] for row in data_lines(summary)
            }
            self.assertAlmostEqual(populations[(0, 1)], 0.5)
            self.assertAlmostEqual(populations[(1, 0)], 0.5)

            stationary = root / "stationary.xyz"
            stationary.write_text(
                "".join(f"2\nframe {frame}\nX 0 0 0\nX 1 0 0\n" for frame in range(4)),
                encoding="utf8",
            )
            q_file = root / "q.dat"
            q_file.write_text(f"{math.pi / 3.0:.16f} 0 0\n", encoding="utf8")
            intermediate = root / "intermediate.dat"
            dynamic = root / "dynamic.dat"
            run(
                str(TRAJANA),
                "dsf",
                "--input",
                str(stationary),
                "--output",
                str(intermediate),
                "--spectrum",
                str(dynamic),
                "--q-file",
                str(q_file),
                "--dt-fs",
                "1",
            )
            for row in data_lines(intermediate):
                # |1 + exp(i*pi/3)|^2 / N = 3/2. This includes the
                # coherent cross term that the original dynStruct.f90 missed.
                self.assertAlmostEqual(row[1], 1.5, places=12)
                self.assertAlmostEqual(row[2], 0.0, places=12)
                self.assertAlmostEqual(row[3], 1.0, places=12)
            dynamic_values = data_lines(dynamic)
            self.assertGreater(dynamic_values[0][1], 0.0)
            self.assertAlmostEqual(
                sum(row[1] for row in dynamic_values[1:]), 0.0, places=12
            )


if __name__ == "__main__":
    unittest.main()
