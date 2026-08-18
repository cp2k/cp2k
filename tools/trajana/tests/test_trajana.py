#!/usr/bin/env python3

import math
import subprocess
import tempfile
import unittest
from pathlib import Path

TOOL_DIR = Path(__file__).resolve().parents[1]
TRAJANA = TOOL_DIR / "trajana.x"
TRAJCONVERT = TOOL_DIR / "trajconvert.x"
TWOPT = TOOL_DIR / "twopt.x"
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


def named_data_lines(path):
    return {
        fields[0]: [float(value) for value in fields[1:]]
        for line in Path(path).read_text(encoding="utf8").splitlines()
        if line and not line.startswith("#")
        for fields in [line.split()]
    }


class TrajectoryToolTests(unittest.TestCase):
    def test_twopt_atomic_fluidicity_and_sum_rule(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            velocities = root / "argon-vel.xyz"
            frames = 32
            velocities.write_text(
                "".join(
                    "2\nframe {frame}\nAr {speed:.16f} 0 0\nAr {opposite:.16f} 0 0\n".format(
                        frame=frame,
                        speed=0.5 + math.sin(2.0 * math.pi * frame / frames),
                        opposite=-0.5 - math.sin(2.0 * math.pi * frame / frames),
                    )
                    for frame in range(frames)
                ),
                encoding="utf8",
            )
            masses = root / "masses.dat"
            masses.write_text("Ar 39.948\n", encoding="utf8")
            thermo = root / "thermo.dat"
            spectrum = root / "spectrum.dat"
            run(
                str(TWOPT),
                "--velocity",
                str(velocities),
                "--velocity-unit",
                "angstrom/ps",
                "--mass-file",
                str(masses),
                "--dt-fs",
                "1000",
                "--temperature",
                "300",
                "--volume",
                "1000",
                "--output",
                str(thermo),
                "--spectrum",
                str(spectrum),
            )
            channels = named_data_lines(thermo)
            translation = channels["translation"]
            self.assertAlmostEqual(translation[0], 6.0, places=11)
            self.assertGreater(translation[3], 0.0)
            self.assertGreater(translation[4], 0.0, translation)
            self.assertLessEqual(translation[4], 1.0)
            expected_k = (
                (translation[2] / 2.99792458e10 * 1.0e-2)
                / 2.0
                * math.sqrt(
                    math.pi * 6.02214076e23 * 1.380649e-23 * 300.0 / (39.948e-3)
                )
                * 2.0
                / 9.0
                * (2.0 / 1000.0) ** (1.0 / 3.0)
                * 1.0e10
                * (6.0 / math.pi) ** (2.0 / 3.0)
            )
            self.assertAlmostEqual(translation[3], expected_k, places=12)
            k_value = translation[3]
            fluidicity = translation[4]
            residual = (
                2.0 * k_value ** (-4.5) * fluidicity**7.5
                - 6.0 * k_value ** (-3.0) * fluidicity**5.0
                - k_value ** (-1.5) * fluidicity**3.5
                + 6.0 * k_value ** (-1.5) * fluidicity**2.5
                + 2.0 * fluidicity
                - 2.0
            )
            self.assertAlmostEqual(residual, 0.0, places=8)
            self.assertTrue(all(math.isfinite(value) for value in channels["total"]))

            values = data_lines(spectrum)
            frequency_step = values[1][0] - values[0][0]
            integrated = frequency_step * (
                0.5 * (values[0][1] + values[-1][1])
                + sum(row[1] for row in values[1:-1])
            )
            self.assertAlmostEqual(integrated, 6.0, places=10)

    def test_twopt_molecular_decomposition(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            positions = root / "water-pos.xyz"
            velocities = root / "water-vel.xyz"
            frame_count = 32
            position_text = []
            velocity_text = []
            coordinates = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
            masses_value = (15.9994, 1.0079, 1.0079)
            center = tuple(
                sum(
                    mass * coordinate[axis]
                    for mass, coordinate in zip(masses_value, coordinates)
                )
                / sum(masses_value)
                for axis in range(3)
            )
            for frame in range(frame_count):
                phase = 2.0 * math.pi * frame / frame_count
                translation = (math.cos(phase), 0.5 * math.sin(phase), 0.2)
                omega = (
                    0.1 * math.sin(phase),
                    0.2 * math.cos(phase),
                    0.05 * math.sin(2.0 * phase),
                )
                internal_h1 = (0.3 * math.cos(3.0 * phase), 0.0, 0.0)
                internal_h2 = (0.0, -0.3 * math.cos(3.0 * phase), 0.0)
                internal_o = tuple(
                    -(
                        masses_value[1] * internal_h1[axis]
                        + masses_value[2] * internal_h2[axis]
                    )
                    / masses_value[0]
                    for axis in range(3)
                )
                internal = (internal_o, internal_h1, internal_h2)
                atom_velocities = []
                for coordinate, extra in zip(coordinates, internal):
                    relative = tuple(
                        coordinate[axis] - center[axis] for axis in range(3)
                    )
                    rotation = (
                        omega[1] * relative[2] - omega[2] * relative[1],
                        omega[2] * relative[0] - omega[0] * relative[2],
                        omega[0] * relative[1] - omega[1] * relative[0],
                    )
                    atom_velocities.append(
                        tuple(
                            translation[axis] + rotation[axis] + extra[axis]
                            for axis in range(3)
                        )
                    )
                position_text.append(
                    '3\nLattice="10 0 0 0 10 0 0 0 10"\n'
                    + "".join(
                        f"{label} {coordinate[0]} {coordinate[1]} {coordinate[2]}\n"
                        for label, coordinate in zip(("O", "H", "H"), coordinates)
                    )
                )
                velocity_text.append(
                    f"3\nframe {frame}\n"
                    + "".join(
                        f"{label} {velocity[0]:.16f} {velocity[1]:.16f} {velocity[2]:.16f}\n"
                        for label, velocity in zip(("O", "H", "H"), atom_velocities)
                    )
                )
            positions.write_text("".join(position_text), encoding="utf8")
            velocities.write_text("".join(velocity_text), encoding="utf8")
            groups = root / "groups.dat"
            groups.write_text("W 1 2 3\n", encoding="utf8")
            masses = root / "masses.dat"
            masses.write_text("O 15.9994\nH 1.0079\n", encoding="utf8")
            thermo = root / "thermo.dat"
            spectrum = root / "spectrum.dat"
            run(
                str(TWOPT),
                "--velocity",
                str(velocities),
                "--position",
                str(positions),
                "--velocity-unit",
                "angstrom/ps",
                "--mass-file",
                str(masses),
                "--groups",
                str(groups),
                "--dt-fs",
                "2",
                "--temperature",
                "300",
                "--rotational-symmetry",
                "2",
                "--keep-system-drift",
                "--output",
                str(thermo),
                "--spectrum",
                str(spectrum),
            )
            channels = named_data_lines(thermo)
            self.assertAlmostEqual(channels["translation"][0], 3.0, places=10)
            self.assertAlmostEqual(channels["rotation"][0], 3.0, places=10)
            self.assertAlmostEqual(channels["vibration"][0], 3.0, places=10)
            self.assertAlmostEqual(channels["total"][0], 9.0, places=10)

            values = data_lines(spectrum)
            frequency_step = values[1][0] - values[0][0]
            for column, reference in ((2, 3.0), (5, 3.0), (8, 3.0)):
                integrated = frequency_step * (
                    0.5 * (values[0][column] + values[-1][column])
                    + sum(row[column] for row in values[1:-1])
                )
                self.assertAlmostEqual(integrated, reference, places=10)

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
            self.assertGreater(spectrum_values[0][3], 0.0)
            self.assertGreater(sum(row[3] for row in spectrum_values[1:]), 0.0)
            self.assertAlmostEqual(spectrum_values[1][2], 4.135667696 * 500.0)

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
            self.assertGreater(dynamic_values[0][3], 0.0)
            self.assertAlmostEqual(
                sum(row[3] for row in dynamic_values[1:]), 0.0, places=12
            )

    def test_self_scattering_and_collective_currents(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            positions = root / "positions.xyz"
            velocities = root / "velocities.xyz"
            positions.write_text(
                "".join(f"1\nframe {frame}\nH {frame} 0 0\n" for frame in range(8)),
                encoding="utf8",
            )
            velocities.write_text(
                "".join(f"1\nframe {frame}\nH 1 2 3\n" for frame in range(8)),
                encoding="utf8",
            )
            q_file = root / "q.dat"
            q_file.write_text(f"{math.pi / 2.0:.16f} 0 0\n", encoding="utf8")
            coherent = root / "coherent.dat"
            incoherent = root / "incoherent.dat"
            incoherent_spectrum = root / "incoherent-spectrum.dat"
            currents = root / "currents.dat"
            current_spectrum = root / "current-spectrum.dat"
            peaks = root / "peaks.dat"
            summary = root / "summary.dat"
            run(
                str(TRAJANA),
                "dsf",
                "--input",
                str(positions),
                "--velocity",
                str(velocities),
                "--q-file",
                str(q_file),
                "--dt-fs",
                "1",
                "--output",
                str(coherent),
                "--self-output",
                str(incoherent),
                "--self-spectrum",
                str(incoherent_spectrum),
                "--current",
                str(currents),
                "--current-spectrum",
                str(current_spectrum),
                "--summary",
                str(summary),
                "--mass-density",
                "1",
                "--peaks",
                str(peaks),
            )

            self_values = data_lines(incoherent)
            self.assertAlmostEqual(self_values[0][1], 1.0, places=12)
            self.assertAlmostEqual(self_values[1][1], 0.0, places=12)
            self.assertAlmostEqual(abs(self_values[1][2]), 1.0, places=12)
            self.assertAlmostEqual(self_values[2][1], -1.0, places=12)
            incoherent_spectrum_values = data_lines(incoherent_spectrum)
            self.assertAlmostEqual(
                incoherent_spectrum_values[2][3],
                incoherent_spectrum_values[2][4],
                places=12,
            )

            current_values = data_lines(currents)
            self.assertAlmostEqual(current_values[0][1], 1.0, places=12)
            self.assertAlmostEqual(current_values[0][4], 6.5, places=12)
            spectrum_values = data_lines(current_spectrum)
            peak_bin = spectrum_values[2]
            self.assertGreater(peak_bin[3], 0.0)
            self.assertAlmostEqual(peak_bin[3], peak_bin[5], places=12)

            summary_values = data_lines(summary)[0]
            self.assertAlmostEqual(summary_values[0], math.pi / 2.0, places=12)
            self.assertAlmostEqual(summary_values[1], 1.0, places=12)
            self.assertAlmostEqual(summary_values[2], 1.0, places=12)
            self.assertAlmostEqual(summary_values[3], 6.5, places=12)
            self.assertAlmostEqual(summary_values[4], (math.pi / 2.0) ** 2, places=12)
            self.assertAlmostEqual(summary_values[5], (math.pi / 2.0) ** 2, places=12)
            self.assertAlmostEqual(summary_values[6], (math.pi / 2.0) ** 2, places=12)
            self.assertAlmostEqual(summary_values[7], 1.0, places=12)
            self.assertAlmostEqual(summary_values[8], 1.0, places=12)
            self.assertAlmostEqual(summary_values[9], 10000.0, places=8)
            self.assertAlmostEqual(summary_values[10], 10000.0, places=8)

            peak_rows = [
                line.split()
                for line in peaks.read_text(encoding="utf8").splitlines()
                if line and not line.startswith("#")
            ]
            peak_channels = {row[0] for row in peak_rows}
            self.assertEqual(
                peak_channels,
                {"coherent", "incoherent", "longitudinal", "transverse"},
            )
            for row in peak_rows:
                self.assertAlmostEqual(float(row[5]), 100.0, places=12)


if __name__ == "__main__":
    unittest.main()
