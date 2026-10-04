#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
#   CP2K: A general program to perform molecular dynamics simulations
#   Copyright 2000-2026 CP2K developers group <https://cp2k.org>
#
#   SPDX-License-Identifier: GPL-2.0-or-later
# --------------------------------------------------------------------------------------------------

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from rixs_postprocess import get_sigma_tensor


class RIXSPostprocessTests(unittest.TestCase):
    def test_core_state_weights_accumulate_independently(self):
        core_energies = np.array([1.0, 3.0])
        emission_energies = np.array([[2.0], [2.0]])
        dipoles = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        transition_dipoles = dipoles[:, None, :]
        tensor = get_sigma_tensor(
            0, 1.0, core_energies, emission_energies, dipoles, transition_dipoles, 1.0
        )
        self.assertAlmostEqual(tensor[0, 0, 0, 0], 4.0 + 36.0 / 5.0)

    def test_normalized_spectrum_uses_each_core_resonance(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "sample.rixs"
            source.write_text(
                "Excitation from ground-state\n"
                "5 1 0 0 1\n"
                "7 0 1 1 2\n"
                "Emission from core-excited state\n"
                "5 1 1 1 0 2\n"
                "5 2 0 1 0 1\n"
                "7 1 1 0 1 2\n"
                "7 2 1 1 1 3\n",
                encoding="utf8",
            )
            subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).with_name("rixs_postprocess.py")),
                    "--filename",
                    source.name,
                    "--v_states",
                    "2",
                    "--gamma",
                    "1",
                ],
                cwd=directory,
                check=True,
                capture_output=True,
                text=True,
            )
            actual = np.loadtxt(Path(directory) / "sample_rixs.dat")
            incident = np.linspace(4.8, 5.2, 50)[:, None]
            loss = np.linspace(0.0, 2.0, 50)[None, :]
            expected = np.zeros((50, 50))
            ground_dipoles = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]])
            emission_dipoles = np.array(
                [[[1.0, 1.0, 0.0], [0.0, 1.0, 0.0]], [[1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]]
            )
            cosine_squared = np.cos(np.radians(60.0)) ** 2
            a = 3.0 + cosine_squared
            b = 0.5 * (1.0 - 3.0 * cosine_squared)
            for final_index, final in enumerate((1.0, 2.0)):
                for core_index, core in enumerate((5.0, 7.0)):
                    ground = ground_dipoles[core_index]
                    transition = emission_dipoles[core_index, final_index]
                    resonance = (
                        (core - final) ** 2 * core**2 / ((incident - core) ** 2 + 1.0)
                    )
                    emission = 1.0 / ((final - loss) ** 2 + 1.0)
                    angular = (
                        a * np.dot(ground, ground) * np.dot(transition, transition)
                        + 2.0 * b * np.dot(ground, transition) ** 2
                    )
                    expected += resonance * emission * angular
            expected /= expected.max()
            np.testing.assert_allclose(actual, expected, atol=5.1e-7, rtol=0)


if __name__ == "__main__":
    unittest.main()
