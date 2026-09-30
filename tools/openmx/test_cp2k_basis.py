# SPDX-License-Identifier: GPL-2.0-or-later
"""Optional CP2K integral/SCF test; set CP2K_EXE to enable it."""

import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

import numpy as np


def printed_matrix(output, label, nao):
    """Read only the small AO matrices printed by the supplied smoke input."""
    start = output.rindex("\n " + label + "\n")
    matrix = np.full((nao, nao), np.nan)
    columns = None
    for line in output[start:].splitlines():
        if re.fullmatch(r"\s+(?:[0-9]+\s*)+", line):
            columns = np.array([int(x) - 1 for x in line.split()])
        elif columns is not None and re.match(
            r"\s*\d+\s+\d+\s+[A-Za-z]+\s+\S+\s+", line
        ):
            fields = line.split()
            row = int(fields[0]) - 1
            if not 0 <= row < nao or np.any(columns < 0) or np.any(columns >= nao):
                raise ValueError("Unexpected AO indices in smoke output")
            matrix[row, columns] = [float(x) for x in fields[4:]]
            if np.isfinite(matrix).all():
                return matrix
    raise ValueError(f"Incomplete {label}")


@unittest.skipUnless(os.environ.get("CP2K_EXE"), "Set CP2K_EXE for the SCF smoke test")
class Cp2kBasisTests(unittest.TestCase):
    def test_basis_integrals_and_electron_count(self):
        tool_dir = Path(__file__).resolve().parent
        data_dir = tool_dir.parents[1] / "data"
        cp2k = str(Path(os.environ["CP2K_EXE"]).resolve())
        text = (tool_dir / "examples" / "H2O_basis_smoke.inp").read_text()
        text = text.replace("../../../data/", str(data_dir) + "/")
        with tempfile.TemporaryDirectory(prefix="openmx-basis-scf-") as temporary:
            root = Path(temporary)
            (root / "h2o.inp").write_text(text)
            run = subprocess.run(
                [cp2k, "-i", "h2o.inp", "-o", "h2o.out"],
                cwd=root,
                env=dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1"),
                capture_output=True,
                text=True,
                timeout=180,
            )
            output = (root / "h2o.out").read_text()
            self.assertEqual(run.returncode, 0, run.stderr + output[-6000:])
            self.assertIn("SCF run converged", output)
            self.assertIn("PROGRAM ENDED", output)
            self.assertNotIn("WARNING", output)
            self.assertRegex(output, r"Number of electrons:\s+8\b")
            self.assertRegex(output, r"Number of orbital functions:\s+14\b")
            s = printed_matrix(output, "OVERLAP MATRIX", 14)
            t = printed_matrix(output, "KINETIC ENERGY MATRIX", 14)
            p = printed_matrix(output, "DENSITY MATRIX", 14)
            report = json.loads((data_dir / "BASIS_OMX.json").read_text())
            by_name = {b["cp2k_name"]: b for b in report["basis_sets"]}
            selections = [
                "OMX-FIT-O7.0-s1p1",
                "OMX-FIT-H7.0-s2p1",
                "OMX-FIT-H7.0-s2p1",
            ]
            expected_t = [
                kinetic
                for name in selections
                for channel in by_name[name]["radial_fits"]
                for kinetic in channel["kinetic_fit_hartree"]
                for _ in range(2 * channel["l"] + 1)
            ]
            np.testing.assert_allclose(s, s.T, atol=1e-10, rtol=0)
            self.assertGreater(np.linalg.eigvalsh(s)[0], 1e-6)
            np.testing.assert_allclose(np.diag(s), 1, atol=1e-8, rtol=0)
            np.testing.assert_allclose(np.diag(t), expected_t, atol=1e-8, rtol=0)
            self.assertLess(abs(np.trace(p @ s) - 8), 1e-7)


if __name__ == "__main__":
    unittest.main()
