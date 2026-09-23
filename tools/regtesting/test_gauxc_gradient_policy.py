#!/usr/bin/env python3
"""Test GauXC build capabilities and compiled gradient selection without SCFs.

Run with python3 tools/regtesting/test_gauxc_gradient_policy.py.
Requires CMake and GNU Fortran (FC may select its executable).
"""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class GradientPolicyTest(unittest.TestCase):
    def test_dependency_patch_is_repeatable_and_rejects_drift(self):
        cmake = shutil.which("cmake")
        if not cmake or not shutil.which("patch"):
            self.skipTest("CMake and patch are required")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "dependency.txt"
            patch = root / "dependency.patch"
            patch.write_text(
                "--- a/dependency.txt\n+++ b/dependency.txt\n"
                "@@ -1 +1 @@\n-original\n+corrected\n"
            )
            command = [
                cmake,
                f"-DPATCH_FILE={patch}",
                "-P",
                str(
                    ROOT
                    / "tools/toolchain/scripts/stage6/gauxc-apply-dependency-patch.cmake"
                ),
            ]
            target.write_text("original\n")
            for _ in range(2):
                result = subprocess.run(
                    command, cwd=root, capture_output=True, text=True
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(target.read_text(), "corrected\n")
            target.write_text("incompatible upstream change\n")
            result = subprocess.run(command, cwd=root, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(target.read_text(), "incompatible upstream change\n")

    def test_compiled_policy(self):
        compiler = shutil.which(os.environ.get("FC", "gfortran"))
        if not compiler:
            self.skipTest("GNU Fortran not available")
        for defines in ([], ["__GAUXC"], ["__GAUXC", "__GAUXC_ONEDFT_GRADIENT_FIX"]):
            with self.subTest(defines=defines), tempfile.TemporaryDirectory() as tmp:
                executable = Path(tmp) / "gradient-policy-test"
                subprocess.run(
                    [compiler, "-cpp", "-ffree-form", "-fcheck=all", "-Wall", "-Werror"]
                    + ["-D" + define for define in defines]
                    + [
                        str(ROOT / "src/xc/xc_gauxc_gradient_policy.F"),
                        str(ROOT / "tools/regtesting/gauxc_gradient_policy_test.F"),
                        "-o",
                        str(executable),
                    ],
                    cwd=tmp,
                    check=True,
                    capture_output=True,
                    text=True,
                )
                result = subprocess.run(
                    [str(executable)], check=True, capture_output=True, text=True
                )
                self.assertIn("Gradient policy checks passed", result.stdout)

    def test_cmake_capability(self):
        cmake = shutil.which("cmake")
        if not cmake:
            self.skipTest("CMake not available")
        cases = [
            ("unmarked", {}, False, True),
            ("marked", {"GAUXC_HAS_ONEDFT_GRADIENT_FIX": "ON"}, True, True),
            ("override", {"CP2K_GAUXC_ASSUME_ONEDFT_GRADIENT_FIX": "ON"}, True, True),
            (
                "no-model",
                {"GAUXC_HAS_ONEDFT": "OFF", "GAUXC_HAS_ONEDFT_GRADIENT_FIX": "ON"},
                False,
                True,
            ),
            (
                "invalid-override",
                {
                    "GAUXC_HAS_ONEDFT": "OFF",
                    "CP2K_GAUXC_ASSUME_ONEDFT_GRADIENT_FIX": "ON",
                },
                False,
                False,
            ),
            (
                "disabled",
                {"CP2K_USE_GAUXC": "OFF", "GAUXC_HAS_ONEDFT_GRADIENT_FIX": "ON"},
                False,
                True,
            ),
        ]
        for name, overrides, enabled, success in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmp:
                variables = {
                    "CP2K_USE_GAUXC": "ON",
                    "GAUXC_HAS_ONEDFT": "ON",
                    **overrides,
                }
                script = Path(tmp) / "check.cmake"
                expected = "ON" if enabled else "OFF"
                script.write_text(
                    "cmake_minimum_required(VERSION 3.22)\n"
                    + "\n".join(
                        f"set({key} {value})" for key, value in variables.items()
                    )
                    + f'\ninclude("{ROOT / "cmake/GauXCGradientSupport.cmake"}")\n'
                    + f'if(NOT "${{CP2K_GAUXC_HAS_ONEDFT_GRADIENT_FIX}}" STREQUAL "{expected}")\n'
                    + '  message(FATAL_ERROR "wrong capability")\nendif()\n'
                    # Rechecking a different, unmarked package must not retain the marker.
                    + "unset(GAUXC_HAS_ONEDFT_GRADIENT_FIX)\n"
                    + "set(CP2K_GAUXC_ASSUME_ONEDFT_GRADIENT_FIX OFF)\n"
                    + f'include("{ROOT / "cmake/GauXCGradientSupport.cmake"}")\n'
                    + "if(CP2K_GAUXC_HAS_ONEDFT_GRADIENT_FIX)\n"
                    + '  message(FATAL_ERROR "stale capability")\nendif()\n'
                )
                result = subprocess.run(
                    [cmake, "-P", str(script)], capture_output=True, text=True
                )
                self.assertEqual(
                    result.returncode == 0, success, result.stdout + result.stderr
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
