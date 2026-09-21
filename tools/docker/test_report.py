#!/usr/bin/env python3

# author: Thomas D. Kühne

import json
import os
import subprocess
import tempfile
import time
import unittest
from pathlib import Path


class ReportTest(unittest.TestCase):
    def test_report_output(self) -> None:
        dockerfile = Path(__file__).with_name("Dockerfile.test_conventions")
        command = next(
            line.removeprefix("CMD ")
            for line in dockerfile.read_text(encoding="utf8").splitlines()
            if line.startswith("CMD ")
        )
        shell = not command.startswith("[")
        args = command if shell else json.loads(command)
        with tempfile.TemporaryDirectory(prefix="cp2k-ci-report-") as workdir:
            report = Path(workdir) / "report.log"
            for age in (0, 60, 660, 86400):
                for status in ("OK", "FAILED"):
                    with self.subTest(age=age, status=status):
                        content = (
                            f"Test output\nSummary: Test result\nStatus: {status}\n"
                        )
                        report.write_text(content, encoding="utf8")
                        timestamp = time.time() - age
                        os.utime(report, (timestamp, timestamp))
                        result = subprocess.run(
                            args,
                            shell=shell,
                            cwd=workdir,
                            input="",
                            text=True,
                            capture_output=True,
                            timeout=5,
                            check=True,
                        )
                        self.assertEqual(result.stdout, content)
            report.unlink()
            result = subprocess.run(
                args,
                shell=shell,
                cwd=workdir,
                input="",
                text=True,
                capture_output=True,
                timeout=5,
            )
            self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()

# EOF
