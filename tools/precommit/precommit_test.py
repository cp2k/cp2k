#!/usr/bin/env python3

import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import precommit
from check_file_properties import check_file


class TestFilePropertyPaths(unittest.TestCase):
    def test_explicit_and_discovered_paths(self) -> None:
        for relative in ("src/grpp/grpp_binomial.c", "src/wignernj/gaunt.c"):
            for filename in (relative, "./" + relative, os.path.abspath(relative)):
                with self.subTest(filename=filename):
                    with patch.object(precommit, "run_local_tool") as run:
                        precommit.run_check_file_properties(filename)
                    run.assert_called_once_with(
                        "./tools/precommit/check_file_properties.py", relative
                    )

    def test_bundled_source_banner(self) -> None:
        root = Path(__file__).resolve().parents[2]
        old_cwd = Path.cwd()
        try:
            os.chdir(root)
            for relative in (
                "src/grpp/grpp_binomial.c",
                "src/wignernj/gaunt.c",
                "LICENSE",
                "src/grid/LICENSE",
                "src/dbm/LICENSE",
                "src/offload/LICENSE",
                "src/grpp/LICENSE",
                "src/wignernj/LICENSE",
            ):
                with self.subTest(filename=relative):
                    precommit.run_check_file_properties(str(root / relative))
        finally:
            os.chdir(old_cwd)

    def test_license_ownership(self) -> None:
        old_cwd = Path.cwd()
        with tempfile.TemporaryDirectory() as directory:
            try:
                os.chdir(directory)
                license_text = (
                    "Copyright (C) 1989, 1991 Free Software Foundation, Inc.\n"
                )
                for filename in ("LICENSE", "src/grpp/LICENSE", "src/wignernj/LICENSE"):
                    path = Path(filename)
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(license_text, encoding="utf8")
                    with self.subTest(filename=filename):
                        self.assertEqual(check_file(path), [])
                for filename in (
                    "src/dbm/LICENSE",
                    "src/grid/LICENSE",
                    "src/offload/LICENSE",
                ):
                    path = Path(filename)
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(
                        "Copyright 2000-2001 CP2K developers group\n", encoding="utf8"
                    )
                    with self.subTest(filename=filename):
                        self.assertIn(
                            f"{path}: Copyright banner malformed", check_file(path)
                        )
            finally:
                os.chdir(old_cwd)


if __name__ == "__main__":
    unittest.main()
