#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
#   CP2K: A general program to perform molecular dynamics simulations
#   Copyright 2000-2026 CP2K developers group <https://cp2k.org>
#
#   SPDX-License-Identifier: GPL-2.0-or-later
# --------------------------------------------------------------------------------------------------

import asyncio
import sys
import unittest
from types import SimpleNamespace
from typing import List, Optional, Tuple
from unittest.mock import Mock, patch

import do_regtest


class StopBeforeRunningTests(Exception):
    pass


class CpuCountTests(unittest.TestCase):
    def test_restricted_affinity(self) -> None:
        host_count = Mock(return_value=128)
        affinity = Mock(return_value=set(range(64, 96)))
        with patch.object(
            do_regtest,
            "os",
            SimpleNamespace(
                getenv=Mock(return_value=None),
                cpu_count=host_count,
                sched_getaffinity=affinity,
            ),
        ):
            self.assertEqual(do_regtest.cpu_count(), 32)
        affinity.assert_called_once_with(0)
        host_count.assert_not_called()

    def test_noncontiguous_affinity(self) -> None:
        affinity = Mock(return_value={2, 7, 31})
        host_count = Mock(return_value=128)
        with patch.object(
            do_regtest,
            "os",
            SimpleNamespace(
                getenv=Mock(return_value=None),
                cpu_count=host_count,
                sched_getaffinity=affinity,
            ),
        ):
            self.assertEqual(do_regtest.cpu_count(), 3)
        host_count.assert_not_called()

    def test_without_affinity_api(self) -> None:
        with patch.object(
            do_regtest,
            "os",
            SimpleNamespace(
                getenv=Mock(return_value=None), cpu_count=Mock(return_value=12)
            ),
        ):
            self.assertEqual(do_regtest.cpu_count(), 12)

    def test_affinity_query_failure(self) -> None:
        with patch.object(
            do_regtest,
            "os",
            SimpleNamespace(
                getenv=Mock(return_value=None),
                cpu_count=Mock(return_value=8),
                sched_getaffinity=Mock(side_effect=OSError("unavailable")),
            ),
        ):
            self.assertEqual(do_regtest.cpu_count(), 8)

    def test_unknown_cpu_count(self) -> None:
        with patch.object(
            do_regtest,
            "os",
            SimpleNamespace(
                getenv=Mock(return_value=None), cpu_count=Mock(return_value=None)
            ),
        ):
            self.assertEqual(do_regtest.cpu_count(), 1)

    def test_python_cpu_count_override(self) -> None:
        for override in ("6", "128"):
            with self.subTest(override=override):
                getenv = Mock(return_value=override)
                host_count = Mock(return_value=128)
                affinity = Mock(return_value=set(range(32)))
                with patch.object(
                    do_regtest,
                    "os",
                    SimpleNamespace(
                        getenv=getenv,
                        cpu_count=host_count,
                        sched_getaffinity=affinity,
                    ),
                ):
                    self.assertEqual(do_regtest.cpu_count(), int(override))
                getenv.assert_called_once_with("PYTHON_CPU_COUNT")
                affinity.assert_not_called()
                host_count.assert_not_called()

    def test_empty_python_cpu_count_uses_affinity(self) -> None:
        with patch.object(
            do_regtest,
            "os",
            SimpleNamespace(
                getenv=Mock(return_value=""),
                cpu_count=Mock(return_value=128),
                sched_getaffinity=Mock(return_value=set(range(32))),
            ),
        ):
            self.assertEqual(do_regtest.cpu_count(), 32)

    def test_cli_default_and_override(self) -> None:
        cases: Tuple[Tuple[Optional[str], List[str], int], ...] = (
            (None, [], 32),
            ("6", [], 6),
            ("6", ["--maxtasks", "8"], 8),
        )
        for override, options, expected in cases:
            with self.subTest(override=override, options=options):
                with (
                    patch.object(
                        do_regtest,
                        "os",
                        SimpleNamespace(
                            getenv=Mock(return_value=override),
                            cpu_count=Mock(return_value=128),
                            sched_getaffinity=Mock(return_value=set(range(32))),
                        ),
                    ),
                    patch.object(
                        do_regtest, "Config", side_effect=StopBeforeRunningTests
                    ) as config,
                    patch.object(sys, "argv", ["do_regtest.py", *options, ".", "pdbg"]),
                ):
                    with self.assertRaises(StopBeforeRunningTests):
                        asyncio.run(do_regtest.main())
                self.assertEqual(config.call_args.args[0].maxtasks, expected)


if __name__ == "__main__":
    unittest.main()
