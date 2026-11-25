#!/usr/bin/env python3
"""Small argument, topology and flock tests; real MPI checks use affinity_probe.c."""

from contextlib import ExitStack, redirect_stderr
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from cp2k_local import argument_parser, reserve, topology


class LocalRunTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory()))

    def claim(self, allowed, count):
        groups = [(0, 4), (1, 5), (2, 6), (3, 7)]
        return self.stack.enter_context(
            reserve(self.root / "pool", groups, allowed, count)
        )

    def test_arguments(self):
        cases = (
            ([], {}, (1, 1)),
            (["-n", "2", "-t", "3"], {}, (2, 3)),
            ([], {"MPI_RANKS": "2", "OMP_NUM_THREADS": "3"}, (2, 3)),
            (["-n", "2"], {"MPI_RANKS": "bad", "OMP_NUM_THREADS": "3"}, (2, 3)),
            (["-t", "3"], {"MPI_RANKS": "2", "OMP_NUM_THREADS": "4,2"}, (2, 3)),
            ([], {"MPI_RANKS": "2"}, (2, 1)),
            ([], {"OMP_NUM_THREADS": "3"}, (1, 3)),
            ([], {"OMPI_COMM_WORLD_SIZE": "8"}, (1, 1)),
        )
        for cli, env, expected in cases:
            with (
                self.subTest(cli=cli, env=env),
                patch.dict(os.environ, env, clear=True),
            ):
                args = argument_parser().parse_args([*cli, "--", "cp2k.psmp"])
                self.assertEqual((args.ranks, args.threads), expected)

    def test_invalid_values_and_overrides(self):
        for option, variable in (
            ("--ranks", "MPI_RANKS"),
            ("--threads", "OMP_NUM_THREADS"),
        ):
            for value in ("", "abc", "0", "-1", "2,3", "1.5"):
                with (
                    self.subTest(option=option, value=value),
                    patch.dict(os.environ, {variable: value}, clear=True),
                ):
                    parser = argument_parser()
                    with (
                        redirect_stderr(io.StringIO()),
                        self.assertRaises(SystemExit) as error,
                    ):
                        parser.parse_args(["--", "cp2k.psmp"])
                    self.assertEqual(error.exception.code, 2)
                    args = parser.parse_args([option, "2", "--", "cp2k.psmp"])
                    self.assertEqual(getattr(args, option[2:]), 2)
                    with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                        parser.parse_args([option, value, "--", "cp2k.psmp"])

    def test_siblings_and_rollback(self):
        a = self.claim({0, 1, 2, 3}, 2)
        self.assertEqual(a[1], [0, 1])
        b = self.claim({4, 5, 6, 7}, 2)
        self.assertEqual(b[1], [6, 7])
        with self.assertRaisesRegex(ValueError, "only 0 available"):
            self.claim(set(range(8)), 1)
        for lock in a[2] + b[2]:
            lock.close()
        with self.assertRaisesRegex(ValueError, "only 4 available"):
            self.claim(set(range(8)), 5)
        self.assertEqual(self.claim(set(range(8)), 4)[1], [0, 1, 2, 3])

    def test_single_allowed_core(self):
        self.assertEqual(self.claim({4}, 1)[1], [4])
        with self.assertRaisesRegex(ValueError, "only 0 available"):
            self.claim({0}, 1)

    def test_private_pool(self):
        pool = self.root / "pool"
        pool.mkdir()
        pool.chmod(0o755)
        with self.assertRaisesRegex(ValueError, "0700"):
            self.claim({0, 1}, 1)
        pool.chmod(0o700)
        with patch("cp2k_local.os.getuid", return_value=os.getuid() + 1):
            with self.assertRaisesRegex(ValueError, "owned by this user"):
                self.claim({0, 1}, 1)
        (pool / "allocation.lock").symlink_to(self.root / "victim")
        with self.assertRaises(OSError):
            self.claim({0, 1}, 1)
        self.assertFalse((self.root / "victim").exists())

    def test_nonadjacent_and_offline_siblings(self):
        (self.root / "online").write_text("0,1,4")
        for cpu, siblings in ((0, "0,4"), (4, "0,4"), (1, "1,5")):
            path = self.root / f"cpu{cpu}/topology"
            path.mkdir(parents=True)
            (path / "thread_siblings_list").write_text(siblings)
        online, groups = topology(self.root)
        self.assertEqual(online, {0, 1, 4})
        self.assertEqual(groups, [(0, 4), (1, 5)])
        with reserve(self.root / "pool", groups, online & {4, 5}, 1) as allocation:
            self.assertEqual(allocation[1], [4])
        (self.root / "cpu4/topology/thread_siblings_list").write_text("4")
        with self.assertRaisesRegex(ValueError, "Inconsistent"):
            topology(self.root)


if __name__ == "__main__":
    unittest.main()
