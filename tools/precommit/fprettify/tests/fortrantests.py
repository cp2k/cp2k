#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
#    This file is part of fprettify.
#    Copyright (C) 2016-2019 Patrick Seewald, CP2K developers group
#
#    fprettify is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    fprettify is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with fprettify. If not, see <http://www.gnu.org/licenses/>.
###############################################################################


import configparser
import difflib
import hashlib
import io
import logging
import os
import re
import shlex
import shutil
import sys
from datetime import datetime

import fprettify
from fprettify.tests.test_common import _MYPATH, FprettifyTestCase, joinpath

_TIMESTAMP = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

# main directory for running tests
TEST_MAIN_DIR = joinpath(_MYPATH, r"../../fortran_tests")

# directory for external Fortran code
TEST_EXT_DIR = joinpath(TEST_MAIN_DIR, r"test_code")

# directory containing Fortran examples
EXAMPLE_DIR = joinpath(_MYPATH, r"../../examples/in")

# backup directory
BACKUP_DIR = joinpath(TEST_MAIN_DIR, r"test_code_in_" + _TIMESTAMP)

# where to store summarized results
RESULT_DIR = joinpath(TEST_MAIN_DIR, r"test_results")

# expected hash-sums
RESULT_FILE = joinpath(RESULT_DIR, r"expected_results")

# test failures
FAILED_FILE = joinpath(RESULT_DIR, r"failed_results")


fprettify.set_fprettify_logger(logging.ERROR)


class FprettifyIntegrationTestCase(FprettifyTestCase):
    def shortDescription(self):
        """don't print doc string of testmethod"""
        return None

    def setUp(self):
        """
        setUp to be recognized by unittest.
        We have large files to compare, raise the limit
        """
        self.maxDiff = None

    @classmethod
    def setUpClass(cls):
        """
        setUpClass to be recognized by unittest.
        """

        cls.n_success = 0
        cls.n_parsefail = 0
        cls.n_internalfail = 0
        cls.n_unexpectedfail = 0

        FprettifyIntegrationTestCase.eprint("-" * 70)
        FprettifyIntegrationTestCase.eprint("recognized Fortran files")
        FprettifyIntegrationTestCase.eprint(", ".join(fprettify.FORTRAN_EXTENSIONS))
        FprettifyIntegrationTestCase.eprint("-" * 70)
        FprettifyIntegrationTestCase.eprint(
            "Applying fprettify to Fortran files in " + TEST_EXT_DIR
        )
        FprettifyIntegrationTestCase.eprint(
            "Writing backup of original files to " + BACKUP_DIR
        )
        FprettifyIntegrationTestCase.eprint(
            "Storing expected results in " + RESULT_FILE
        )
        FprettifyIntegrationTestCase.eprint("Storing failed results in " + FAILED_FILE)
        FprettifyIntegrationTestCase.eprint("-" * 70)

    @classmethod
    def tearDownClass(cls):
        """
        tearDownClass to be recognized by unittest. Used for test summary
        output.
        """
        if cls.n_parsefail + cls.n_internalfail > 0:
            format_str = "{:<20}{:<6}"
            FprettifyIntegrationTestCase.eprint("\n" + "=" * 70)
            FprettifyIntegrationTestCase.eprint(
                "IGNORED errors: invalid or old Fortran"
            )
            FprettifyIntegrationTestCase.eprint("-" * 70)
            FprettifyIntegrationTestCase.eprint(
                format_str.format("parse errors: ", cls.n_parsefail)
            )
            FprettifyIntegrationTestCase.eprint(
                format_str.format("internal errors: ", cls.n_internalfail)
            )

    @staticmethod
    def write_result(filename, content, sep_str):  # pragma: no cover
        with io.open(filename, "a", encoding="utf-8") as outfile:
            outfile.write(sep_str.join(content) + "\n")

    @staticmethod
    def eprint(*args, **kwargs):
        """
        Print to stderr - to print output compatible with default unittest output.
        """

        print(*args, file=sys.stderr, flush=True, **kwargs)


def generate_suite(suite=None, name=None):
    # make sure all directories exist
    if not os.path.exists(TEST_EXT_DIR):  # pragma: no cover
        os.makedirs(TEST_EXT_DIR)
    if not os.path.exists(BACKUP_DIR):  # pragma: no cover
        os.makedirs(BACKUP_DIR)
    if not os.path.exists(RESULT_DIR):  # pragma: no cover
        os.makedirs(RESULT_DIR)
    if not os.path.exists(RESULT_FILE):  # pragma: no cover
        io.open(RESULT_FILE, "w", encoding="utf-8").close()
    if os.path.exists(FAILED_FILE):  # pragma: no cover
        # erase failures from previous testers
        io.open(FAILED_FILE, "w", encoding="utf-8").close()

    import git

    config = configparser.ConfigParser()
    config.read(joinpath(TEST_MAIN_DIR, "testsuites.config"))

    if suite is None and name is None:
        return None

    for key in config.sections():
        code = config[key]
        if code["suite"] == suite or key == name:
            orig = os.getcwd()
            try:
                os.chdir(TEST_EXT_DIR)

                if not os.path.isdir(code["path"]):
                    # code is reinitialized only if path doesn't exist
                    # this allows iterative debugging to compare changes between different versions of fprettify
                    print(f"obtaining {key} ...")
                    exec(code["obtain"])
            finally:
                os.chdir(orig)

            add_test_code(code["path"], code["options"])
    return FprettifyIntegrationTestCase


def normalize_line(line):
    """
    Normalize fortran line in a way that resulting string should be the same
    whether fprettify has been applied or not.
    """
    # fprettify might add missing ampersands when splitting string:
    line_out = re.sub("^\s*&", "", line.lower(), flags=re.MULTILINE)
    line_out = re.sub("&\s*$", "", line_out, flags=re.MULTILINE)
    # remove all whitespace characters (including newline)
    line_out = re.sub(r"\s", r"", line_out)
    return line_out


def add_test_code(code_path, options):
    print(f"creating test cases from {code_path} ...")
    # dynamically create test cases from fortran files in test directory

    parser = fprettify.get_arg_parser()
    args = parser.parse_args(shlex.split(options))
    fprettify_args = fprettify.process_args(args)

    for dirpath, _, filenames in os.walk(joinpath(TEST_EXT_DIR, code_path)):
        for example in [
            f
            for f in filenames
            if any(f.endswith(_) for _ in fprettify.FORTRAN_EXTENSIONS)
        ]:
            rel_dirpath = os.path.relpath(dirpath, start=TEST_EXT_DIR)

            include_file = True
            if args.exclude_max_lines is not None:
                line_count = 0
                with open(joinpath(dirpath, example)) as f:
                    for i in f:
                        line_count += 1
                        if line_count > args.exclude_max_lines:
                            include_file = False
                            break

            if include_file:
                add_test_method(
                    FprettifyIntegrationTestCase, rel_dirpath, example, fprettify_args
                )


def add_test_method(testcase, fpath, ffile, args):
    """add a test method for each example."""

    def testmethod(testcase):
        """this is the test method invoked for each example."""

        example_path = joinpath(TEST_EXT_DIR, fpath)
        backup_path = joinpath(BACKUP_DIR, fpath)
        if not os.path.exists(backup_path):
            os.makedirs(backup_path)

        example = joinpath(example_path, ffile)
        example_backup = joinpath(backup_path, ffile)

        def test_result(path, info):
            return [os.path.relpath(path, TEST_EXT_DIR), info]

        with io.open(example, "r", encoding="utf-8") as infile:
            instring = infile.read()

        # write backup of original file
        with io.open(example_backup, "w", encoding="utf-8") as outfile:
            outfile.write(instring)

        # initialize outstring containing reformatted file content
        outstring = instring

        # apply fprettify
        try:
            fprettify.reformat_inplace(example, **args)

            # update outstring
            with io.open(example, "r", encoding="utf-8") as outfile:
                outstring = outfile.read()

            m = hashlib.sha256()
            m.update(outstring.encode("utf-8"))

            test_info = "checksum"
            test_content = test_result(example, m.hexdigest())

            FprettifyIntegrationTestCase.n_success += 1
        except fprettify.FprettifyParseException as e:
            test_info = "parse error"
            fprettify.log_exception(e, test_info, level="warning")
            test_content = test_result(example, test_info)
            FprettifyIntegrationTestCase.n_parsefail += 1
        except fprettify.FprettifyInternalException as e:
            test_info = "internal error"
            fprettify.log_exception(e, test_info, level="warning")
            test_content = test_result(example, test_info)
            FprettifyIntegrationTestCase.n_internalfail += 1
        except:  # pragma: no cover
            FprettifyIntegrationTestCase.n_unexpectedfail += 1
            raise

        # check that no changes other than whitespace changes or lower/upper case occured
        orig_stripped = normalize_line(instring)
        new_stripped = normalize_line(outstring)

        testcase.assertEqual(
            orig_stripped,
            new_stripped,
            "fprettify caused changes other than whitespace or lower/upper case",
        )

        sep_str = " : "
        with io.open(RESULT_FILE, "r", encoding="utf-8") as infile:
            found = False
            for line in infile:
                line_content = line.strip().split(sep_str)
                if line_content[0] == test_content[0]:
                    found = True
                    FprettifyIntegrationTestCase.eprint(test_info, end=" ")
                    msg = "{} (old) != {} (new)".format(
                        line_content[1], test_content[1]
                    )
                    if test_info == "checksum" and outstring.count("\n") < 10000:
                        # difflib can not handle large files
                        result = list(
                            difflib.unified_diff(
                                instring.splitlines(True),
                                outstring.splitlines(True),
                                fromfile=test_content[0],
                                tofile=line_content[0],
                            )
                        )
                        msg += "\n" + "".join(result)
                    try:
                        testcase.assertEqual(line_content[1], test_content[1], msg)
                    except AssertionError:  # pragma: no cover
                        FprettifyIntegrationTestCase.write_result(
                            FAILED_FILE, test_content, sep_str
                        )
                        raise
                    break

        if not found:  # pragma: no cover
            FprettifyIntegrationTestCase.eprint(test_info + " new", end=" ")
            FprettifyIntegrationTestCase.write_result(
                RESULT_FILE, test_content, sep_str
            )

    # not sure why this even works, using "test something" (with a space) as function name...
    # however it gives optimal test output
    testmethod.__name__ = "test " + joinpath(fpath, ffile)

    setattr(testcase, testmethod.__name__, testmethod)
