#!/usr/bin/env python3

import unittest
import os
import tempfile
import shutil
import sys

from format_fortran import main
from prettify_cp2k import selftest

interface_cpp_content = """\
MODULE interface_cpp_test
   IMPLICIT NONE

CONTAINS

   SUBROUTINE initialize()
      INTERFACE
         SUBROUTINE initialize_aux() BIND(C, name='initialize_aux')
         END SUBROUTINE initialize_aux
      END INTERFACE

#if defined(__DLAF)
      CALL initialize_aux()
#endif
   END SUBROUTINE initialize
END MODULE interface_cpp_test
"""


class TestSingleFileFolder(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.fname = os.path.join(self.tempdir, "prettify_selftest.F")

        # create temporary file with example code
        with open(self.fname, "w", encoding="utf8") as fhandle:
            fhandle.write(selftest.content)

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_prettify(self):
        # call prettify, the return value should be 0 (OK)
        self.assertEqual(main([self.fname]), 0)

        # check if file was altered (it shouldn't)
        with open(self.fname, encoding="utf8") as fhandle:
            result = fhandle.read()

        self.assertEqual(result.splitlines(), selftest.content.splitlines())


class TestInterfaceCppDirective(unittest.TestCase):
    def test_prettify_preserves_interface(self):
        with tempfile.TemporaryDirectory() as tempdir:
            fname = os.path.join(tempdir, "prettify_interface.F")
            with open(fname, "w", encoding="utf8") as fhandle:
                fhandle.write(interface_cpp_content)

            self.assertEqual(main([fname]), 0)

            with open(fname, encoding="utf8") as fhandle:
                result = fhandle.read()

        self.assertEqual(result, interface_cpp_content)


if __name__ == "__main__":
    unittest.main()
