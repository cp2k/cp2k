#!/usr/bin/env python3

import unittest
from tempfile import NamedTemporaryFile

import cp2k

TEST_FILE_CONTENT = """
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
    POTENTIAL_FILE_NAME POTENTIAL
    LSD
    &MGRID
      CUTOFF 140
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-8
    &END QS
    &SCF
      EPS_DIIS 0.1
      EPS_SCF 1.0E-4
      MAX_DIIS 4
      MAX_SCF 3
      SCF_GUESS atomic
      &PRINT
        &RESTART OFF
        &END
      &END
    &END SCF
    &XC
      &XC_FUNCTIONAL Pade
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 8.0 4.0 4.0
    &END CELL
    &COORD
    H     0.000000  0.000000  0.000000
    H     1.000000  0.000000  0.000000
    &END COORD
    &KIND H
      BASIS_SET DZV-GTH-PADE
      POTENTIAL GTH-PADE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
&GLOBAL
  PROJECT libcp2k_unittest_H2
&END GLOBAL
"""


def setUpModule():
    cp2k.init()


def tearDownModule():
    cp2k.finalize()


class TestBasic(unittest.TestCase):
    def setUp(self):
        self._input = NamedTemporaryFile(mode="w+")
        self._input.write(TEST_FILE_CONTENT)
        self._input.flush()
        self._output = NamedTemporaryFile()

    def tearDown(self):
        self._input.close()
        self._output.close()

    def test_version_string(self):
        self.assertIn("CP2K version", cp2k.get_version_string())

    def test_run_input(self):
        cp2k.run_input(self._input.name, self._output.name)


if __name__ == "__main__":
    unittest.main()
