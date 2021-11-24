#!/usr/bin/env python3
#  vim: set ts=4 sw=4 tw=0 :

import unittest
from tempfile import NamedTemporaryFile

import numpy as np

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


class TestForceEnvironment(unittest.TestCase):
    def setUp(self):
        self._input = NamedTemporaryFile(mode="w+")
        self._input.write(TEST_FILE_CONTENT)
        self._input.flush()
        self._output = NamedTemporaryFile()
        self._fenv = cp2k.ForceEnvironment(self._input.name, self._output.name)

    def tearDown(self):
        self._fenv.destroy()
        self._input.close()
        self._output.close()

    def test_numbers(self):
        self.assertEqual(self._fenv.natom, 2)
        self.assertEqual(self._fenv.nparticle, 2)

    def test_calc_energy(self):
        self._fenv.calc_energy()
        self.assertLess(self._fenv.potential_energy, 0)

    def test_calc_forces(self):
        self._fenv.calc_energy_force()
        self.assertTrue(np.any(abs(self._fenv.forces) > 0.0))

    def test_positions(self):
        old_pos = self._fenv.positions
        zeros = np.zeros(old_pos.shape, dtype=np.double)
        self._fenv.positions = zeros
        self.assertTrue(np.all(self._fenv.positions == zeros))
        self._fenv.positions = old_pos


if __name__ == "__main__":
    unittest.main()
