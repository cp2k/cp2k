#!/usr/bin/env python3
#  vim: set ts=4 sw=4 tw=0 :

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


def sample_run(input_file_path):

    with cp2k.ForceEnvironment(input_file_path, "/dev/null") as fenv:
        print("potential energy: {:e}".format(fenv.potential_energy))
        print("calculating energy..")
        fenv.calc_energy()
        print(".. done!")
        print("potential energy: {:e}".format(fenv.potential_energy))
        print("positions:")
        pos = fenv.positions
        print(pos)
        print("zeroify positions..")
        fenv.positions = np.zeros(pos.shape)
        print("positions:")
        print(fenv.positions)


if __name__ == "__main__":
    cp2k.init()

    tempfile = NamedTemporaryFile(mode="w+")
    tempfile.write(TEST_FILE_CONTENT)
    tempfile.flush()

    sample_run(tempfile.name)

    cp2k.finalize()
