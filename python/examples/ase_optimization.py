"""A complete, deliberately small ASE optimization through libcp2k."""

# SPDX-License-Identifier: GPL-2.0-or-later

from ase import Atoms
from ase.optimize import BFGS

from cp2k import CP2K
from cp2k.ase import CP2KCalculator

atoms = Atoms("H2", positions=[[3.6, 4, 4], [4.4, 4, 4]], cell=[8] * 3, pbc=True)
inp = {
    "FORCE_EVAL": {
        "METHOD": "Quickstep",
        "DFT": {
            "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
            "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
            "MGRID": {"CUTOFF": 200},
            "SCF": {
                "EPS_SCF": 1e-8,
                "MAX_SCF": 100,
                "OT": {"MINIMIZER": "DIIS"},
                "PRINT": {"RESTART": {"_": "OFF"}},
            },
            "XC": {"XC_FUNCTIONAL": {"_": "PADE"}},
        },
        "SUBSYS": {
            "KIND": {
                "_": "H",
                "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                "POTENTIAL": "GTH-PADE-q1",
            }
        },
    }
}
with CP2K() as cp:
    with CP2KCalculator(cp, inp, label="h2-optimization") as calc:
        atoms.calc = calc
        optimizer = BFGS(atoms, trajectory="h2.traj")
        optimizer.run(fmax=0.05, steps=20)
        print("Energy [eV]:", atoms.get_potential_energy())
        print("Bond length [angstrom]:", atoms.get_distance(0, 1))
