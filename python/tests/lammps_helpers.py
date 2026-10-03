"""Shared two-atom setup and reference checks for native and socket LAMMPS tests."""

# SPDX-License-Identifier: GPL-2.0-or-later

from contextlib import contextmanager

import numpy as np

from cp2k._units import BOHR_TO_ANGSTROM, HARTREE_TO_EV, HARTREE_TO_KCALMOL


@contextmanager
def lammps_system(units, *, comm=None):
    from lammps import lammps

    lmp = lammps(comm=comm, cmdargs=["-log", "none", "-screen", "none"])
    try:
        lmp.commands_string(f"""
units {units}
atom_style atomic
boundary p p p
region box prism 0 20 0 21 0 22 1 2 3
create_box 1 box
create_atoms 1 single 7.2 4.8 4.5
create_atoms 1 single 4 4 4
mass 1 39.948
pair_style zero 8
pair_coeff * *
compute cp_pressure all pressure NULL virial
thermo_style custom step pe c_cp_pressure[1] c_cp_pressure[2] c_cp_pressure[3] c_cp_pressure[4] c_cp_pressure[5] c_cp_pressure[6]
thermo_modify norm no
""")
        yield lmp
    finally:
        lmp.close()


def check_lammps_reference(lmp, reference, cell, units):
    """Check unit conversion, reversed atom IDs and all six virial components."""
    factor = HARTREE_TO_EV if units == "metal" else HARTREE_TO_KCALMOL
    np.testing.assert_allclose(
        lmp.get_thermo("pe"), reference.energy * factor, rtol=1e-7, atol=1e-10
    )
    tags = lmp.numpy.extract_atom("id")[:2]
    indices = [1 if tag == 1 else 0 for tag in tags]
    np.testing.assert_allclose(
        lmp.numpy.extract_atom("f")[:2],
        reference.forces[indices] * factor / BOHR_TO_ANGSTROM,
        rtol=1e-8,
        atol=1e-10,
    )
    pressure = (
        reference.virial.flat[[0, 4, 8, 1, 2, 5]]
        * factor
        / np.linalg.det(cell * BOHR_TO_ANGSTROM)
        * lmp.extract_global("nktv2p")
    )
    np.testing.assert_allclose(
        lmp.numpy.extract_compute("cp_pressure", 0, 1), pressure, rtol=1e-8, atol=1e-10
    )


def run_lammps_npt(coupling, units):
    coupling.command("velocity all create 10 731 mom yes rot no dist gaussian")
    coupling.command("fix thermostat all npt temp 10 10 100 iso 0 0 1000")
    coupling.command("timestep " + ("0.0001" if units == "metal" else "0.1"))
    coupling.run(3)
    assert np.isfinite(coupling.lammps.get_thermo("pe"))
    coupling.command("unfix thermostat")
