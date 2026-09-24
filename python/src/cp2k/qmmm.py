"""Subtractive, mechanically embedded QM/MM corrections with explicit regions."""

# SPDX-License-Identifier: GPL-2.0-or-later

import numpy as np

from .library import CalculationResult, _array


class SubtractiveQMMM:
    """Supply E_QM(region) - E_MM(region) on top of a full-system MM potential.

    ``qm`` and ``mm`` are caller-owned environments for the *same* selected
    atoms, in the order of zero-based ``qm_atoms``. ``mm`` must reproduce exactly
    the internal MM terms of this region in the host's full-system force field.
    Both implement the CP2K environment protocol in atomic units; independent
    CP2K servers can be used for concurrent QM and MM environments.

    The returned force is a correction, zero outside the region. Keep the
    full-system MM force enabled in OpenMM/LAMMPS. QM-MM interactions remain MM
    (mechanical embedding); this does NOT polarize QM with external charges.
    Select whole molecules: covalent cuts/link atoms, adaptive regions and
    mismatched electrostatic/exclusion conventions are not handled here. Region,
    charge, spin and reference force-field parameters are explicit user choices.
    """

    def __init__(self, qm, mm, *, qm_atoms, positions, cell):
        if qm is mm:
            raise ValueError("QM and MM need distinct environments")
        qm._check()
        mm._check()
        self.qm, self.mm = qm, mm
        self.communicator = qm.communicator
        if self.communicator is not None or mm.communicator is not None:
            from mpi4py import MPI

            if (
                self.communicator is None
                or mm.communicator is None
                or MPI.Comm.Compare(self.communicator, mm.communicator)
                not in (MPI.IDENT, MPI.CONGRUENT)
            ):
                raise ValueError("QM and MM require congruent client communicators")
        positions = np.asarray(positions)
        if positions.ndim != 2 or positions.shape[1] != 3 or not len(positions):
            raise ValueError("positions must have shape (N, 3), N > 0")
        self.natom = self.nparticle = len(positions)
        selected = np.asarray(qm_atoms)
        if (
            selected.ndim != 1
            or selected.dtype.kind not in "iu"
            or not len(selected)
            or np.any(selected < 0)
            or np.any(selected >= self.natom)
            or len(np.unique(selected)) != len(selected)
        ):
            raise ValueError("qm_atoms must be distinct in-range integer indices")
        if any(
            env.natom != len(selected) or env.nparticle != len(selected)
            for env in (qm, mm)
        ):
            raise ValueError(
                "QM and MM environments must contain exactly the selected atoms"
            )
        self._indices = selected.copy()
        self.positions = positions
        self.cell = cell

    def _check(self):
        self.qm._check()
        self.mm._check()

    @property
    def positions(self):
        self._check()
        return self._positions.copy()

    @positions.setter
    def positions(self, values):
        self._check()
        self._positions = _array(values, (self.nparticle, 3), "positions")

    @property
    def cell(self):
        self._check()
        return self._cell.copy()

    @cell.setter
    def cell(self, values):
        self._check()
        cell = _array(values, (3, 3), "cell")
        determinant = np.linalg.det(cell)
        if not np.isfinite(determinant) or determinant <= 0:
            raise ValueError("cell must be nonsingular and right-handed")
        self._cell = cell

    def calculate(self, *, forces=True, stress=False):
        self._check()
        results = []
        for environment in (self.qm, self.mm):
            environment.cell = self.cell
            environment.positions = self.positions[self._indices]
            results.append(environment.calculate(forces=forces, stress=stress))
        high, low = results
        energy = high.energy - low.energy
        if not np.isfinite(energy):
            raise RuntimeError("Non-finite QM/MM correction energy")
        correction = None
        if forces:
            correction = np.zeros((self.nparticle, 3))
            correction[self._indices] = _array(
                high.forces - low.forces, (len(self._indices), 3), "QM/MM forces"
            )
        virial = None
        if stress:
            virial = _array(high.virial - low.virial, (3, 3), "QM/MM virial")
        return CalculationResult(
            energy,
            correction,
            None if virial is None else virial / np.linalg.det(self.cell),
            virial,
        )
