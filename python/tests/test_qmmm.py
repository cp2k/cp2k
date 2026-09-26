"""Check subtraction, region mapping and variational energy/force/virial signs."""

# SPDX-License-Identifier: GPL-2.0-or-later

import numpy as np
import pytest

from cp2k import CalculationResult
from cp2k.qmmm import SubtractiveQMMM


class HarmonicModel:
    communicator = None
    natom = nparticle = 2

    def __init__(self, k):
        self.k = k
        self.closed = False

    def _check(self):
        if self.closed:
            raise RuntimeError("closed model")

    def calculate(self, *, forces=True, stress=False):
        virial = -self.k * self.positions.T @ self.positions
        return CalculationResult(
            self.k * np.sum(self.positions**2) / 2,
            -self.k * self.positions if forces else None,
            virial / np.linalg.det(self.cell) if stress else None,
            virial if stress else None,
        )


def correction():
    return SubtractiveQMMM(
        HarmonicModel(3),
        HarmonicModel(1),
        qm_atoms=[2, 0],
        positions=[[1, 2, 3], [2, 4, 6], [3, 2, 1]],
        cell=np.eye(3) * 10,
    )


def test_qmmm_derivatives():
    model = correction()
    positions, cell = model.positions, model.cell
    result = model.calculate(stress=True)
    np.testing.assert_array_equal(result.forces[1], 0)
    np.testing.assert_allclose(result.forces[[2, 0]], -2 * positions[[2, 0]])
    assert result.energy == pytest.approx(np.sum(positions[[2, 0]] ** 2))
    delta = 1e-5
    for i in range(3):
        for j in range(3):
            values = []
            for sign in (1, -1):
                perturbed = positions.copy()
                perturbed[i, j] += sign * delta
                model.positions = perturbed
                values.append(model.calculate(forces=False).energy)
            assert result.forces[i, j] == pytest.approx(
                -(values[0] - values[1]) / (2 * delta), abs=1e-8
            )
    for i, j in ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)):
        values = []
        for sign in (1, -1):
            deformation = np.eye(3)
            deformation[i, j] += sign * delta
            model.cell = cell @ deformation.T
            model.positions = positions @ deformation.T
            values.append(model.calculate(forces=False).energy)
        assert result.virial[i, j] == pytest.approx(
            -(values[0] - values[1]) / (2 * delta), abs=1e-8
        )


@pytest.mark.parametrize("selection", [[0, 0], [0, 3], [-1, 0], [0.5, 1], [], [0]])
def test_qmmm_selection(selection):
    with pytest.raises(ValueError):
        SubtractiveQMMM(
            HarmonicModel(3),
            HarmonicModel(1),
            qm_atoms=selection,
            positions=np.ones((3, 3)),
            cell=np.eye(3),
        )


def test_qmmm_ownership():
    model = correction()
    model.qm.closed = True
    with pytest.raises(RuntimeError, match="closed"):
        model.calculate()


@pytest.mark.parametrize("platform", ["Reference", "CPU"])
def test_openmm_subtractive(platform):
    openmm = pytest.importorskip("openmm", minversion="8.6.1")
    from cp2k.openmm import create_force
    from cp2k._units import BOHR_TO_NM, HARTREE_TO_KJMOL

    model = correction()
    system = openmm.System()
    # Host MM: k=1 for all three atoms; replace k=1 by k=3 only in the QM region.
    mm = openmm.CustomExternalForce("0.5*k*(x*x+y*y+z*z)")
    mm.addGlobalParameter("k", HARTREE_TO_KJMOL / BOHR_TO_NM**2)
    for i in range(3):
        system.addParticle(1)
        mm.addParticle(i, [])
    system.addForce(mm)
    system.addForce(create_force(model, periodic=False))
    integrator = openmm.VerletIntegrator(0.000001)
    context = openmm.Context(
        system, integrator, openmm.Platform.getPlatformByName(platform)
    )
    context.setPositions(model.positions * BOHR_TO_NM)
    state = context.getState(getEnergy=True, getForces=True)
    expected_energy = np.sum(model.positions**2) / 2 + model.calculate().energy
    expected_forces = -model.positions + model.calculate().forces
    assert state.getPotentialEnergy()._value == pytest.approx(
        expected_energy * HARTREE_TO_KJMOL, rel=1e-6
    )
    np.testing.assert_allclose(
        state.getForces(asNumpy=True)._value,
        expected_forces * HARTREE_TO_KJMOL / BOHR_TO_NM,
        rtol=1e-6,
    )
    del context, integrator
