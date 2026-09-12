#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
#   CP2K: A general program to perform molecular dynamics simulations
#   Copyright 2000-2026 CP2K developers group <https://cp2k.org>
#
#   SPDX-License-Identifier: GPL-2.0-or-later
# --------------------------------------------------------------------------------------------------

"""Canonical-ensemble reference, independent of CP2K trajectory sampling.

Three particles have U = k*(r12**2 + r32**2)/2 and xi = r12 - r32.
After removing the irrelevant translation/rotation factors, the unconstrained
configurational partition density is the one-dimensional radial integral
P(xi) = integral r12**2*r32**2*exp(-beta*U) ds, with positive bond lengths.
The reference is -kBT*d(log(P))/dxi, with no Blue Moon formula in P.

Constrained phase-space integration instead has configurational weight
sqrt(Z)*r12**2*r32**2*exp(-beta*U) ds d(cos(theta)). Velocities are integrated
analytically using their tangent-space Maxwell covariance. The multiplier follows
from twice differentiating xi: -lambda = (grad(xi).M^-1.grad(U) - u.H.u)/Z,
where u is the Cartesian velocity vector.
Neither its construction nor the reference partition integral calls the tool's
automatic differentiation or metric routines. Gaussian quadrature removes the
statistical/timestep uncertainty of an MD benchmark, but does not test CP2K's
integrator. The separate recorded-trajectory regression tests the CP2K I/O path.
"""

import math
import unittest

import numpy as np

from blue_moon import KB_HARTREE_PER_K, FloatArray, analyze


def radial_quadrature(
    target: float, order: int
) -> tuple[FloatArray, FloatArray, FloatArray]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    # k = kBT in atomic length units: exp(-beta*U) = exp(-(r1**2+r2**2)/2).
    # At s=12 the omitted tail is negligible for the tested targets.
    s = 6 * (nodes + 1)
    r1, r2 = s + max(target, 0), s + max(-target, 0)
    density = r1**2 * r2**2 * np.exp(-(r1**2 + r2**2) / 2)
    return r1, r2, 6 * weights * density


def reference_gradient(target: float, order: int, step: float) -> float:
    def log_partition(xi: float) -> float:
        return math.log(float(radial_quadrature(xi, order)[2].sum()))

    return -(log_partition(target + step) - log_partition(target - step)) / (2 * step)


def pair_derivatives(
    unit: FloatArray, radius: float, a: int, b: int
) -> tuple[FloatArray, FloatArray]:
    gradient = np.zeros((3, 3))
    gradient[a], gradient[b] = unit, -unit
    hessian = np.zeros((9, 9))
    block = (np.eye(3) - np.outer(unit, unit)) / radius
    for i, sign_i in ((a, 1), (b, -1)):
        for j, sign_j in ((a, 1), (b, -1)):
            hessian[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] = sign_i * sign_j * block
    return gradient.ravel(), hessian


class EnsembleTests(unittest.TestCase):
    def test_partition_function_gradient(self) -> None:
        kt = KB_HARTREE_PER_K * 300
        cosines, angular_weights = np.polynomial.legendre.leggauss(48)
        for target in (-1.0, 0.75):
            expected = reference_gradient(target, 64, 1e-4)
            # Check both quadrature order and the independent derivative step.
            self.assertAlmostEqual(
                expected, reference_gradient(target, 80, 5e-5), delta=1e-8
            )
            for masses in ([12.0, 1.0, 16.0], [1.0, 2.0, 3.0]):
                with self.subTest(target=target, masses=masses):
                    cfg = {
                        "cv": {
                            "type": "linear_combination",
                            "terms": [
                                {
                                    "coefficient": 1,
                                    "cv": {"type": "distance", "atoms": [1, 2]},
                                },
                                {
                                    "coefficient": -1,
                                    "cv": {"type": "distance", "atoms": [3, 2]},
                                },
                            ],
                        },
                        "cell_angstrom": (np.eye(3) * 40).tolist(),
                        "masses_amu": {str(i + 1): m for i, m in enumerate(masses)},
                        "temperature_kelvin": 300,
                        "target_au": target,
                        "target_tolerance_au": 1e-10,
                    }
                    inverse_mass = np.diag(np.repeat(1 / np.array(masses), 3))
                    frames: list[tuple[int, FloatArray]] = []
                    multipliers: list[float] = []
                    measure: list[float] = []
                    r1s, r2s, radial_weights = radial_quadrature(target, 48)
                    for r1, r2, wr in zip(r1s, r2s, radial_weights):
                        for cosine, wc in zip(cosines, angular_weights):
                            u1 = np.array([1.0, 0.0, 0.0])
                            u2 = np.array([cosine, math.sqrt(1 - cosine**2), 0.0])
                            grad1, hess1 = pair_derivatives(u1, float(r1), 0, 1)
                            grad2, hess2 = pair_derivatives(u2, float(r2), 2, 1)
                            gradient = grad1 - grad2
                            normal = inverse_mass @ gradient
                            z = float(gradient @ normal)
                            covariance = kt * (
                                inverse_mass - np.outer(normal, normal) / z
                            )
                            grad_u = kt * (r1 * grad1 + r2 * grad2)
                            curvature = float(np.sum((hess1 - hess2) * covariance))
                            minus_lambda = (float(normal @ grad_u) - curvature) / z
                            frames.append(
                                (
                                    len(frames) + 1,
                                    np.array([r1 * u1, np.zeros(3), r2 * u2]),
                                )
                            )
                            multipliers.append(-minus_lambda)
                            measure.append(float(wr * wc * math.sqrt(z)))

                    samples = np.array(list(analyze(cfg, frames, multipliers, 1)))
                    probability = np.array(measure)
                    weights = probability * samples[:, 5]
                    corrected = float(np.average(samples[:, 6], weights=weights)) / kt
                    without_g = float(np.average(-samples[:, 2], weights=weights)) / kt
                    without_z = (
                        float(np.average(samples[:, 6], weights=probability)) / kt
                    )
                    # The thermodynamic answer is mass independent, although Z,
                    # G, and the constrained multiplier distribution are not.
                    self.assertAlmostEqual(corrected, expected, delta=2e-7)
                    self.assertGreater(abs(without_g - expected), 1e-3)
                    self.assertGreater(abs(without_z - expected), 1e-3)
                    self.assertGreater(float(np.ptp(samples[:, 3])), 0.1)
                    self.assertGreater(float(np.max(np.abs(samples[:, 4]))), 1e-3)


if __name__ == "__main__":
    unittest.main()
