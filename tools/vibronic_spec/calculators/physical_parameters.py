"""
Physical parameter calculations for spectrum methods
"""

import numpy as np
from typing import Dict, List

from constants import kB


def calculate_alpha_parameter(
    requested_states: List[int],
    mode_count: int,
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
) -> List[float]:
    """
    Calculate alpha parameter for LQ methods
    alpha = sum (displacement_j^2 * frequency_j^2) / 4
    """
    alphas: List[float] = []

    for state in requested_states:
        alpha = (
            sum(
                (displacements[state][j] ** 2) * (frequencies[j] ** 2)
                for j in range(1, mode_count + 1)
            )
            / 4
        )
        alphas.append(alpha)

    return alphas


def calculate_gamma_parameter(
    requested_states: List[int],
    mode_count: int,
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
) -> List[float]:
    """
    Calculate gamma parameter for LQ3 method
    gamma = (1/12) * sum (displacement_j^2 * frequency_j^3)
    """
    gammas: List[float] = []

    for state in requested_states:
        gamma = (1 / 12) * sum(
            (displacements[state][j] ** 2) * (frequencies[j] ** 3)
            for j in range(1, mode_count + 1)
        )
        gammas.append(gamma)

    return gammas


def calculate_adiabatic_energies(
    requested_states: List[int],
    displacements: Dict[int, Dict[int, float]],
    mode_count: int,
    frequencies: Dict[int, float],
    vertical_energies: List[float],
) -> List[float]:
    """
    Calculate adiabatic excitation energies
    E_adiabatic = E_vertical - sum (frequency_j/2 * displacement_j^2)
    """
    adiabatic_list: List[float] = []

    for i, state in enumerate(requested_states):
        total_relaxation = sum(
            (frequencies[j] / 2) * (displacements[state][j] ** 2)
            for j in range(1, mode_count + 1)
        )
        adiabatic_energy = vertical_energies[i] - total_relaxation
        adiabatic_list.append(adiabatic_energy)

    return adiabatic_list


def calculate_huang_rhys_factors(
    state: int, displacements: Dict[int, Dict[int, float]], mode_count: int
) -> List[float]:
    """
    Calculate Huang-Rhys factors for IMDHO method
    S_j = displacement_j^2 / 2
    """
    huang_rhys_list: List[float] = []

    for j in range(1, mode_count + 1):
        huang_rhys = (displacements[state][j] ** 2) / 2
        huang_rhys_list.append(huang_rhys)

    return huang_rhys_list


def calculate_thermal_factors(
    frequencies: List[float], mode_count: int, temperature: float
) -> List[float]:
    """
    Calculate thermal occupation factors for IMDHO method
    a_j = 2 / [exp(frequency_j / (k_B * T)) - 1] + 1
    """

    thermal_list: List[float] = []

    for j in range(mode_count):
        if temperature > 0:
            exponent = frequencies[j] / (kB * temperature)
            thermal_factor = 2 / (np.exp(exponent) - 1) + 1
        else:
            # at T=0, all factors are 1
            thermal_factor = 1.0

        thermal_list.append(thermal_factor)

    return thermal_list
