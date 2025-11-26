"""
Physical parameter calculations for spectrum methods
"""

import numpy as np

from utils.constants import kB


def calculate_alpha_parameter(requested_states, mode_count, displacements, frequencies):
    """
    Calculate alpha parameter for LQ methods.
    alpha = sum (displacement_j^2 * frequency_j^2) / 4
    """
    alpha_list = []
    
    for state in requested_states:
        alpha = sum(
            (displacements[state][j] ** 2) * (frequencies[j] ** 2)
            for j in range(1, mode_count + 1)
        ) / 4
        alpha_list.append(alpha)
    
    return alpha_list


def calculate_gamma_parameter(requested_states, mode_count, displacements, frequencies):
    """
    Calculate gamma parameter for LQ3 method.
    gamma = (1/12) * sum (displacement_j^2 * frequency_j^3)
    """
    gamma_list = []
    
    for state in requested_states:
        gamma = (1 / 12) * sum(
            (displacements[state][j] ** 2) * (frequencies[j] ** 3)
            for j in range(1, mode_count + 1)
        )
        gamma_list.append(gamma)
    
    return gamma_list


def calculate_adiabatic_energies(requested_states, displacements, mode_count, frequencies, vertical_energies):
    """
    Calculate adiabatic excitation energies.
    E_adiabatic = E_vertical - sum (frequency_j/2 * displacement_j^2)
    """
    adiabatic_list = []
    
    for i, state in enumerate(requested_states):
        total_relaxation = sum(
            (frequencies[j] / 2) * (displacements[state][j] ** 2)
            for j in range(1, mode_count + 1)
        )
        adiabatic_energy = vertical_energies[i] - total_relaxation
        adiabatic_list.append(adiabatic_energy)
    
    return adiabatic_list


def calculate_huang_rhys_factors(state, displacements, mode_count):
    """
    Calculate Huang-Rhys factors for IMDHO method.    
    S_j = displacement_j^2 / 2
    """
    huang_rhys_list = []
    
    for j in range(1, mode_count + 1):
        huang_rhys = (displacements[state][j] ** 2) / 2
        huang_rhys_list.append(huang_rhys)
    
    return huang_rhys_list


def calculate_thermal_factors(frequencies, mode_count, temperature):
    """
    Calculate thermal occupation factors for IMDHO method.    
    a_j = 2 / [exp(frequency_j / (k_B * T)) - 1] + 1
    """
    
    thermal_list = []
    
    for j in range(mode_count):
        if temperature > 0:
            exponent = frequencies[j] / (kB * temperature)
            thermal_factor = 2 / (np.exp(exponent) - 1) + 1
        else:
            # at T=0, all factors are 1
            thermal_factor = 1.0 
        
        thermal_list.append(thermal_factor)
    
    return thermal_list