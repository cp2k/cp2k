"""
Independent Mode Displaced Harmonic Oscillator (IMDHO) Method
"""

import numpy as np
from typing import Any, Dict, List, Optional

from .physical_parameters import calculate_huang_rhys_factors, calculate_thermal_factors


def calculate_imdho_spectrum_point(
    energy: float,
    state_idx: int,
    spectrum_type: str,
    oscillator_strengths: Dict[int, Dict[str, float]],
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
    mode_count: int,
    gamma: float,
    theta: float,
    temperature: float,
    vertical_energies: List[float],
    stokes_shift: float,
    adiabatic_energies: List[float],
    requested_states: List[int],
    integration_params: Dict[str, Any],
) -> float:
    """
    Calculate IMDHO spectrum point for absorption or fluorescence

    Parameters:
    - energy: Energy point (atomic units)
    - state_idx: Index in state list (1-based)
    - spectrum_type: 'absorption' or 'fluorescence'
    - oscillator_strengths: Oscillator strength data
    - displacements: Displacement vectors
    - frequencies: Vibrational frequencies
    - mode_count: Number of vibrational modes
    - gamma: Broadening parameter (atomic units)
    - theta: Inhomogeneous broadening parameter (atomic units)
    - temperature: Temperature (Kelvin)
    - vertical_energies: Vertical excitation energies
    - stokes_shift: Stokes shift
    - adiabatic_energies: Adiabatic energies (only needed for fluorescence)
    - requested_states: List of state numbers that will be included in the final spectrum
    - integration_params: Parameters for time integration

    Returns:
    - Spectrum intensity value
    """
    actual_state = requested_states[state_idx - 1]
    oscillator_strength = oscillator_strengths[actual_state]["oscillator_strength"]

    if spectrum_type == "absorption":
        prefactor = (4 * energy * np.pi**2 / (3 * 137)) * (
            oscillator_strength / vertical_energies[state_idx - 1]
        )
        integral = _integrate_imdho_time_domain(
            energy,
            state_idx,
            spectrum_type,
            displacements,
            frequencies,
            mode_count,
            gamma,
            theta,
            temperature,
            stokes_shift,
            adiabatic_energies,
            requested_states,
            integration_params,
        )

    elif spectrum_type == "fluorescence":
        prefactor = (4 * energy**3 / (3 * np.pi * 137**3)) * (
            oscillator_strength / vertical_energies[state_idx - 1]
        )
        integral = _integrate_imdho_time_domain(
            energy,
            state_idx,
            spectrum_type,
            displacements,
            frequencies,
            mode_count,
            gamma,
            theta,
            temperature,
            stokes_shift,
            adiabatic_energies,
            requested_states,
            integration_params,
        )
    else:
        raise ValueError(f"Unknown spectrum_type: {spectrum_type}")

    return float(prefactor * integral)


def _integrate_imdho_time_domain(
    energy: float,
    state_idx: int,
    spectrum_type: str,
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
    mode_count: int,
    gamma: float,
    theta: float,
    temperature: float,
    stokes_shift: float,
    adiabatic_energies: List[float],
    requested_states: List[int],
    integration_params: Dict[str, Any],
) -> float:
    """
    Time integration for IMDHO
    """
    max_time_slices = integration_params.get("max_time_slices", 5000)
    slice_size = integration_params.get("slice_size", 10000)
    time_step = integration_params.get("time_step", 30)
    convergence = integration_params.get("convergence", 0.000000001)

    if spectrum_type == "absorption":
        energy_difference = (
            energy - adiabatic_energies[state_idx - 1] - stokes_shift / 2
        )
    elif spectrum_type == "fluorescence":
        # if adiabatic_energies is None:
        #     raise ValueError("Adiabatic energies must be provided for fluorescence.")
        energy_difference = (
            adiabatic_energies[state_idx - 1] - stokes_shift / 2 - energy
        )

    actual_state = requested_states[state_idx - 1]
    huang_rhys_factors = calculate_huang_rhys_factors(
        actual_state, displacements, mode_count
    )
    frequencies_list = [frequencies[j] for j in range(1, mode_count + 1)]

    if temperature > 0:
        thermal_factors = calculate_thermal_factors(
            frequencies_list, mode_count, temperature
        )
    else:
        thermal_factors = [1.0] * mode_count

    if energy_difference == 0:
        delta_time = 2 * np.pi / time_step
    else:
        time_scale_1 = abs(2 * np.pi / energy_difference)
        time_scale_2 = 2 * np.pi / max(frequencies)
        delta_time = min(time_scale_1, time_scale_2) / time_step

    slice_index = 0
    integral_value = 0.0
    convergence_reached = False

    while slice_index < max_time_slices and not convergence_reached:
        time_array = np.arange(
            slice_index * slice_size * delta_time,
            (slice_index + 1) * slice_size * delta_time,
            delta_time,
        )

        # Complex sum over vibrational modes
        cosine_sum = np.zeros(time_array.size)
        sine_sum = np.zeros(time_array.size)

        for mode_idx in range(mode_count):
            frequency = frequencies_list[mode_idx]
            huang_rhys = huang_rhys_factors[mode_idx]
            thermal_factor = thermal_factors[mode_idx]

            # Common cosine term for both
            if temperature > 0:
                cosine_sum += (
                    huang_rhys * thermal_factor * (1 - np.cos(frequency * time_array))
                )
            else:
                cosine_sum += huang_rhys * (1 - np.cos(frequency * time_array))

            # Different sign for sine term
            if spectrum_type == "absorption":
                sine_sum += -huang_rhys * np.sin(frequency * time_array)
            elif spectrum_type == "fluorescence":
                sine_sum += huang_rhys * np.sin(frequency * time_array)

        cosine_term = np.cos(energy_difference * time_array)
        exponential_term = np.exp(
            -gamma * time_array - 0.5 * theta**2 * time_array**2 - cosine_sum
        )
        # real part of exp(1j * sine_sum)
        phase_term = np.cos(sine_sum)

        integrand = cosine_term * exponential_term * phase_term

        if spectrum_type == "absorption":
            integrand *= 1 / np.pi

        slice_integral = float(np.trapz(integrand, dx=delta_time))

        if abs(slice_integral) < convergence:
            convergence_reached = True
            integral_value += slice_integral
        else:
            integral_value += slice_integral
            slice_index += 1

    if not convergence_reached:
        print(f"WARNING: IMDHO time integral for {spectrum_type} did not converge.")

    return float(integral_value)


# def calculate_imdho_absorption(
#     energy: float,
#     state_idx: int,
#     oscillator_strengths: Dict[int, Dict[str, float]],
#     displacements: Dict[int, Dict[int, float]],
#     frequencies: Dict[int, float],
#     mode_count: int,
#     gamma: float,
#     theta: float,
#     temperature: float,
#     vertical_energies: List[float],
#     stokes_shift: float,
#     adiabatic_energies: List[float],
#     requested_states: List[int],
#     integration_params: Dict[str, Any],
# ) -> float:
#     """Wrapper for IMDHO absorption"""
#     return calculate_imdho_spectrum_point(
#         energy,
#         state_idx,
#         "absorption",
#         oscillator_strengths,
#         displacements,
#         frequencies,
#         mode_count,
#         gamma,
#         theta,
#         temperature,
#         vertical_energies,
#         stokes_shift,
#         adiabatic_energies,
#         requested_states,
#         integration_params,
#     )


# def calculate_imdho_fluorescence(
#     energy: float,
#     state_idx: int,
#     oscillator_strengths: Dict[int, Dict[str, float]],
#     displacements: Dict[int, Dict[int, float]],
#     frequencies: Dict[int, float],
#     mode_count: int,
#     gamma: float,
#     theta: float,
#     temperature: float,
#     vertical_energies: List[float],
#     stokes_shift: float,
#     adiabatic_energies: List[float],
#     requested_states: List[int],
#     integration_params: Dict[str, Any],
# ) -> float:
#     """Wrapper for IMDHO fluorescence"""
#     return calculate_imdho_spectrum_point(
#         energy,
#         state_idx,
#         "fluorescence",
#         oscillator_strengths,
#         displacements,
#         frequencies,
#         mode_count,
#         gamma,
#         theta,
#         temperature,
#         vertical_energies,
#         stokes_shift,
#         adiabatic_energies,
#         requested_states,
#         integration_params,
#     )
