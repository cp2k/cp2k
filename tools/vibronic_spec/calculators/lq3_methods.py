"""
LQ3 approximation methods
"""

import numpy as np
from typing import Any, Dict, List, Optional

from .lq2_methods import _calculate_vibrational_relaxation


def calculate_lq3_spectrum_point(
    energy: float,
    state_idx: int,
    spectrum_type: str,
    oscillator_strengths: Dict[int, Dict[str, float]],
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
    mode_count: int,
    alphas: List[float],
    gammas: List[float],
    vertical_energies: List[float],
    stokes_shift: float,
    adiabatic_energies: Optional[List[float]],
    requested_states: List[int],
    integration_params: Dict[str, Any],
) -> float:
    """
    Calculate LQ3 spectrum point for absorption or fluorescence

    Parameters:
    - energy: Energy point (atomic units)
    - state_idx: Index in state list (1-based)
    - spectrum_type: 'absorption' or 'fluorescence'
    - oscillator_strengths: Oscillator strength data
    - displacements: Displacement vectors
    - frequencies: Vibrational frequencies
    - mode_count: Number of vibrational modes
    - alphas: Alpha parameters for each state
    - gammas: Gamma parameters for each state
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

    integral = _integrate_lq3_time_domain(
        energy,
        state_idx,
        spectrum_type,
        displacements,
        frequencies,
        mode_count,
        alphas,
        gammas,
        vertical_energies,
        stokes_shift,
        adiabatic_energies,
        requested_states,
        integration_params,
    )

    return integral * oscillator_strength


def _integrate_lq3_time_domain(
    energy: float,
    state_idx: int,
    spectrum_type: str,
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
    mode_count: int,
    alphas: List[float],
    gammas: List[float],
    vertical_energies: List[float],
    stokes_shift: float,
    adiabatic_energies: Optional[List[float]],
    requested_states: List[int],
    integration_params: Dict[str, Any],
) -> float:
    """
    Time integration for LQ3
    """

    max_time_slices = integration_params.get("max_time_slices", 5000)
    slice_size = integration_params.get("slice_size", 30000)
    time_step = integration_params.get("time_step", 30)
    convergence = integration_params.get("convergence", 0.0000001)

    if spectrum_type == "absorption":
        # Absorption: E - E_vertical - Stokes/2
        energy_difference = energy - vertical_energies[state_idx - 1] - stokes_shift / 2
    elif spectrum_type == "fluorescence":
        if adiabatic_energies is None:
            raise ValueError("Adiabatic energies must be provided for fluorescence.")
        # Fluorescence: E_adiabatic - Stokes/2 - E - relaxation_energy
        adiabatic_energy = adiabatic_energies[state_idx - 1]
        relaxation_energy = _calculate_vibrational_relaxation(
            state_idx, mode_count, displacements, frequencies, requested_states
        )
        energy_difference = (
            adiabatic_energy - stokes_shift / 2 - energy - relaxation_energy
        )
    else:
        raise ValueError(f"Unknown spectrum_type: {spectrum_type}")

    # Determine optimal time step
    time_scale_1 = (
        abs(2 * np.pi / energy_difference) if energy_difference != 0 else float("inf")
    )
    time_scale_2 = (2 * np.pi / gammas[state_idx - 1]) ** (1 / 3)
    min_time_scale = min(time_scale_1, time_scale_2)
    delta_time = min_time_scale / time_step

    slice_index = 0
    integral_value = 0.0
    convergence_reached = False

    while slice_index < max_time_slices and not convergence_reached:
        time_array = np.arange(
            slice_index * slice_size * delta_time,
            (slice_index + 1) * slice_size * delta_time,
            delta_time,
        )

        integrand = np.exp(-alphas[state_idx - 1] * time_array**2) * np.cos(
            energy_difference * time_array + gammas[state_idx - 1] * time_array**3
        )

        slice_integral = float(np.trapz(integrand, dx=delta_time))

        if abs(slice_integral) < convergence:
            convergence_reached = True
            integral_value += slice_integral
        else:
            integral_value += slice_integral
            slice_index += 1

    if not convergence_reached:
        print(f"WARNING: LQ3 time integral for {spectrum_type} did not converge.")

    return integral_value


def calculate_lq3_absorption(
    energy: float,
    state_idx: int,
    oscillator_strengths: Dict[int, Dict[str, float]],
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
    mode_count: int,
    alphas: List[float],
    gammas: List[float],
    vertical_energies: List[float],
    stokes_shift: float,
    requested_states: List[int],
    integration_params: Dict[str, Any],
) -> float:
    """Wrapper for absorption"""
    return calculate_lq3_spectrum_point(
        energy,
        state_idx,
        "absorption",
        oscillator_strengths,
        displacements,
        frequencies,
        mode_count,
        alphas,
        gammas,
        vertical_energies,
        stokes_shift,
        None,
        requested_states,
        integration_params,
    )


def calculate_lq3_fluorescence(
    energy: float,
    state_idx: int,
    oscillator_strengths: Dict[int, Dict[str, float]],
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
    mode_count: int,
    alphas: List[float],
    gammas: List[float],
    vertical_energies: List[float],
    stokes_shift: float,
    adiabatic_energies: List[float],
    requested_states: List[int],
    integration_params: Dict[str, Any],
) -> float:
    """Wrapper for fluorescence"""
    return calculate_lq3_spectrum_point(
        energy,
        state_idx,
        "fluorescence",
        oscillator_strengths,
        displacements,
        frequencies,
        mode_count,
        alphas,
        gammas,
        vertical_energies,
        stokes_shift,
        adiabatic_energies,
        requested_states,
        integration_params,
    )
