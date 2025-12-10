"""
LQ2 approximation methods
"""

import numpy as np
from typing import Dict, List, Optional


def calculate_lq2_spectrum_point(
    energy: float,
    state_idx: int,
    spectrum_type: str,
    oscillator_strengths: Dict[int, Dict[str, float]],
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
    mode_count: int,
    alphas: List[float],
    vertical_energies: List[float],
    stokes: float,
    adiabatic_energies: Optional[List[float]],
    requested_states: List[int],
) -> float:
    """
    Calculate LQ2 spectrum point for absorption or fluorescence

    Parameters:
    - energy: Energy point (atomic units)
    - state_idx: Index in state list (1-based)
    - spectrum_type: 'absorption' or 'fluorescence'
    - oscillator_strengths: Oscillator strength data
    - displacements: Displacement vectors
    - frequencies: Vibrational frequencies
    - mode_count: Number of vibrational modes
    - alphas: Alpha parameters for each state
    - vertical_energies: Vertical excitation energies
    - stokes: Stokes shift
    - adiabatic_energies: Adiabatic energies (only needed for fluorescence)
    - requested_states: List of state numbers that will be included in the final spectrum

    Returns:
    - Spectrum intensity value
    """
    state = requested_states[state_idx - 1]
    oscillator_strength = oscillator_strengths[state]["oscillator_strength"]
    alpha = alphas[state_idx - 1]

    if spectrum_type == "absorption":
        # Absorption: E - E_vertical - Stokes/2
        energy_difference = energy - vertical_energies[state_idx - 1] - stokes / 2
    elif spectrum_type == "fluorescence":
        if adiabatic_energies is None:
            raise ValueError("Adiabatic energies must be provided for fluorescence.")
        # Fluorescence: E_adiabatic - Stokes/2 - E - relaxation_energy
        relaxation = _calculate_vibrational_relaxation(
            state_idx, mode_count, displacements, frequencies, requested_states
        )
        energy_difference = (
            adiabatic_energies[state_idx - 1] - stokes / 2 - energy - relaxation
        )
    else:
        raise ValueError(f"Unknown spectrum_type: {spectrum_type}")

    if alpha > 0:
        prefactor = 1.0 / np.sqrt(4 * np.pi * alpha)
    else:
        prefactor = 1.0

    sigma = float(
        oscillator_strength * prefactor * np.exp(-(energy_difference**2) / (4 * alpha))
    )
    return sigma


def _calculate_vibrational_relaxation(
    state_idx: int,
    mode_count: int,
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
    requested_states: List[int],
) -> float:
    """Calculate vibrational relaxation energy"""
    state = requested_states[state_idx - 1]
    frequency_array = np.array([frequencies[j] for j in range(1, mode_count + 1)])
    displacement_array = np.array(
        [displacements[state][j] for j in range(1, mode_count + 1)]
    )
    return float(np.sum((displacement_array**2 * frequency_array) / 2))


def calculate_lq2_absorption(
    energy: float,
    state_idx: int,
    oscillator_strengths: Dict[int, Dict[str, float]],
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
    mode_count: int,
    alphas: List[float],
    vertical_energies: List[float],
    stokes: float,
    requested_states: List[int],
) -> float:
    """Wrapper for absorption spectrum calculation"""
    return calculate_lq2_spectrum_point(
        energy,
        state_idx,
        "absorption",
        oscillator_strengths,
        displacements,
        frequencies,
        mode_count,
        alphas,
        vertical_energies,
        stokes,
        None,
        requested_states,
    )


def calculate_lq2_fluorescence(
    energy: float,
    state_idx: int,
    oscillator_strengths: Dict[int, Dict[str, float]],
    displacements: Dict[int, Dict[int, float]],
    frequencies: Dict[int, float],
    mode_count: int,
    alphas: List[float],
    vertical_energies: List[float],
    stokes: float,
    adiabatic_energies: List[float],
    requested_states: List[int],
) -> float:
    """Wrapper for fluorescence spectrum calculation"""
    return calculate_lq2_spectrum_point(
        energy,
        state_idx,
        "fluorescence",
        oscillator_strengths,
        displacements,
        frequencies,
        mode_count,
        alphas,
        vertical_energies,
        stokes,
        adiabatic_energies,
        requested_states,
    )
