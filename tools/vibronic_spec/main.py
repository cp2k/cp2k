#!/usr/bin/env python3

"""
Main entry point for Vibronic Spectroscopy Tools (VibronicSpec)

Author: Beliz Sertcan
"""

import sys
import numpy as np
import time
import os
from typing import Dict, Any, List, Union, Optional

from file_parsers import parse_excited_state_forces, parse_molden_file

from calculators.lq2_methods import calculate_lq2_spectrum_point
from calculators.lq3_methods import calculate_lq3_spectrum_point
from calculators.imdho_methods import calculate_imdho_spectrum_point
from calculators.physical_parameters import (
    calculate_alpha_parameter,
    calculate_gamma_parameter,
    calculate_adiabatic_energies,
)

from output_formatters import write_spectrum_output
from constants import EV_TO_AU, AU_TO_EV, ATOMIC_MASSES, E_MASS


def run_spectrum_calculation(config: Dict[str, Any]) -> None:
    """
    Main function that runs the complete spectrum calculation
    """
    time_start = time.time()

    print("Starting spectrum calculation...")

    data = load_calculation_data(config)

    data = calculate_displacement_vectors(data, config)

    data = calculate_physical_parameters(data, config)

    results = calculate_spectrum(data, config)

    write_spectrum_output(results, config)

    time_end = time.time()

    print("Calculation completed successfully!")
    print(f"INFO: Output written to: {config['output_filename']}")
    print(f"INFO: Spectrum calculation finished in {time_end - time_start:.2f} seconds")


def load_calculation_data(config: Dict[str, Any]) -> Dict[str, Any]:
    """Load all required data from input files"""
    print("Loading input data...")

    data: Dict[str, Any] = {}

    config_dir = config.get("config_dir", ".")
    force_file = os.path.join(config_dir, config["force_filename"])
    vibrations_file = os.path.join(config_dir, config["vibrations_filename"])

    tdforce_data = parse_excited_state_forces(force_file)
    data["oscillator_strengths"] = tdforce_data
    forces = {
        state: tdforce_data[state]["force"]
        for state in tdforce_data
        if "force" in tdforce_data[state]
    }

    if "states" not in config:
        raise ValueError(
            "'states' configuration is required. Accepted formats are 'all', [1,2,3], and 'threshold:0.01'."
        )

    final_states = parse_states_specification(config["states"], tdforce_data, forces)

    data["forces"] = {state: forces[state] for state in final_states}
    data["requested_states"] = final_states
    data["state_count"] = len(final_states)

    molden_data = parse_molden_file(vibrations_file)
    data.update(molden_data)

    print(f"INFO: Found {data['atom_count']} atoms")
    print(f"INFO: Found {data['mode_count']} vibrational modes")
    if data["negative_freq_warnings"]:
        print(
            f"\tWARNING: Removed {len(data['negative_freq_warnings'])} negative frequencies"
        )
        for freq in data["negative_freq_warnings"]:
            print(f"\tWARNING: Negative frequency {freq:.1f} cm^-1")

    return data


def parse_states_specification(
    states_config: Union[str, List[int]],
    oscillator_strengths: Dict[int, Dict[str, Any]],
    available_forces: Dict[int, Any],
) -> List[int]:
    """
    Parse states specification
    Returns list of state numbers to calculate spectrum
    """
    states_with_forces = list(available_forces.keys())

    if isinstance(states_config, list):
        requested = states_config
        missing_states = [
            state for state in requested if state not in states_with_forces
        ]
        if missing_states:
            raise ValueError(
                f"Requested states {missing_states} do not have force data available. "
                f"Available states with forces: {states_with_forces}"
            )
        print(f"Using explicitly specified states: {requested}")
        return requested

    elif isinstance(states_config, str):
        if states_config.lower() == "all":
            states_with_forces.sort()
            print(
                f"Using all {len(states_with_forces)} states with available force data: {states_with_forces}"
            )
            return states_with_forces

        elif states_config.startswith("threshold:"):
            try:
                threshold = float(states_config.split(":")[1])
            except (IndexError, ValueError):
                raise ValueError(f"Invalid threshold format. Use 'threshold:0.01'")

            selected_states = []
            for state_num, state_info in oscillator_strengths.items():
                if (
                    state_info["oscillator_strength"] >= threshold
                    and state_num in states_with_forces
                ):
                    selected_states.append(state_num)

            selected_states.sort()
            if selected_states:
                print(
                    f"INFO: Selected {len(selected_states)} states with oscillator strength >= {threshold}"
                )
                print(f"INFO: Selected states: {selected_states}")
                return selected_states
            else:
                raise ValueError(
                    f"No states meet oscillator threshold {threshold} with available forces"
                )

        else:
            raise ValueError(
                f"Invalid states specification: {states_config}. Use 'all', [1,2,5], or 'threshold:0.01'"
            )

    else:
        raise ValueError("States configuration must be a list or string")


def calculate_displacement_vectors(
    data: Dict[str, Any], config: Dict[str, Any]
) -> Dict[str, Any]:
    """Calculate displacement vectors between ground and excited states"""
    print("Calculating displacement vectors...")

    data["displacements"] = {}

    mode_count = data["mode_count"]
    atom_count = data["atom_count"]
    frequencies = data["frequencies"]
    normal_modes = data["normal_modes"]
    geometry = data["geometry"]

    for state in data["requested_states"]:
        if state == 0:
            data["displacements"][state] = {
                mode: 0.0 for mode in range(1, mode_count + 1)
            }
            continue

        if state not in data["forces"]:
            raise ValueError(f"No forces found for state {state}")

        forces = data["forces"][state]

        displacement_vector = calculate_displacement_vector(
            mode_count, atom_count, frequencies, normal_modes, forces, geometry
        )

        data["displacements"][state] = {}
        for mode_idx, displacement in enumerate(displacement_vector):
            mode_number = mode_idx + 1
            data["displacements"][state][mode_number] = displacement

    return data


def calculate_frequency_matrix(
    mode_count: int, frequencies: Dict[int, float]
) -> np.ndarray:
    """Calculate diagonal matrix containing squares of vibrational frequencies"""
    frequency_squares_list: List[float] = []
    for mode in range(1, mode_count + 1):
        frequency_squares_list.append((frequencies[mode] ** 2))
    frequency_squares = np.array(frequency_squares_list)
    frequency_matrix = np.diag(frequency_squares)
    return frequency_matrix


def calculate_inverse_frequency_matrix(
    mode_count: int, frequencies: Dict[int, float]
) -> np.ndarray:
    """Calculate inverse of the frequency eigenvalue matrix"""
    inverse_frequency_list = []
    frequency_matrix = calculate_frequency_matrix(mode_count, frequencies)
    diagonal_frequencies = np.diag(frequency_matrix)

    for j in range(mode_count):
        inverse_frequency_list.append(1.0 / diagonal_frequencies[j])

    inverse_frequency_matrix = np.diag(inverse_frequency_list)
    return inverse_frequency_matrix


def calculate_normal_mode_matrix(
    normal_modes: Dict[int, Dict[str, List[float]]],
    mode_count: int,
    atom_count: int,
) -> np.ndarray:
    """Calculate the normalized mode matrix"""
    normal_mode_matrix = np.zeros((mode_count, 3 * atom_count))

    for mode in range(1, mode_count + 1):
        if mode in normal_modes:
            mode_data = normal_modes[mode]
            x_coords = mode_data.get("x", [0.0] * atom_count)
            y_coords = mode_data.get("y", [0.0] * atom_count)
            z_coords = mode_data.get("z", [0.0] * atom_count)

            flattened = np.zeros(3 * atom_count)
            for atom_idx in range(atom_count):
                flattened[3 * atom_idx] = x_coords[atom_idx]
                flattened[3 * atom_idx + 1] = y_coords[atom_idx]
                flattened[3 * atom_idx + 2] = z_coords[atom_idx]

            flattened = flattened / np.linalg.norm(flattened)

            normal_mode_matrix[mode - 1, :] = flattened
        else:
            raise ValueError(f"Normal mode {mode} not found in normal_modes")

    return normal_mode_matrix


def calculate_gradient_vector(
    forces: np.ndarray,
    atom_count: int,
    geometry: Dict[int, Dict[str, Any]],
) -> np.ndarray:
    """Calculate the mass-weighed gradient vector from forces"""
    gradient_vector = []

    for atom_idx in range(atom_count):
        element = geometry[atom_idx + 1]["element"]
        mass = ATOMIC_MASSES[element]
        mass_factor = np.sqrt(mass / E_MASS)

        fx, fy, fz = forces[atom_idx]
        gradient_vector.extend(
            [-fx / mass_factor, -fy / mass_factor, -fz / mass_factor]
        )

    return np.array(gradient_vector)


def calculate_displacement_vector(
    mode_count: int,
    atom_count: int,
    frequencies: Dict[int, float],
    normal_modes: Dict[int, Dict[str, List[float]]],
    forces: np.ndarray,
    geometry: Dict[int, Dict[str, Any]],
) -> np.ndarray:
    """Calculate displacement vector for excited state"""
    displacements = []

    inverse_frequency_matrix = calculate_inverse_frequency_matrix(
        mode_count, frequencies
    )
    normal_mode_matrix = calculate_normal_mode_matrix(
        normal_modes, mode_count, atom_count
    )
    gradient_vector = calculate_gradient_vector(forces, atom_count, geometry)

    intermediate = np.dot(inverse_frequency_matrix, normal_mode_matrix)
    unscaled_displacements = np.dot(intermediate, gradient_vector)
    unscaled_displacements_list = list(unscaled_displacements)

    for mode_index in range(mode_count):
        frequency_squared = frequencies[mode_index + 1] ** 2
        scaling_factor = frequency_squared ** (1 / 4)
        displacements.append(scaling_factor * unscaled_displacements_list[mode_index])

    return np.array(displacements)


def calculate_physical_parameters(
    data: Dict[str, Any], config: Dict[str, Any]
) -> Dict[str, Any]:
    """Calculate alpha, gamma, adiabatic energies, etc."""
    print("Calculating physical parameters...")

    vertical_energies = []
    for state in sorted(data["oscillator_strengths"].keys()):
        vertical_energies.append(
            data["oscillator_strengths"][state]["excitation_energy"]
        )
    data["vertical_energies"] = vertical_energies

    data["alphas"] = calculate_alpha_parameter(
        data["requested_states"],
        data["mode_count"],
        data["displacements"],
        data["frequencies"],
    )

    data["gammas"] = calculate_gamma_parameter(
        data["requested_states"],
        data["mode_count"],
        data["displacements"],
        data["frequencies"],
    )

    data["adiabatic_energies"] = calculate_adiabatic_energies(
        data["requested_states"],
        data["displacements"],
        data["mode_count"],
        data["frequencies"],
        data["vertical_energies"],
    )

    return data


def calculate_spectrum(data: Dict[str, Any], config: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate spectrum using the chosen method"""
    method = config["method"].lower()
    spectrum_type = config["spectrum_type"].lower()

    print(f"Calculating {spectrum_type} spectrum using {method.upper()} method...")

    energy_min_au = config["energy_min"] * EV_TO_AU
    energy_max_au = config["energy_max"] * EV_TO_AU
    energies_au = np.linspace(energy_min_au, energy_max_au, config["energy_points"])
    energies_ev = energies_au * AU_TO_EV

    integration_params = get_integration_parameters(config, method)

    combined_intensities = np.zeros(len(energies_au))
    individual_intensities = []

    for state_idx, state_number in enumerate(data["requested_states"], 1):
        print(
            f"  Processing state {state_number} ({state_idx}/{len(data['requested_states'])})..."
        )
        state_spectrum = np.zeros(len(energies_au))

        for energy_idx, energy_au in enumerate(energies_au):
            intensity = calculate_spectrum_point(
                energy_au,
                state_idx,
                spectrum_type,
                method,
                data,
                config,
                integration_params,
            )
            state_spectrum[energy_idx] = intensity
            combined_intensities[energy_idx] += intensity

            if (
                energy_idx == 0
                or energy_idx == len(energies_au) - 1
                or (len(energies_au) > 1 and energy_idx == len(energies_au) // 2)
            ):
                percentage = (energy_idx + 1) / len(energies_au) * 100
                print(
                    f"    Energy point {energy_idx + 1}/{len(energies_au)} "
                    f"({percentage:.1f}%)"
                )

        individual_intensities.append(state_spectrum)

    return {
        "energies_ev": energies_ev,
        "intensities": combined_intensities,
        "individual_intensities": individual_intensities,
        "method": method,
        "spectrum_type": spectrum_type,
        "requested_states": data["requested_states"],
    }


def calculate_spectrum_point(
    energy_au: float,
    state_idx: int,
    spectrum_type: str,
    method: str,
    data: Dict[str, Any],
    config: Dict[str, Any],
    integration_params: Dict[str, Any],
) -> float:
    """Calculate spectrum intensity at a single energy point"""

    if method == "lq2":
        return calculate_lq2_spectrum_point(
            energy_au,
            state_idx,
            spectrum_type,
            data["oscillator_strengths"],
            data["displacements"],
            data["frequencies"],
            data["mode_count"],
            data["alphas"],
            data["vertical_energies"],
            config["stokes_shift"],
            data["adiabatic_energies"],
            data["requested_states"],
        )

    elif method == "lq3":
        return calculate_lq3_spectrum_point(
            energy_au,
            state_idx,
            spectrum_type,
            data["oscillator_strengths"],
            data["displacements"],
            data["frequencies"],
            data["mode_count"],
            data["alphas"],
            data["gammas"],
            data["vertical_energies"],
            config["stokes_shift"],
            data["adiabatic_energies"],
            data["requested_states"],
            integration_params,
        )

    elif method == "imdho":
        return calculate_imdho_spectrum_point(
            energy_au,
            state_idx,
            spectrum_type,
            data["oscillator_strengths"],
            data["displacements"],
            data["frequencies"],
            data["mode_count"],
            config["gamma_broadening"],
            config["theta_broadening"],
            config["temperature"],
            data["vertical_energies"],
            config["stokes_shift"],
            data["adiabatic_energies"],
            data["requested_states"],
            integration_params,
        )

    else:
        raise ValueError(f"Unknown calculation method: {method}")


def get_integration_parameters(config: Dict[str, Any], method: str) -> Dict[str, Any]:
    """Get integration parameters for the chosen method"""
    if method == "lq3":
        return {
            "max_time_slices": config.get("max_time_slices", 5000),
            "time_step": config.get("lq3_time_step", 30),
            "convergence": config.get("lq3_convergence", 0.0000001),
        }
    elif method == "imdho":
        return {
            "max_time_slices": config.get("max_time_slices", 5000),
            "time_step": config.get("imdho_time_step", 30),
            "convergence": config.get("imdho_convergence", 0.000000001),
        }
    else:
        return {}


def load_configuration(config_file: str) -> Dict[str, Any]:
    """Load configuration from TOML file"""
    if config_file.endswith(".toml"):

        config_dir = os.path.dirname(os.path.abspath(config_file))
        with open(config_file, "r", encoding="utf-8") as f:
            content = f.read()

        # Try importing toml from various places.
        try:
            import tomllib  # not available before Python 3.11
        except ImportError:
            try:
                import pip._vendor.tomli as tomllib  # type: ignore
            except ImportError:
                try:
                    import pip._vendor.toml as tomllib  # type: ignore
                except ImportError:
                    import toml as tomllib  # type: ignore

        config = tomllib.loads(content)

    else:
        raise ValueError(
            "Only TOML configuration files are supported. Use .toml extension"
        )

    flattened = flatten_config(config)
    flattened["config_dir"] = config_dir

    return flattened


def flatten_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """Flatten the configuration file into a simple dictionary"""
    flattened = {}

    sections = ["files", "calculation", "imdho", "output", "integration"]
    for section in sections:
        if section in config:
            flattened.update(config[section])

    flattened.setdefault("print_individual_states", False)
    flattened.setdefault("gradient_files", [])
    flattened.setdefault("stokes_shift", 0.0)
    flattened.setdefault("temperature", 298.15)
    flattened.setdefault("gamma_broadening", 0.01)
    flattened["gamma_broadening"] *= EV_TO_AU
    flattened.setdefault("theta_broadening", 0.0)
    flattened["theta_broadening"] *= EV_TO_AU

    return flattened


def main() -> None:
    """Main entry point"""
    if len(sys.argv) != 2:
        print("Usage: python main.py <config_file.toml>")
        print("Example: python main.py calculation_setup.toml")
        sys.exit(1)

    config = load_configuration(sys.argv[1])

    required_params = [
        "vibrations_filename",
        "output_filename",
        "states",
        "energy_min",
        "energy_max",
        "energy_points",
        "method",
        "spectrum_type",
    ]

    missing_params = [param for param in required_params if param not in config]
    if missing_params:
        raise ValueError(f"Missing required configuration parameters: {missing_params}")

    run_spectrum_calculation(config)


if __name__ == "__main__":
    main()
