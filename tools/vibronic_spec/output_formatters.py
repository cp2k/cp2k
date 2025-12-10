"""
Functions for output formatting
"""

import os
from typing import Any, Dict, List


def write_spectrum_output(results: Dict[str, Any], config: Dict[str, Any]) -> None:
    """
    Main output function
    Writes combined spectrum to main file and individual states
    to separate file if requested
    """
    output_file = config["output_filename"]

    write_combined_spectrum(results, output_file)

    if config.get("print_individual_states", False):
        write_individual_states_file(results, config)


def write_combined_spectrum(results: Dict[str, Any], output_file: str) -> None:
    """
    Write combined spectrum (sum of all states) to main output file
    Format: energy(eV) intensity(au)
    """
    with open(output_file, "w") as f:
        f.write("# energy(eV) intensity(au)\n")

        for energy, intensity in zip(results["energies_ev"], results["intensities"]):
            f.write(f"{energy:.6f} {intensity:.6e}\n")


def write_individual_states_file(
    results: Dict[str, Any], config: Dict[str, Any]
) -> None:
    """
    Write individual state spectra to separate file in CP2K-like format
    File: states_{method}.txt
    Format:
        # STATE NR.   1
        energy1 intensity1
        energy2 intensity2
        ...
        # STATE NR.   2
        energy1 intensity1
        energy2 intensity2
        ...
    """
    method = results["method"].lower()
    output_dir = os.path.dirname(config["output_filename"]) or "."
    states_file = f"{output_dir}/states_{method}.txt"

    with open(states_file, "w") as f:
        energies_ev = results["energies_ev"]
        individual_intensities = results["individual_intensities"]

        for state_idx, state_number in enumerate(results["requested_states"]):
            f.write(f"# STATE NR.   {state_number}\n")

            state_intensities = individual_intensities[state_idx]
            for energy, intensity in zip(energies_ev, state_intensities):
                f.write(f"{energy:.6f} {intensity:.6e}\n")

            if state_idx < len(results["requested_states"]) - 1:
                f.write("\n")
