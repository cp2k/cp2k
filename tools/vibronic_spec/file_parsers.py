"""
Parsers
"""

from typing import Dict, List, Any, Optional
import numpy as np
import re
from constants import EV_TO_AU, WAVENUMBER_TO_AU


def parse_molden_file(file_path: str) -> Dict[str, Any]:
    """
    Parse
    - geometry from [Atoms] section
    - normal modes from [FR-NORM-COORD] section
    - frequencies from [FREQ] section
    of VIBRATIONS Molden file

    Returns dict with 'geometry', 'normal_modes', 'mode_count', 'frequencies',
    'atom_count' and 'negative_freq_warnings'
    """
    data: Dict[str, Any] = {}
    geometry: Dict[int, Dict[str, Any]] = {}
    normal_modes: Dict[int, Dict[str, List[float]]] = {}
    frequencies: Dict[int, float] = {}
    negative_freq_warnings: List[float] = []

    with open(file_path, "rt") as file:
        content = file.read()

    # Parse geometry from [Atoms] section
    if "[Atoms]" in content:
        atoms_section = content.split("[Atoms]")[1]
        if "[" in atoms_section:
            atoms_section = atoms_section.split("[")[0]

        lines = atoms_section.strip().split("\n")

        atom_idx = 1
        for line in lines:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) >= 5:
                try:
                    element = parts[0]
                    x, y, z = map(float, parts[3:6])

                    geometry[atom_idx] = {"element": element, "coordinates": [x, y, z]}
                    atom_idx += 1

                except (ValueError, IndexError):
                    continue

    # Parse frequencies from [FREQ] section
    if "[FREQ]" in content:
        freq_section = content.split("[FREQ]")[1].split("[FR-COORD]")[0]
        mode_count = 0

        for line in freq_section.strip().split("\n"):
            line = line.strip()
            if line:
                try:
                    frequency = float(line)
                    # Only keep positive frequencies
                    if frequency > 0:
                        mode_count += 1
                        frequencies[mode_count] = frequency * WAVENUMBER_TO_AU
                    else:
                        negative_freq_warnings.append(frequency)
                except ValueError:
                    continue

    # Parse normal modes from [FR-NORM-COORD] section
    if "[FR-NORM-COORD]" in content:
        norm_section = content.split("[FR-NORM-COORD]")[1].split("[INT]")[0]
        lines = norm_section.strip().split("\n")

        current_mode: Optional[int] = None
        collected_coords: List[float] = []

        for line in lines:
            line = line.strip()
            if not line:
                continue

            if "vibration" in line.lower():
                # Save previous mode
                if current_mode is not None and collected_coords:
                    atom_count = len(geometry)
                    x_coords = collected_coords[0::3][:atom_count]
                    y_coords = collected_coords[1::3][:atom_count]
                    z_coords = collected_coords[2::3][:atom_count]

                    normal_modes[current_mode] = {
                        "x": x_coords,
                        "y": y_coords,
                        "z": z_coords,
                    }

                # New mode
                try:
                    current_mode = int(line.split()[1])
                    collected_coords = []
                except (ValueError, IndexError):
                    current_mode = None

            elif current_mode is not None:
                try:
                    coords = list(map(float, line.split()[:3]))
                    collected_coords.extend(coords)
                except ValueError:
                    continue

        # Save last mode
        if current_mode is not None and collected_coords:
            atom_count = len(geometry)
            x_coords = collected_coords[0::3][:atom_count]
            y_coords = collected_coords[1::3][:atom_count]
            z_coords = collected_coords[2::3][:atom_count]
            normal_modes[current_mode] = {"x": x_coords, "y": y_coords, "z": z_coords}

    data["geometry"] = geometry
    data["normal_modes"] = normal_modes
    data["mode_count"] = len(frequencies)
    data["frequencies"] = frequencies
    data["atom_count"] = len(geometry)
    data["negative_freq_warnings"] = negative_freq_warnings

    return data


def parse_excited_state_forces(file_path: str) -> Dict[int, Dict[str, Any]]:
    """Parse
    - excitation energies
    - oscillator strengths
    - excited state forces
     from CP2K TDFORCE file"""
    data: Dict[int, Dict[str, Any]] = {}

    with open(file_path, "r") as f:
        content = f.read()

    state_blocks = content.split("# STATE NR.")[1:]

    for block in state_blocks:
        lines = block.strip().split("\n")

        if not lines:
            continue

        first_line = lines[0].strip()
        parts = first_line.split()

        if len(parts) < 6:
            continue

        state_number = int(parts[0])
        energy_ev = float(parts[1])
        excitation_energy = energy_ev * EV_TO_AU

        osc_str = None
        for i, part in enumerate(parts):
            if part == "strength:" and i + 1 < len(parts):
                osc_str = float(parts[i + 1])
                break

        if osc_str is None:
            continue

        data[state_number] = {
            "oscillator_strength": osc_str,
            "excitation_energy": excitation_energy,
        }

        if len(lines) >= 3:
            try:
                atom_count = int(lines[1].strip())
                if len(lines) >= 3 + atom_count:
                    force_array = np.zeros((atom_count, 3))

                    for i in range(atom_count):
                        force_line = lines[3 + i].strip()
                        fx, fy, fz = map(float, force_line.split())
                        force_array[i] = [fx, fy, fz]

                    data[state_number]["force"] = force_array

            except (ValueError, IndexError):
                # state doesn't have forces
                continue

    print(f"INFO: Parsed {len(data)} total states")

    states_with_forces = [s for s in data if "force" in data[s]]
    print(f"INFO: Found {len(states_with_forces)} states with forces")

    return data
