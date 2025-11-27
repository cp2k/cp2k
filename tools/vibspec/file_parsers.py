"""
Parsers
"""

from typing import Dict, List, Any
import numpy as np
import re
from constants import EV_TO_AU, WAVENUMBER_TO_AU, ANGSTROM_TO_BOHR


def parse_vibrational_frequencies(file_path: str) -> Dict[int, float]:
    """
    Extract vibrational frequencies from Molden file
    """
    data = {}
    mode_count = 0
    warnings = []
    
    with open(file_path, "rt") as file:
        content = file.read()
        
        if "[FREQ]" in content:
            freq_section = content.split("[FREQ]")[1].split("[FR-COORD]")[0]
            
            for line in freq_section.strip().split('\n'):
                line = line.strip()
                if line:
                    try:
                        frequency = float(line)
                        # Only keep positive frequencies
                        if frequency > 0:
                            mode_count += 1
                            data[mode_count] = frequency * WAVENUMBER_TO_AU
                        else:
                            warnings.append(frequency)
                    except ValueError:
                        continue
    
    return data, warnings


def parse_normal_modes(file_path: str, atom_count: int) -> Dict[int, Dict[str, List[float]]]:
    """
    Extract normal mode vectors from Molden format file
    Format: [FR-NORM-COORD] section with vibration headers
    """
    data = {}
    
    with open(file_path, "rt") as file:
        content = file.read()
        
        if "[FR-NORM-COORD]" in content:
            norm_section = content.split("[FR-NORM-COORD]")[1].split("[INT]")[0]
            lines = norm_section.strip().split('\n')
            
            current_mode = None
            collected_coords = []
            
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                    
                if "vibration" in line.lower():
                    # Save previous mode
                    if current_mode is not None and collected_coords:
                        x_coords = collected_coords[0::3][:atom_count]
                        y_coords = collected_coords[1::3][:atom_count]
                        z_coords = collected_coords[2::3][:atom_count]
                        
                        data[current_mode] = {"x": x_coords, "y": y_coords, "z": z_coords}
                    
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
                x_coords = collected_coords[0::3][:atom_count]
                y_coords = collected_coords[1::3][:atom_count]
                z_coords = collected_coords[2::3][:atom_count]
                data[current_mode] = {"x": x_coords, "y": y_coords, "z": z_coords}
    
    return data


def parse_excited_state_forces(file_path):
    """Parse 
    - excitation energies
    - oscillator strengths
    - excited state forces
     from CP2K TDFORCE file"""
    data = {}
    
    with open(file_path, "r") as f:
        content = f.read()
    
    state_blocks = content.split("# STATE NR.")[1:]
    
    for block in state_blocks:
        lines = block.strip().split('\n')
        
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
            if part == "strength:" and i+1 < len(parts):
                osc_str = float(parts[i+1])
                break
        
        if osc_str is None:
            continue
            
        data[state_number] = {
            "oscillator_strength": osc_str,
            "excitation_energy": excitation_energy
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

def parse_geometry_from_xyz(file_path: str) -> Dict[int, Dict[str, Any]]:
    """
    Extract molecular geometry from standard XYZ format file
    """
    data = {}
    
    with open(file_path, "rt") as file:
        lines = file.readlines()
        
        if len(lines) >= 2:
            try:
                atom_count = int(lines[0].strip())
                
                for i in range(2, 2 + atom_count):
                    if i < len(lines):
                        parts = lines[i].split()
                        if len(parts) >= 4:
                            element = parts[0]
                            x, y, z = map(float, parts[1:4])
                            
                            atom_number = i - 1
                            data[atom_number] = {
                                "element": element,
                                "coordinates": [x * ANGSTROM_TO_BOHR, y * ANGSTROM_TO_BOHR, z * ANGSTROM_TO_BOHR]
                            }
            except ValueError:
                print(f"Invalid XYZ format in {file_path}")
    
    return data


def get_atom_count(geometry_data: Dict[int, Any]) -> int:
    """Get number of atoms from geometry data"""
    return len(geometry_data)