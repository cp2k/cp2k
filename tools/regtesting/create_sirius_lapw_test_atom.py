#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Write small LAPW integration fixtures, not production atom setups."""

import argparse
import json
import math
from pathlib import Path


def helium_setup():
    radius = [1e-7 + (20.0 - 1e-7) * (i / 2999) ** 3 for i in range(3000)]
    basis = [{"enu": -0.3, "dme": dme, "auto": 0} for dme in (0, 1)]
    return {
        "name": "Helium integration fixture",
        "symbol": "He",
        "number": 2,
        "mass": 4.002602,
        "rmin": 1e-6,
        "rmt": 1.8,
        "nrmt": 801,
        "core": "",
        "valence": [{"basis": basis}]
        + [{"n": l + 1, "l": l, "basis": basis} for l in range(4)],
        "lo": [
            {
                "l": 0,
                "basis": [
                    {"n": 1, "enu": -0.7, "dme": dme, "auto": 0} for dme in (0, 1)
                ],
            }
        ],
        "free_atom": {
            "radial_grid": radius,
            "density": [16.0 / math.pi * math.exp(-4.0 * r) for r in radius],
        },
    }


def hydrogen_setup():
    atom = helium_setup()
    atom.update(name="Hydrogen integration fixture", symbol="H", number=1, mass=1.00794)
    atom["free_atom"]["density"] = [
        math.exp(-2.0 * r) / math.pi for r in atom["free_atom"]["radial_grid"]
    ]
    return atom


def lithium_setup():
    atom = helium_setup()
    atom.update(
        name="Lithium integration fixture with an explicit 1s core",
        symbol="Li",
        number=3,
        mass=6.94,
        rmt=2.5,
        nrmt=1601,
        core="1s",
    )
    for channel in atom["valence"]:
        if "n" in channel:
            channel["n"] = max(2, channel["n"])
    for orbital in atom["lo"][0]["basis"]:
        orbital.update(n=2, enu=-0.2)
    atom["free_atom"]["density"] = [
        2 * 2.7**3 / math.pi * math.exp(-5.4 * r)
        + 0.7**3 / math.pi * math.exp(-1.4 * r)
        for r in atom["free_atom"]["radial_grid"]
    ]
    return atom


def neon_setup():
    atom = helium_setup()
    atom.update(
        name="Neon integration fixture with an explicit 1s core",
        symbol="Ne",
        number=10,
        mass=20.1797,
        nrmt=3201,
        core="1s",
    )
    for channel in atom["valence"]:
        if "n" in channel:
            channel["n"] = max(2, channel["n"])
    for orbital in atom["lo"][0]["basis"]:
        orbital["n"] = 2
    # Local 2s/2p orbitals resolve the valence shell at the small test PW cutoff.
    atom["lo"].extend(
        [
            {
                "l": 0,
                "basis": [
                    {"n": 2, "enu": energy, "dme": 0, "auto": 0}
                    for energy in (-2.0, -1.0, -0.3)
                ],
            },
            {
                "l": 1,
                "basis": [
                    {"n": 2, "enu": -1.0, "dme": dme, "auto": 0} for dme in (0, 1)
                ],
            },
        ]
    )
    atom["free_atom"]["density"] = [
        2 * 9.0**3 / math.pi * math.exp(-18 * r)
        + 8 * 1.6**3 / math.pi * math.exp(-3.2 * r)
        for r in atom["free_atom"]["radial_grid"]
    ]
    return atom


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    setups = {
        "H": hydrogen_setup,
        "He": helium_setup,
        "Li": lithium_setup,
        "Ne": neon_setup,
    }
    parser.add_argument("--element", choices=tuple(setups), default="He")
    args = parser.parse_args()
    atom = setups[args.element]()
    args.output.write_text(json.dumps(atom, indent=2) + "\n")
