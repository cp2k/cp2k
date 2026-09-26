"""Cheap input shared by the standalone PLUMED coupling tests."""

# SPDX-License-Identifier: GPL-2.0-or-later


def plumed_input():
    """Cheap, analytic, periodic test potential with nonzero off-diagonal stress."""
    return {
        "GLOBAL": {"PROJECT": "argon", "PRINT_LEVEL": "LOW"},
        "FORCE_EVAL": {
            "METHOD": "FIST",
            "STRESS_TENSOR": "ANALYTICAL",
            "MM": {
                "FORCEFIELD": {
                    "CHARGE": {"ATOM": "Ar", "CHARGE": 0},
                    "NONBONDED": {
                        "LENNARD-JONES": {
                            "ATOMS": ["Ar", "Ar"],
                            "EPSILON": "[hartree] 0.001",
                            "SIGMA": 3,
                            "RCUT": 8,
                        }
                    },
                },
                "POISSON": {"EWALD": {"EWALD_TYPE": "NONE"}},
            },
            "SUBSYS": {
                "CELL": {"A": [20, 0, 0], "B": [1, 21, 0], "C": [2, 3, 22]},
                "COORD": {"_lines": ["Ar 4 4 4", "Ar 7.2 4.8 4.5"]},
                "KIND": {"_": "Ar", "ELEMENT": "Ar"},
                "TOPOLOGY": {"CONN_FILE_FORMAT": "OFF"},
            },
        },
    }
