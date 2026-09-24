# SPDX-License-Identifier: GPL-2.0-or-later

import ctypes as ct
import os
from types import SimpleNamespace

import numpy as np
import pytest

import cp2k.library as library


@pytest.fixture(scope="session")
def real_runtime():
    path = os.environ.get("CP2K_TEST_LIBRARY")
    if not path:
        pytest.skip("Set CP2K_TEST_LIBRARY to run native integration tests")
    comm = None
    if os.environ.get("CP2K_TEST_EXTERNAL_MPI"):
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
    with library.CP2K(library=path, comm=comm) as runtime:
        yield runtime


class FakeFunction:
    def __init__(self, name, lib):
        self.name, self.lib = name, lib

    def __call__(self, *args):
        self.lib.calls.append((self.name, args))
        operation = self.name.removeprefix("cp2k_")
        if operation == "get_version":
            args[0].value = b"CP2K test library"
        elif operation.startswith("create_force_env"):
            self.lib.counter += 1
            args[0]._obj.value = self.lib.counter
        elif operation in ("get_natom", "get_nparticle"):
            args[1]._obj.value = 2
        elif operation == "get_potential_energy":
            args[1]._obj.value = -1.0
        elif operation == "get_scf_convergence":
            args[1]._obj.value = self.lib.scf_status
        elif operation in ("get_positions", "get_cell", "get_forces"):
            data = getattr(self.lib, operation[4:])
            np.ctypeslib.as_array(args[1], shape=(data.size,))[:] = data.ravel()
        elif operation in ("set_positions", "set_cell", "set_velocities"):
            n = args[2] if len(args) == 3 else 9
            setattr(
                self.lib,
                operation[4:],
                np.ctypeslib.as_array(args[1], shape=(n,)).copy().reshape(-1, 3),
            )


@pytest.fixture
def fake_library(monkeypatch):
    lib = SimpleNamespace(
        calls=[],
        counter=0,
        scf_status=1,
        positions=np.arange(6.0).reshape(2, 3),
        cell=np.array([[5.0, 0, 0], [0.4, 6, 0], [0.2, 0.3, 7]]),
        forces=np.ones((2, 3)),
    )
    names = (
        "get_version init init_without_mpi init_without_mpi_comm finalize "
        "finalize_without_mpi create_force_env create_force_env_comm destroy_force_env "
        "get_natom get_nparticle get_potential_energy get_positions get_cell get_forces "
        "set_positions set_cell set_velocities calc_energy calc_energy_force "
        "run_input run_input_comm get_scf_convergence"
    ).split()
    for name in names:
        setattr(lib, "cp2k_" + name, FakeFunction("cp2k_" + name, lib))
    monkeypatch.setattr(library.ct, "CDLL", lambda _: lib)
    monkeypatch.setattr(library, "_runtime", None)
    monkeypatch.delitem(library.sys.modules, "mpi4py.MPI", raising=False)
    return lib


@pytest.fixture
def runtime(fake_library):
    with library.CP2K(library="test-library") as session:
        yield session


@pytest.fixture
def h2_input():
    return {
        "GLOBAL": {"PROJECT": "python-h2", "PRINT_LEVEL": "SILENT"},
        "FORCE_EVAL": {
            "METHOD": "Quickstep",
            "DFT": {
                "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                "MGRID": {"CUTOFF": 200},
                "QS": {"EPS_DEFAULT": 1e-10},
                "SCF": {
                    "EPS_SCF": 1e-9,
                    "MAX_SCF": 100,
                    "OT": {"MINIMIZER": "DIIS"},
                    "PRINT": {"RESTART": {"_": "OFF"}},
                },
                "XC": {"XC_FUNCTIONAL": {"_": "PADE"}},
            },
            "SUBSYS": {
                "CELL": {"ABC": [8, 8, 8]},
                "COORD": {"_lines": ["H 3.6 4 4", "H 4.4 4 4"]},
                "KIND": {
                    "_": "H",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH",
                    "POTENTIAL": "GTH-PADE-q1",
                },
            },
        },
    }
