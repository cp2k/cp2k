# SPDX-License-Identifier: GPL-2.0-or-later

from concurrent.futures import ThreadPoolExecutor
import ctypes as ct
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pytest

from cp2k import CP2K
import cp2k.library as library


def test_abi_signatures(runtime, fake_library):
    assert runtime.version == "CP2K test library"
    assert fake_library.cp2k_get_positions.argtypes == [ct.c_int, library._DP, ct.c_int]
    assert fake_library.cp2k_get_positions.restype is None


def test_missing_library(monkeypatch):
    monkeypatch.delenv("CP2K_LIBRARY", raising=False)
    monkeypatch.setattr(library, "find_library", lambda _: None)
    with pytest.raises(FileNotFoundError, match="CP2K_LIBRARY"):
        library._load_library(None)


def test_missing_symbol(fake_library):
    del fake_library.cp2k_get_nparticle
    with pytest.raises(RuntimeError, match="get_nparticle"):
        CP2K(library="test-library")
    assert not fake_library.calls


def test_environment_arrays_and_results(runtime, fake_library, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with runtime.create_force_env("&GLOBAL\n&END\n") as env:
        assert env.natom == env.nparticle == 2
        with pytest.raises(RuntimeError, match="calculate"):
            _ = env.potential_energy
        np.testing.assert_array_equal(env.cell, fake_library.cell)
        positions = env.positions
        positions[0, 0] = -30
        assert env.positions[0, 0] == 0
        # Fortran-contiguous, noncontiguous, and float32 inputs are copied.
        env.positions = np.asfortranarray(np.ones((2, 3), dtype=np.float32))
        env.cell = [[4, 0, 0], [0.2, 5, 0], [0.1, 0.3, 6]]
        env.set_velocities(np.zeros((2, 6))[:, ::2])
        assert fake_library.cell[1, 0] == 0.2
        result = env.calculate()
        assert result.energy == -1.0
        np.testing.assert_array_equal(result.forces, np.ones((2, 3)))
        env.positions = env.positions + 0.1
        with pytest.raises(RuntimeError, match="calculate"):
            _ = env.forces
        assert env.calculate(forces=False).forces is None
        with pytest.raises(RuntimeError, match="forces=True"):
            _ = env.forces
    assert not list(tmp_path.glob("cp2k-python-*.inp"))
    env.close()
    with pytest.raises(RuntimeError, match="closed"):
        _ = env.positions


@pytest.mark.parametrize(
    "bad",
    [np.zeros(6), np.zeros((3, 3)), [[0, 0, 0], [0, 0, np.nan]], np.ones((2, 3)) * 1j],
)
def test_bad_positions(runtime, bad, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with runtime.create_force_env("&GLOBAL\n&END\n") as env:
        with pytest.raises(ValueError):
            env.positions = bad


def test_bad_cell(runtime, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with runtime.create_force_env("&GLOBAL\n&END\n") as env:
        for cell in (np.zeros((3, 3)), -np.eye(3)):
            with pytest.raises(ValueError, match="right-handed"):
                env.cell = cell


def test_file_validation(runtime, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for kwargs in (
        {},
        {"inp": " "},
        {"inp": "x\0"},
        {"input_file": "absent"},
        {"inp": "x", "input_file": "x"},
    ):
        with pytest.raises((ValueError, FileNotFoundError)):
            runtime.create_force_env(**kwargs)
    with pytest.raises(IsADirectoryError):
        runtime.create_force_env("x", output_file=tmp_path)
    with pytest.raises(ValueError):
        library._path_bytes("x\0")
    source = tmp_path / "input.inp"
    source.write_text("&GLOBAL\n&END\n")
    with runtime.create_force_env(input_file=source):
        with pytest.raises(RuntimeError, match="Close all force environments"):
            runtime.run_input(input_file=source)
    runtime.run_input(input_file=source)
    assert source.exists()


def test_cleanup_and_single_initialization(
    runtime, fake_library, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    env1 = runtime.create_force_env("x")
    with pytest.raises(RuntimeError, match="current force environment"):
        runtime.create_force_env("x")
    env1.close()
    env2 = runtime.create_force_env("x")
    with pytest.raises(RuntimeError, match="once"):
        CP2K(library="test-library")
    runtime.close()
    runtime.close()
    with pytest.raises(RuntimeError, match="once"):
        CP2K(library="test-library")
    assert fake_library.calls[-2:] == [
        ("cp2k_destroy_force_env", (env2._handle,)),
        ("cp2k_finalize", ()),
    ]


def test_threads_and_forks_rejected(runtime, monkeypatch):
    with ThreadPoolExecutor(1) as pool:
        with pytest.raises(RuntimeError, match="main thread"):
            pool.submit(lambda: runtime.version).result()
    with monkeypatch.context() as patch:
        patch.setattr(library.os, "getpid", lambda: -1)
        with pytest.raises(RuntimeError, match="fork"):
            _ = runtime.version


def test_nonfinite_results(runtime, fake_library, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with runtime.create_force_env("x") as env:
        fake_library.forces[0, 0] = np.inf
        with pytest.raises(RuntimeError, match="non-finite"):
            env.calculate()
        with pytest.raises(RuntimeError, match="calculate"):
            _ = env.potential_energy


def test_external_mpi_ownership(fake_library, monkeypatch, tmp_path):
    class Comm:
        def py2f(self):
            return 42

    comm = Comm()
    mpi = SimpleNamespace(
        Is_initialized=lambda: True,
        Is_finalized=lambda: False,
        Query_thread=lambda: 3,
        THREAD_MULTIPLE=3,
        Intracomm=Comm,
        COMM_NULL=None,
        COMM_WORLD=comm,
    )
    monkeypatch.setitem(sys.modules, "mpi4py.MPI", mpi)
    monkeypatch.chdir(tmp_path)
    with CP2K(library="test-library") as session:
        with session.create_force_env("x"):
            pass
        session.run_input("x")
    assert fake_library.calls[0] == ("cp2k_init_without_mpi_comm", (42,))
    assert fake_library.calls[-1] == ("cp2k_finalize_without_mpi", ())
    assert any(name == "cp2k_create_force_env_comm" for name, _ in fake_library.calls)
    assert any(name == "cp2k_run_input_comm" for name, _ in fake_library.calls)


def test_insufficient_mpi_threads(fake_library, monkeypatch):
    mpi = SimpleNamespace(
        Is_initialized=lambda: True,
        Is_finalized=lambda: False,
        Query_thread=lambda: 1,
        THREAD_MULTIPLE=3,
    )
    monkeypatch.setitem(sys.modules, "mpi4py.MPI", mpi)
    with pytest.raises(RuntimeError, match="THREAD_MULTIPLE"):
        CP2K(library="test-library")
    assert not fake_library.calls
