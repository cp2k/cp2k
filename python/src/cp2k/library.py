"""Typed ctypes bindings with explicit process, MPI, and environment ownership."""

# SPDX-License-Identifier: GPL-2.0-or-later

import atexit
from collections.abc import Mapping
from contextlib import contextmanager
import ctypes as ct
from ctypes.util import find_library
from dataclasses import dataclass, field
import os
from pathlib import Path
import sys
import tempfile
import threading

import numpy as np

from .input import input_to_string

_runtime = None
_C_INT_MAX = 2 ** (8 * ct.sizeof(ct.c_int) - 1) - 1
_DP = ct.POINTER(ct.c_double)


def _load_library(path):
    name = os.fspath(path) if path is not None else os.environ.get("CP2K_LIBRARY")
    if not name:
        name = find_library("cp2k")
    if not name:
        raise FileNotFoundError(
            "libcp2k was not found. Build CP2K with BUILD_SHARED_LIBS=ON and set "
            "CP2K_LIBRARY to the shared library, or pass library= explicitly."
        )
    try:
        lib = ct.CDLL(name)
    except OSError as error:
        raise OSError(f"Cannot load libcp2k {name!r}: {error}") from error
    signatures = {
        "get_version": [ct.POINTER(ct.c_char), ct.c_int],
        "init": [],
        "init_without_mpi": [],
        "finalize": [],
        "finalize_without_mpi": [],
        "create_force_env": [ct.POINTER(ct.c_int), ct.c_char_p, ct.c_char_p],
        "create_force_env_comm": [
            ct.POINTER(ct.c_int),
            ct.c_char_p,
            ct.c_char_p,
            ct.c_int,
        ],
        "destroy_force_env": [ct.c_int],
        "get_natom": [ct.c_int, ct.POINTER(ct.c_int)],
        "get_nparticle": [ct.c_int, ct.POINTER(ct.c_int)],
        "get_potential_energy": [ct.c_int, _DP],
        "get_positions": [ct.c_int, _DP, ct.c_int],
        "get_forces": [ct.c_int, _DP, ct.c_int],
        "set_positions": [ct.c_int, _DP, ct.c_int],
        "set_velocities": [ct.c_int, _DP, ct.c_int],
        "get_cell": [ct.c_int, _DP],
        "set_cell": [ct.c_int, _DP],
        "calc_energy": [ct.c_int],
        "calc_energy_force": [ct.c_int],
        "run_input": [ct.c_char_p, ct.c_char_p],
        "run_input_comm": [ct.c_char_p, ct.c_char_p, ct.c_int],
    }
    if hasattr(lib, "cp2k_init_without_mpi_comm"):
        signatures["init_without_mpi_comm"] = [ct.c_int]
    if hasattr(lib, "cp2k_get_stress_tensor"):
        signatures["get_stress_tensor"] = [ct.c_int, _DP, ct.POINTER(ct.c_int)]
    if hasattr(lib, "cp2k_get_scf_convergence"):
        signatures["get_scf_convergence"] = [ct.c_int, ct.POINTER(ct.c_int)]
    for function, arguments in signatures.items():
        try:
            symbol = getattr(lib, "cp2k_" + function)
        except AttributeError as error:
            raise RuntimeError(f"libcp2k is missing cp2k_{function}") from error
        symbol.argtypes = arguments
        symbol.restype = None
    return lib


def _path_bytes(path, *, output=False):
    if output and os.fspath(path) == "__STD_OUT__":
        return b"__STD_OUT__"
    path = Path(path).expanduser().absolute()
    encoded = os.fsencode(path)
    if b"\0" in encoded or len(encoded) >= 1024:
        raise ValueError(
            "CP2K paths must contain no NUL and be shorter than 1024 bytes"
        )
    if output:
        if not path.parent.is_dir():
            raise FileNotFoundError(path.parent)
        if path.is_dir():
            raise IsADirectoryError(path)
    elif not path.is_file():
        raise FileNotFoundError(path)
    return encoded


def _array(values, shape, name):
    if np.iscomplexobj(values):
        raise ValueError(f"{name} must be real")
    result = np.array(values, dtype=np.float64, order="C", copy=True)
    if result.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {result.shape}")
    if result.size > _C_INT_MAX or not np.isfinite(result).all():
        raise ValueError(f"{name} must have a C-int-sized array of finite values")
    return result


@dataclass(frozen=True)
class CalculationResult:
    """Atomic-unit results; stress/virial are potential-only, pressure-positive.

    Stress has units hartree/bohr**3; virial = stress * volume has units hartree.
    Both are Cartesian (3, 3) tensors, with no kinetic contribution.
    SCF status is None when unavailable.
    """

    energy: float
    forces: np.ndarray | None
    stress: np.ndarray | None = None
    virial: np.ndarray | None = None
    scf_converged: bool | None = field(default=None, kw_only=True)


class SCFConvergenceError(RuntimeError):
    """CP2K returned from an SCF that did not meet its convergence criteria."""


class CP2K:
    """One libcp2k runtime per Python process, with one live environment at a time.

    ``library`` is the shared library path (default: CP2K_LIBRARY, then system
    loader). ``comm`` may be a live mpi4py intracommunicator; all ranks must
    execute the same calls. MPI must provide THREAD_MULTIPLE and use the same
    implementation as CP2K. If mpi4py has already initialized MPI, it retains
    ownership even when comm is omitted. Otherwise CP2K manages MPI.

    Use a context manager or close() explicitly. Closing is terminal: native
    CP2K cannot safely be reinitialized. Keep this runtime open in notebooks,
    and create/close force environments as needed. Only use the main thread.
    Native CP2K errors may abort Python; they are not Python exceptions.
    """

    def __init__(self, library=None, *, comm=None):
        global _runtime
        if threading.current_thread() is not threading.main_thread():
            raise RuntimeError("CP2K must be used from Python's main thread")
        if _runtime is not None:
            raise RuntimeError(
                "CP2K may be initialized only once per process. Reuse the existing "
                "runtime, or start a fresh Python process after close()/fork()."
            )
        self._lib = _load_library(library)
        self._pid = os.getpid()
        self._closed = False
        self._environments = []
        self._comm = None
        self._mpi = sys.modules.get("mpi4py.MPI")
        self._external_mpi = False
        if comm is not None:
            from mpi4py import MPI

            self._mpi = MPI
        if self._mpi is not None:
            mpi = self._mpi
            if mpi.Is_finalized():
                raise RuntimeError("MPI has already been finalized")
            if comm is not None and not mpi.Is_initialized():
                raise RuntimeError("Initialize MPI before passing a communicator")
            if mpi.Is_initialized():
                if mpi.Query_thread() < mpi.THREAD_MULTIPLE:
                    raise RuntimeError("CP2K requires MPI_THREAD_MULTIPLE")
                comm = mpi.COMM_WORLD if comm is None else comm
                if not isinstance(comm, mpi.Intracomm) or comm == mpi.COMM_NULL:
                    raise ValueError("comm must be a live MPI intracommunicator")
                handle = comm.py2f()
                if not -_C_INT_MAX - 1 <= handle <= _C_INT_MAX:
                    raise OverflowError("MPI Fortran communicator does not fit a C int")
                if not hasattr(self._lib, "cp2k_init_without_mpi_comm"):
                    raise RuntimeError(
                        "This libcp2k lacks communicator-aware initialization"
                    )
                # Keep the communicator alive; the caller must not Free it early.
                self._comm = comm
                self._external_mpi = True
        # Validation above must finish before CP2K changes process-global state.
        _runtime = self
        if self._external_mpi:
            self._lib.cp2k_init_without_mpi_comm(self._comm.py2f())
        else:
            self._lib.cp2k_init()
        atexit.register(self._at_exit)

    def _check(self):
        if os.getpid() != self._pid:
            raise RuntimeError("Do not use CP2K after fork; use a spawned process")
        if threading.current_thread() is not threading.main_thread():
            raise RuntimeError("CP2K must be used from Python's main thread")
        if self._closed:
            raise RuntimeError("The CP2K runtime is closed")
        if self._mpi is not None and self._mpi.Is_finalized():
            raise RuntimeError("MPI was finalized before CP2K was closed")
        if self._comm is not None and self._comm == self._mpi.COMM_NULL:
            raise RuntimeError("The MPI communicator was freed before CP2K was closed")

    @property
    def version(self):
        self._check()
        buffer = ct.create_string_buffer(256)
        self._lib.cp2k_get_version(buffer, len(buffer))
        return buffer.value.decode("utf-8")

    @contextmanager
    def _input_file(self, inp, input_file):
        if (inp is None) == (input_file is None):
            raise ValueError("Specify exactly one of inp= or input_file=")
        if input_file is not None:
            yield _path_bytes(input_file)
            return
        text = input_to_string(inp) if isinstance(inp, Mapping) else inp
        if not isinstance(text, str) or "\0" in text or not text.strip():
            raise ValueError("inp must be nonempty CP2K input text or a mapping")
        # Each rank owns its scratch file. The CP2K reader broadcasts root's
        # input. Keep relative @INCLUDE/data paths relative to the caller's cwd.
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", suffix=".inp", prefix="cp2k-python-", dir="."
        ) as handle:
            handle.write(text)
            handle.flush()
            yield _path_bytes(handle.name)

    def create_force_env(self, inp=None, *, input_file=None, output_file="__STD_OUT__"):
        """Create an environment from CP2K input text, a mapping, or a file.

        Files named in the input (including PROJECT output) are resolved by
        native CP2K in the current working directory, not the input directory.
        Close the previous environment before creating another one.
        """
        self._check()
        if self._environments:
            raise RuntimeError(
                "Close the current force environment before creating another"
            )
        output = _path_bytes(output_file, output=True)
        with self._input_file(inp, input_file) as source:
            handle = ct.c_int()
            if self._comm is None:
                self._lib.cp2k_create_force_env(ct.byref(handle), source, output)
            else:
                self._lib.cp2k_create_force_env_comm(
                    ct.byref(handle), source, output, self._comm.py2f()
                )
        env = ForceEnvironment(self, handle.value)
        self._environments.append(env)
        return env

    def run_input(self, inp=None, *, input_file=None, output_file="__STD_OUT__"):
        """Execute a complete native input, including MD or geometry optimization."""
        self._check()
        if self._environments:
            raise RuntimeError("Close all force environments before run_input()")
        output = _path_bytes(output_file, output=True)
        with self._input_file(inp, input_file) as source:
            if self._comm is None:
                self._lib.cp2k_run_input(source, output)
            else:
                self._lib.cp2k_run_input_comm(source, output, self._comm.py2f())

    def close(self):
        """Destroy environments and finalize CP2K, never caller-owned MPI."""
        if self._closed:
            return
        self._check()
        for env in self._environments[::-1]:
            env.close()
        if self._external_mpi:
            self._lib.cp2k_finalize_without_mpi()
        else:
            self._lib.cp2k_finalize()
        self._closed = True
        atexit.unregister(self._at_exit)

    def _at_exit(self):
        # Last resort only. Collective cleanup belongs in explicit with/close.
        if os.getpid() == self._pid and not self._closed:
            if self._mpi is None or not self._mpi.Is_finalized():
                self.close()

    def __enter__(self):
        self._check()
        return self

    def __exit__(self, *exc):
        self.close()


class ForceEnvironment:
    """A reusable native force environment. Create through CP2K.create_force_env.

    Arrays have shape (number of particles, 3), cells have lattice vectors in
    rows, as in ASE. All values use atomic units. Particle counts may differ
    from atom counts, e.g. for shell models. Getters return independent copies.
    """

    def __init__(self, runtime, handle):
        self._runtime = runtime
        self._handle = handle
        self._closed = False
        self._energy_valid = False
        self._forces_valid = False
        self._stress_valid = False
        self._scf_converged = None

    def _check(self):
        self._runtime._check()
        if self._closed:
            raise RuntimeError("The CP2K force environment is closed")

    @property
    def communicator(self):
        """Caller-side communicator, or None for a single-process caller."""
        return self._runtime._comm

    def _count(self, name):
        self._check()
        value = ct.c_int()
        getattr(self._runtime._lib, "cp2k_get_" + name)(self._handle, ct.byref(value))
        if value.value < 0 or value.value > _C_INT_MAX // 3:
            raise RuntimeError("Invalid particle count returned by libcp2k")
        return value.value

    @property
    def natom(self):
        return self._count("natom")

    @property
    def nparticle(self):
        return self._count("nparticle")

    def _get_array(self, name, shape, sized=True):
        self._check()
        result = np.empty(shape, dtype=np.float64)
        args = [self._handle, result.ctypes.data_as(_DP)]
        if sized:
            args.append(result.size)
        getattr(self._runtime._lib, "cp2k_get_" + name)(*args)
        return result

    def _set_array(self, name, values, shape, sized=True):
        self._check()
        array = _array(values, shape, name)
        if name == "cell":
            determinant = np.linalg.det(array)
            if not np.isfinite(determinant) or determinant <= 0:
                raise ValueError("cell must be nonsingular and right-handed")
        args = [self._handle, array.ctypes.data_as(_DP)]
        if sized:
            args.append(array.size)
        # Invalidate before entering native code; never expose stale results.
        self._energy_valid = self._forces_valid = self._stress_valid = False
        self._scf_converged = None
        getattr(self._runtime._lib, "cp2k_set_" + name)(*args)

    @property
    def positions(self):
        return self._get_array("positions", (self.nparticle, 3))

    @positions.setter
    def positions(self, values):
        self._set_array("positions", values, (self.nparticle, 3))

    @property
    def cell(self):
        return self._get_array("cell", (3, 3), sized=False)

    @cell.setter
    def cell(self, values):
        self._set_array("cell", values, (3, 3), sized=False)

    def set_velocities(self, values):
        """Set particle velocities in bohr per atomic unit of time."""
        self._set_array("velocities", values, (self.nparticle, 3))

    @property
    def potential_energy(self):
        self._check()
        if not self._energy_valid:
            raise RuntimeError(
                "Call calculate() after creating/changing the environment"
            )
        value = ct.c_double()
        self._runtime._lib.cp2k_get_potential_energy(self._handle, ct.byref(value))
        return value.value

    @property
    def forces(self):
        self._check()
        if not self._forces_valid:
            raise RuntimeError("Call calculate(forces=True) before requesting forces")
        return self._get_array("forces", (self.nparticle, 3))

    @property
    def stress(self):
        """Potential pressure-positive Cartesian stress in hartree/bohr**3."""
        self._check()
        if not self._stress_valid:
            raise RuntimeError("Call calculate(stress=True) before requesting stress")
        result = np.empty((3, 3), dtype=np.float64, order="F")
        available = ct.c_int()
        self._runtime._lib.cp2k_get_stress_tensor(
            self._handle, result.ctypes.data_as(_DP), ct.byref(available)
        )
        if not available.value:
            raise RuntimeError("Enable FORCE_EVAL/STRESS_TENSOR in the CP2K input")
        return result

    @property
    def virial(self):
        """Potential pressure-positive Cartesian virial in hartree."""
        return self.stress * np.linalg.det(self.cell)

    @property
    def scf_converged(self):
        """True/False for the last SCF, or None if invalidated/unavailable."""
        self._check()
        return self._scf_converged

    def calculate(self, *, forces=True, stress=False, check_convergence=True):
        """Recalculate, raising on a reported SCF failure by default.

        check_convergence=False permits diagnostic unconverged results. Unknown
        status (older libraries or unsupported solvers) remains None, not True.
        This does not suppress native aborts: to inspect failed SCF results,
        explicitly enable SCF/IGNORE_CONVERGENCE_FAILURE in the CP2K input.
        Stress requires STRESS_TENSOR input and also computes native forces.
        """
        self._check()
        self._energy_valid = self._forces_valid = self._stress_valid = False
        self._scf_converged = None
        if stress and not hasattr(self._runtime._lib, "cp2k_get_stress_tensor"):
            raise RuntimeError("This libcp2k lacks cp2k_get_stress_tensor")
        method = "cp2k_calc_energy_force" if forces or stress else "cp2k_calc_energy"
        getattr(self._runtime._lib, method)(self._handle)
        query = getattr(self._runtime._lib, "cp2k_get_scf_convergence", None)
        if query is not None:
            status = ct.c_int(-1)
            query(self._handle, ct.byref(status))
            if status.value not in (-1, 0, 1):
                raise RuntimeError("Invalid SCF convergence status from libcp2k")
            self._scf_converged = None if status.value == -1 else bool(status.value)
        if check_convergence and self._scf_converged is False:
            raise SCFConvergenceError(
                "CP2K SCF did not converge. Inspect the output and SCF settings; "
                "use check_convergence=False only to retrieve diagnostic results."
            )
        self._energy_valid = True
        self._forces_valid = bool(forces)
        self._stress_valid = bool(stress)
        try:
            result = CalculationResult(
                self.potential_energy,
                self.forces if forces else None,
                self.stress if stress else None,
                self.virial if stress else None,
                scf_converged=self._scf_converged,
            )
            if not np.isfinite(result.energy) or any(
                value is not None and not np.isfinite(value).all()
                for value in (result.forces, result.stress, result.virial)
            ):
                raise RuntimeError("CP2K returned non-finite energy, forces, or stress")
        except Exception:
            self._energy_valid = self._forces_valid = self._stress_valid = False
            raise
        return result

    def close(self):
        if self._closed:
            return
        self._check()
        self._runtime._lib.cp2k_destroy_force_env(self._handle)
        self._closed = True
        self._runtime._environments.remove(self)

    def __enter__(self):
        self._check()
        return self

    def __exit__(self, *exc):
        self.close()
