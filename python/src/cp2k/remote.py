"""Client for an independently launched, persistent CP2K (MPI) server."""

# SPDX-License-Identifier: GPL-2.0-or-later

from __future__ import annotations

import os
import socket
import threading
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from numpy.typing import ArrayLike

from ._transport import receive, root_call, send
from .library import CalculationResult, FloatArray, _array

if TYPE_CHECKING:
    from mpi4py import MPI


class SocketEnvironment:
    """Connect to ``python -m cp2k.server``; never load CP2K in the caller.

    The server owns CP2K and may use a different MPI implementation/rank count.
    ``comm`` is the *client* communicator (e.g. LAMMPS's); only its rank 0 uses
    the socket. All client ranks must call construction, calculate and close in
    the same order with consistent data. For OpenMM omit comm. Geometry setters
    are local, with one combined request at calculate(). Units match libcp2k.

    Connections are authenticated but NOT encrypted: use loopback or a trusted
    SSH tunnel, not an untrusted network. Both peers need the same token (at
    least 32 characters). Server-side exceptions, including SCF convergence
    failures, become RuntimeError with the original error in the message and
    close the connection. Timeouts/connection failures are also terminal.
    close() requests server shutdown but does not manage the external MPI launcher.
    """

    def __init__(
        self,
        host: str,
        port: int,
        *,
        token: str,
        comm: MPI.Intracomm | None = None,
        timeout: float = 600,
    ) -> None:
        if not isinstance(token, str) or len(token) < 32:
            raise ValueError("Use a secret token of at least 32 characters")
        if not np.isfinite(timeout) or timeout <= 0:
            raise ValueError("timeout must be positive and finite")
        if threading.current_thread() is not threading.main_thread():
            raise RuntimeError("Use SocketEnvironment from Python's main thread")
        if comm is not None and (comm.Is_inter() or comm.Get_size() == 0):
            raise ValueError("comm must be a live intracommunicator")
        self.communicator = comm
        self._socket: socket.socket | None = None
        self._closed = False
        self._pid = os.getpid()
        self._result: CalculationResult | None = None

        def connect() -> dict[str, Any]:
            self._socket = socket.create_connection((host, port), timeout)
            send(self._socket, {"protocol": 1, "token": token})
            return receive(self._socket)

        try:
            info = root_call(comm, connect)
            if "error" in info:
                raise RuntimeError(info["error"])
            self.natom, self.nparticle = int(info["natom"]), int(info["nparticle"])
            if not 0 < self.natom <= self.nparticle:
                raise ValueError("Invalid server particle counts")
            self.server_ranks = int(info["ranks"])
            self.cell = info["cell"]
            self.positions = info["positions"]
        except Exception:
            self._disconnect()
            raise

    def _check(self) -> None:
        if os.getpid() != self._pid:
            raise RuntimeError("Do not use a CP2K socket after fork")
        if threading.current_thread() is not threading.main_thread():
            raise RuntimeError("Use SocketEnvironment from Python's main thread")
        if self._closed:
            raise RuntimeError("The CP2K socket environment is closed")

    def _request(self, request: dict[str, Any]) -> dict[str, Any]:
        self._check()

        def exchange() -> dict[str, Any]:
            assert self._socket is not None
            send(self._socket, request)
            return receive(self._socket)

        try:
            reply = root_call(self.communicator, exchange)
            if "error" in reply:
                raise RuntimeError(reply["error"])
            return reply
        except Exception:
            self._disconnect()
            raise

    @property
    def positions(self) -> FloatArray:
        self._check()
        return self._positions.copy()

    @positions.setter
    def positions(self, values: ArrayLike) -> None:
        self._check()
        self._positions = _array(values, (self.nparticle, 3), "positions")
        self._result = None

    @property
    def cell(self) -> FloatArray:
        self._check()
        return self._cell.copy()

    @cell.setter
    def cell(self, values: ArrayLike) -> None:
        self._check()
        cell = _array(values, (3, 3), "cell")
        determinant = np.linalg.det(cell)
        if not np.isfinite(determinant) or determinant <= 0:
            raise ValueError("cell must be nonsingular and right-handed")
        self._cell = cell
        self._result = None

    def calculate(
        self,
        *,
        forces: bool = True,
        stress: bool = False,
        check_convergence: bool = True,
    ) -> CalculationResult:
        self._result = None
        reply = self._request(
            {
                "command": "calculate",
                "positions": self.positions.tolist(),
                "cell": self.cell.tolist(),
                "forces": bool(forces),
                "stress": bool(stress),
                "check_convergence": bool(check_convergence),
            }
        )
        energy = float(reply["energy"])
        if not np.isfinite(energy):
            raise RuntimeError("CP2K returned non-finite energy")
        scf_converged = reply.get("scf_converged")
        if scf_converged is not None and not isinstance(scf_converged, bool):
            raise RuntimeError("Invalid SCF convergence status from CP2K server")
        if check_convergence and scf_converged is False:
            raise RuntimeError("CP2K server returned an unconverged SCF result")
        result = CalculationResult(
            energy,
            _array(reply["forces"], (self.nparticle, 3), "forces") if forces else None,
            _array(reply["stress"], (3, 3), "stress") if stress else None,
            _array(reply["virial"], (3, 3), "virial") if stress else None,
            scf_converged=scf_converged,
        )
        self._result = result
        return CalculationResult(
            result.energy,
            None if result.forces is None else result.forces.copy(),
            None if result.stress is None else result.stress.copy(),
            None if result.virial is None else result.virial.copy(),
            scf_converged=result.scf_converged,
        )

    @property
    def scf_converged(self) -> bool | None:
        self._check()
        return None if self._result is None else self._result.scf_converged

    def _get_result(self, name: str) -> float | FloatArray:
        self._check()
        value = None if self._result is None else getattr(self._result, name)
        if value is None:
            raise RuntimeError(f"Call calculate() requesting {name} first")
        return (
            cast(FloatArray, value.copy())
            if isinstance(value, np.ndarray)
            else cast(float, value)
        )

    @property
    def potential_energy(self) -> float:
        return cast(float, self._get_result("energy"))

    @property
    def forces(self) -> FloatArray:
        return cast(FloatArray, self._get_result("forces"))

    @property
    def stress(self) -> FloatArray:
        return cast(FloatArray, self._get_result("stress"))

    @property
    def virial(self) -> FloatArray:
        return cast(FloatArray, self._get_result("virial"))

    def _disconnect(self) -> None:
        if self._socket is not None:
            self._socket.close()
            self._socket = None
        self._closed = True
        self._result = None

    def close(self) -> None:
        if not self._closed:
            try:
                self._request({"command": "close"})
            finally:
                self._disconnect()

    def __enter__(self) -> SocketEnvironment:
        self._check()
        return self

    def __exit__(self, exc_type: type[BaseException] | None, *exc: object) -> None:
        if exc_type is None:
            self.close()
        else:
            # Disconnecting also makes the server leave its receive loop.
            self._disconnect()
