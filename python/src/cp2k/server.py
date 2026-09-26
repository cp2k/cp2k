"""Serve one CP2K environment to an authenticated client, optionally using MPI."""

# SPDX-License-Identifier: GPL-2.0-or-later

import argparse
from pathlib import Path
import secrets
import socket
import sys

import numpy as np

from ._transport import receive, root_call, send
from .library import CP2K, _array


def serve(environment, *, token, host="127.0.0.1", port=8765, timeout=600):
    """Collective server loop; caller retains environment/runtime ownership.

    One client per server. Loopback by default; use an SSH tunnel for remote
    machines. The secret authenticates but does not encrypt the connection.
    Native aborts close the connection; Python/transport failures are terminal
    and are reported collectively. No client-supplied input or code is executed.
    """
    if not isinstance(token, str) or len(token) < 32:
        raise ValueError("Use a secret token of at least 32 characters")
    if not np.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout must be positive and finite")
    environment._check()
    comm = environment.communicator
    root = comm is None or comm.rank == 0
    connection = listener = None

    def accept():
        nonlocal connection, listener
        listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        listener.bind((host, port))
        listener.listen(1)
        listener.settimeout(timeout)
        print(
            f"CP2K server listening on {host}:{listener.getsockname()[1]}",
            file=sys.stderr,
            flush=True,
        )
        connection, _ = listener.accept()
        connection.settimeout(timeout)
        hello = receive(connection)
        supplied = hello.get("token")
        if (
            hello.get("protocol") != 1
            or not isinstance(supplied, str)
            or not secrets.compare_digest(supplied.encode(), token.encode())
        ):
            send(connection, {"error": "Invalid CP2K protocol or authentication"})
            raise ValueError("Invalid CP2K protocol or authentication")

    try:
        root_call(comm, accept)
        info = {
            "natom": environment.natom,
            "nparticle": environment.nparticle,
            "positions": environment.positions.tolist(),
            "cell": environment.cell.tolist(),
            "ranks": 1 if comm is None else comm.size,
        }
        root_call(comm, lambda: send(connection, info))
        while True:
            request = root_call(comm, lambda: receive(connection))
            if request.get("command") == "close":
                root_call(comm, lambda: send(connection, {"closed": True}))
                return
            error, result = None, None
            try:
                if request.get("command") != "calculate":
                    raise ValueError("Unknown CP2K server command")
                positions = _array(
                    request["positions"], (environment.nparticle, 3), "positions"
                )
                cell = _array(request["cell"], (3, 3), "cell")
                determinant = np.linalg.det(cell)
                if not np.isfinite(determinant) or determinant <= 0:
                    raise ValueError("cell must be nonsingular and right-handed")
                if not all(
                    isinstance(request[key], bool) for key in ("forces", "stress")
                ):
                    raise ValueError("forces and stress must be boolean")
                check_convergence = request.get("check_convergence", True)
                if not isinstance(check_convergence, bool):
                    raise ValueError("check_convergence must be boolean")
                environment.cell = cell
                environment.positions = positions
                result = environment.calculate(
                    forces=request["forces"],
                    stress=request["stress"],
                    check_convergence=check_convergence,
                )
            except Exception as caught:
                error = f"{type(caught).__name__}: {caught}"
            errors = [error] if comm is None else comm.allgather(error)
            failures = [item for item in errors if item is not None]
            if failures:
                root_call(comm, lambda: send(connection, {"error": failures[0]}))
                return
            reply = {"energy": result.energy, "scf_converged": result.scf_converged}
            for name in ("forces", "stress", "virial"):
                value = getattr(result, name)
                reply[name] = None if value is None else value.tolist()
            root_call(comm, lambda: send(connection, reply))
    finally:
        if root:
            if connection is not None:
                connection.close()
            if listener is not None:
                listener.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="cp2k-server.out")
    parser.add_argument("--library", default=None)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("--token-file", required=True)
    parser.add_argument(
        "--mpi",
        action="store_true",
        help="Use mpi4py COMM_WORLD; required with mpiexec",
    )
    args = parser.parse_args()
    token = Path(args.token_file).read_text().strip()
    comm = None
    if args.mpi:
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
    with CP2K(library=args.library, comm=comm) as runtime:
        with runtime.create_force_env(
            input_file=args.input, output_file=args.output
        ) as environment:
            serve(
                environment,
                token=token,
                host=args.host,
                port=args.port,
                timeout=args.timeout,
            )


if __name__ == "__main__":
    main()
