"""Bounded JSON framing; no pickle or native-library objects cross the socket."""

# SPDX-License-Identifier: GPL-2.0-or-later

import json
import struct

_MAX_MESSAGE = 64 * 1024 * 1024


def _read(sock, size):
    data = bytearray()
    while len(data) < size:
        chunk = sock.recv(min(size - len(data), 65536))
        if not chunk:
            raise ConnectionError("CP2K socket closed before a complete response")
        data.extend(chunk)
    return data


def receive(sock):
    size = struct.unpack("!I", _read(sock, 4))[0]
    if not 0 < size <= _MAX_MESSAGE:
        raise ValueError("Invalid CP2K message length")
    value = json.loads(_read(sock, size))
    if not isinstance(value, dict):
        raise ValueError("CP2K messages must be JSON objects")
    return value


def send(sock, value):
    data = json.dumps(value, allow_nan=False, separators=(",", ":")).encode()
    if not 0 < len(data) <= _MAX_MESSAGE:
        raise ValueError("CP2K message exceeds the 64 MiB limit")
    sock.sendall(struct.pack("!I", len(data)) + data)


def root_call(comm, function):
    """Relay root-only I/O failures to all caller ranks before raising."""
    reply = None
    if comm is None or comm.rank == 0:
        try:
            reply = {"value": function()}
        except Exception as error:
            reply = {"error": f"{type(error).__name__}: {error}"}
    if comm is not None:
        reply = comm.bcast(reply, root=0)
    if "error" in reply:
        raise RuntimeError(reply["error"])
    return reply["value"]
