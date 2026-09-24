"""Numerical i-PI socket protocol tests (CP2K_TEST_EXECUTABLE required)."""

# SPDX-License-Identifier: GPL-2.0-or-later

from contextlib import contextmanager
from copy import deepcopy
import os
from pathlib import Path
import socket
import struct
import subprocess
import sys
import time
import uuid

import numpy as np
import pytest

from cp2k import input_to_string

pytestmark = pytest.mark.integration


@pytest.fixture
def lj_input():
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


def receive(sock, length):
    result = b""
    while len(result) < length:
        data = sock.recv(length - len(result))
        if not data:
            raise EOFError("CP2K driver closed the connection")
        result += data
    return result


def header(sock, command):
    # Fragment every command deliberately: TCP does not preserve message boundaries.
    for byte in command.encode().ljust(12):
        sock.sendall(bytes([byte]))


def geometry(sock, cell, positions, nat=None):
    header(sock, "POSDATA")
    sock.sendall(np.asarray(cell.T, dtype="=f8").tobytes())
    sock.sendall(np.asarray(np.linalg.inv(cell.T), dtype="=f8").tobytes())
    sock.sendall(struct.pack("=i", len(positions) if nat is None else nat))
    if nat is None:
        sock.sendall(np.asarray(positions, dtype="=f8").tobytes())


@contextmanager
def driver(tmp_path, inp, unix=False):
    executable = os.environ.get("CP2K_TEST_EXECUTABLE")
    if not executable:
        pytest.skip("Set CP2K_TEST_EXECUTABLE to test the i-PI driver")
    path = None
    with socket.socket(
        socket.AF_UNIX if unix else socket.AF_INET, socket.SOCK_STREAM
    ) as server:
        server.settimeout(30)
        if unix:
            name = "cp2k-" + uuid.uuid4().hex[:12]
            path = Path("/tmp/ipi_" + name)
            server.bind(str(path))
            settings = {"HOST": name, "UNIX": True, "PREFIX": "ipi"}
        else:
            server.bind(("127.0.0.1", 0))
            settings = {"HOST": "127.0.0.1", "PORT": server.getsockname()[1]}
        server.listen(1)
        inp = deepcopy(inp)
        inp["GLOBAL"].update(RUN_TYPE="DRIVER", PROJECT="driver")
        inp["MOTION"] = {"DRIVER": settings}
        source = tmp_path / "driver.inp"
        source.write_text(input_to_string(inp))
        with (tmp_path / "driver.log").open("w") as log:
            process = subprocess.Popen(
                [executable, "-i", str(source)],
                cwd=tmp_path,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
            try:
                connection, _ = server.accept()
                with connection:
                    connection.settimeout(30)
                    yield connection, process
            finally:
                if process.poll() is None:
                    process.terminate()
                process.wait(timeout=20)
                if path is not None:
                    path.unlink(missing_ok=True)


@pytest.mark.parametrize("unix", [False, True])
def test_driver_energy_force_virial(
    real_runtime, lj_input, tmp_path, monkeypatch, unix
):
    monkeypatch.chdir(tmp_path)
    with real_runtime.create_force_env(lj_input, output_file="reference.out") as env:
        cell, positions = env.cell, env.positions
        references = []
        for scale in (1.0, 1.02):
            env.cell = cell * scale
            env.positions = positions * scale
            reference = env.calculate()
            # Only this pair is within the cutoff; use the position-force outer product.
            # This uses the existing force API, not the optional stress getter.
            virial = np.outer(env.positions[1] - env.positions[0], reference.forces[1])
            references.append((reference, virial))
    with driver(tmp_path, lj_input, unix) as (sock, process):
        for scale, (reference, reference_virial) in zip((1.0, 1.02), references):
            header(sock, "STATUS")
            assert receive(sock, 12).strip() == b"READY"
            geometry(sock, cell * scale, positions * scale)
            header(sock, "STATUS")
            assert receive(sock, 12).strip() == b"HAVEDATA"
            header(sock, "GETFORCE")
            assert receive(sock, 12).strip() == b"FORCEREADY"
            energy = struct.unpack("=d", receive(sock, 8))[0]
            nat = struct.unpack("=i", receive(sock, 4))[0]
            assert nat == 2
            forces = np.frombuffer(receive(sock, 24 * nat), dtype="=f8").reshape(nat, 3)
            virial = np.frombuffer(receive(sock, 72), dtype="=f8").reshape(3, 3)
            assert struct.unpack("=i", receive(sock, 4))[0] == 0
            np.testing.assert_allclose(energy, reference.energy, rtol=1e-10)
            np.testing.assert_allclose(forces, reference.forces, rtol=1e-10)
            np.testing.assert_allclose(virial, reference_virial, rtol=1e-10)
        header(sock, "EXIT")
        assert process.wait(timeout=20) == 0


@pytest.mark.parametrize(
    "bad", ["GETFORCE", "UNKNOWN", "nat", "truncated", "positions", "cell", "duplicate"]
)
def test_driver_rejects_bad_protocol(lj_input, tmp_path, bad):
    with driver(tmp_path, lj_input) as (sock, process):
        if bad == "nat":
            geometry(sock, np.eye(3) * 40, np.zeros((2, 3)), nat=2**30)
        elif bad == "truncated":
            sock.sendall(b"POS")
            sock.shutdown(socket.SHUT_WR)
        elif bad == "positions":
            geometry(sock, np.eye(3) * 40, np.full((2, 3), np.nan))
        elif bad == "cell":
            geometry(sock, np.diag([-40.0, 40.0, 40.0]), np.zeros((2, 3)))
        elif bad == "duplicate":
            geometry(
                sock, np.eye(3) * 40, np.array([[8.0, 8.0, 8.0], [14.0, 9.0, 9.0]])
            )
            header(sock, "STATUS")
            assert receive(sock, 12).strip() == b"HAVEDATA"
            header(sock, "POSDATA")
        else:
            header(sock, bad)
        assert process.wait(timeout=20) != 0
    output = (tmp_path / "driver.log").read_text()
    expected = {
        "GETFORCE": "without POSDATA",
        "UNKNOWN": "Unknown i-PI message",
        "nat": "Particle number mismatch",
        "truncated": "Unexpected EOF",
        "positions": "Non-finite positions",
        "cell": "nonsingular and right-handed",
        "duplicate": "before GETFORCE",
    }
    assert expected[bad] in output


@pytest.mark.skipif(
    not os.environ.get("CP2K_TEST_IPI"),
    reason="Set CP2K_TEST_IPI for an installed i-PI server",
)
def test_ipi_server(lj_input, tmp_path):
    pytest.importorskip("ipi")
    executable = os.environ.get("CP2K_TEST_EXECUTABLE")
    if not executable:
        pytest.skip("Set CP2K_TEST_EXECUTABLE")
    name = "cp2k-" + uuid.uuid4().hex[:12]
    sockpath = Path("/tmp/ipi_" + name)
    (tmp_path / "init.xyz").write_text(
        "2\npositions{angstrom}\nAr 4 4 4\nAr 7.2 4.8 4.5\n"
    )
    (tmp_path / "input.xml").write_text(f"""<simulation verbosity='low'>
  <output prefix='simulation'>
    <properties stride='1' flush='1'>[step, potential, conserved]</properties>
  </output>
  <total_steps>3</total_steps>
  <prng><seed>1729</seed></prng>
  <ffsocket mode='unix' name='driver'>
    <address>{name}</address><latency>0.001</latency><timeout>20</timeout>
  </ffsocket>
  <system>
    <initialize nbeads='1'>
      <file mode='xyz'>init.xyz</file>
      <cell mode='manual' units='angstrom'>[20,1,2,0,21,3,0,0,22]</cell>
      <momenta mode='thermal' units='kelvin'>10</momenta>
    </initialize>
    <forces><force forcefield='driver'/></forces>
    <ensemble><temperature units='kelvin'>10</temperature></ensemble>
    <motion mode='dynamics'><dynamics mode='nve'>
      <timestep units='femtosecond'>0.1</timestep>
    </dynamics></motion>
  </system>
</simulation>
""")
    inp = deepcopy(lj_input)
    inp["GLOBAL"].update(RUN_TYPE="DRIVER", PROJECT="driver")
    inp["MOTION"] = {"DRIVER": {"HOST": name, "UNIX": True, "PREFIX": "ipi"}}
    (tmp_path / "driver.inp").write_text(input_to_string(inp))
    client = None
    with (tmp_path / "ipi.log").open("w") as log, (tmp_path / "driver.log").open(
        "w"
    ) as clientlog:
        server = subprocess.Popen(
            [sys.executable, str(Path(sys.executable).with_name("i-pi")), "input.xml"],
            cwd=tmp_path,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        try:
            deadline = time.monotonic() + 20
            while not sockpath.exists():
                assert (
                    server.poll() is None
                ), "i-PI failed before accepting a connection; see ipi.log"
                assert time.monotonic() < deadline, "i-PI startup timed out"
                time.sleep(0.05)
            client = subprocess.Popen(
                [executable, "-i", "driver.inp"],
                cwd=tmp_path,
                stdout=clientlog,
                stderr=subprocess.STDOUT,
            )
            assert server.wait(timeout=40) == 0
            assert client.wait(timeout=20) == 0
        finally:
            for process in (client, server):
                if process is not None and process.poll() is None:
                    process.terminate()
                    process.wait(timeout=10)
            sockpath.unlink(missing_ok=True)
    data = np.atleast_2d(np.loadtxt(tmp_path / "simulation.out"))
    assert data[-1, 0] == 3 and np.isfinite(data).all()
    assert np.ptp(data[:, 2]) < 1e-8
