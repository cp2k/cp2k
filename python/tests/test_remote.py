"""Transport failures, independent MPI servers and native adapter regressions."""

# SPDX-License-Identifier: GPL-2.0-or-later

from contextlib import contextmanager
from copy import deepcopy
import os
import json
from pathlib import Path
import re
import secrets
import socket
import struct
import subprocess
import sys
import threading
import time

import numpy as np
import pytest

from cp2k import SocketEnvironment, input_to_string
from cp2k._transport import receive, send


def test_fragmented_transport():
    left, right = socket.socketpair()
    try:
        data = b'{"energy":1.25}'
        chunks = struct.pack("!I", len(data)) + data

        def write():
            for byte in chunks:
                left.sendall(bytes([byte]))

        thread = threading.Thread(target=write)
        thread.start()
        assert receive(right) == {"energy": 1.25}
        thread.join(timeout=5)
        send(left, {"positions": [[1, 2, 3]]})
        assert receive(right)["positions"] == [[1, 2, 3]]
        left.sendall(struct.pack("!I", 100000000))
        with pytest.raises(ValueError, match="length"):
            receive(right)
        left.close()
        with pytest.raises(ConnectionError, match="closed"):
            receive(right)
    finally:
        left.close()
        right.close()


def test_remote_validation():
    with pytest.raises(ValueError, match="token"):
        SocketEnvironment("127.0.0.1", 1, token="short")
    with pytest.raises(ValueError, match="timeout"):
        SocketEnvironment("127.0.0.1", 1, token="x" * 32, timeout=-1)


def test_openmm_native_qmmm(h2_input, lj_input, server_factory):
    openmm = pytest.importorskip("openmm", minversion="8.6.1")
    from cp2k.openmm import create_force
    from cp2k.qmmm import SubtractiveQMMM
    from cp2k._units import BOHR_TO_NM, HARTREE_TO_KJMOL

    # Actual DFT H2 replaces an intentionally simple H-H MM reference. This
    # checks real backend coupling, not the physical quality of the toy MM model.
    low = deepcopy(lj_input)
    low["GLOBAL"]["PROJECT"] = "h2-mm"
    low["FORCE_EVAL"]["SUBSYS"] = {
        "CELL": {"ABC": [8, 8, 8]},
        "COORD": {"_lines": ["H 3.6 4 4", "H 4.4 4 4"]},
        "KIND": {"_": "H", "ELEMENT": "H"},
        "TOPOLOGY": {"CONN_FILE_FORMAT": "OFF"},
    }
    low["FORCE_EVAL"]["MM"]["FORCEFIELD"] = {
        "SPLINE": {"EPS_SPLINE": 1e-12},
        "CHARGE": {"ATOM": "H", "CHARGE": 0},
        "NONBONDED": {
            "LENNARD-JONES": {
                "ATOMS": ["H", "H"],
                "EPSILON": "[hartree] 0.001",
                "SIGMA": 0.7,
                "RCUT": 3,
            }
        },
    }
    with server_factory(h2_input) as (high_config, _), server_factory(low) as (
        low_config,
        _,
    ):
        with SocketEnvironment(**high_config) as high, SocketEnvironment(
            **low_config
        ) as low:
            positions = np.zeros((3, 3))
            positions[[2, 0]] = high.positions
            positions[1] = [1, 2, 3]  # spectator with an independent MM potential
            model = SubtractiveQMMM(
                high, low, qm_atoms=[2, 0], positions=positions, cell=high.cell
            )
            qm = high.calculate()
            reference = low.calculate()
            system = openmm.System()
            for _ in range(3):
                system.addParticle(1)
            system.setDefaultPeriodicBoxVectors(*high.cell * BOHR_TO_NM)
            mm = openmm.CustomBondForce(
                "4*eps*((sigma/r)^12-(sigma/r)^6-(sigma/rc)^12+(sigma/rc)^6)"
            )
            mm.addGlobalParameter("eps", 0.001 * HARTREE_TO_KJMOL)
            mm.addGlobalParameter("sigma", 0.07)  # nm
            mm.addGlobalParameter("rc", 0.3)
            mm.addBond(2, 0, [])
            system.addForce(mm)
            spectator = openmm.CustomExternalForce("0.5*k*(x*x+y*y+z*z)")
            spectator.addGlobalParameter("k", HARTREE_TO_KJMOL / BOHR_TO_NM**2)
            spectator.addParticle(1, [])
            system.addForce(spectator)
            integrator = openmm.VerletIntegrator(0.000001)
            context = openmm.Context(
                system, integrator, openmm.Platform.getPlatformByName("Reference")
            )
            context.setPositions(positions * BOHR_TO_NM)
            state = context.getState(getEnergy=True, getForces=True)
            spectator_energy = np.sum(positions[1] ** 2) / 2
            assert state.getPotentialEnergy()._value == pytest.approx(
                (reference.energy + spectator_energy) * HARTREE_TO_KJMOL, abs=1e-6
            )
            np.testing.assert_allclose(
                state.getForces(asNumpy=True)._value[[2, 0]],
                reference.forces * HARTREE_TO_KJMOL / BOHR_TO_NM,
                rtol=1e-6,
                atol=1e-5,
            )
            system.addForce(create_force(model, periodic=True))
            context.reinitialize(preserveState=True)
            state = context.getState(getEnergy=True, getForces=True)
            assert state.getPotentialEnergy()._value == pytest.approx(
                (qm.energy + spectator_energy) * HARTREE_TO_KJMOL, abs=1e-5
            )
            expected = np.zeros((3, 3))
            expected[[2, 0]] = qm.forces
            expected[1] = -positions[1]
            np.testing.assert_allclose(
                state.getForces(asNumpy=True)._value,
                expected * HARTREE_TO_KJMOL / BOHR_TO_NM,
                rtol=1e-5,
                atol=0.02,
            )
            del context, integrator


@pytest.fixture
def server_factory(tmp_path):
    @contextmanager
    def launch(inp, *, ranks=1):
        library = os.environ.get("CP2K_TEST_LIBRARY")
        if not library:
            pytest.skip("Set CP2K_TEST_LIBRARY for native socket tests")
        if ranks > 1 and not os.environ.get("CP2K_TEST_REMOTE_MPI"):
            pytest.skip("Set CP2K_TEST_REMOTE_MPI for multi-rank CP2K servers")
        directory = tmp_path / f"server-{secrets.token_hex(4)}"
        directory.mkdir()
        (directory / "input.inp").write_text(input_to_string(inp))
        token = secrets.token_hex(32)
        (directory / "token").write_text(token)
        (directory / "token").chmod(0o600)
        command = [
            sys.executable,
            "-m",
            "cp2k.server",
            "--input",
            "input.inp",
            "--token-file",
            "token",
            "--port",
            "0",
            "--timeout",
            "30",
            "--library",
            library,
        ]
        if ranks > 1:
            command = (
                [os.environ.get("MPIEXEC", "mpiexec"), "-n", str(ranks)]
                + command
                + ["--mpi"]
            )
        process = None
        with (directory / "server.log").open("w") as log:
            try:
                process = subprocess.Popen(
                    command, cwd=directory, stdout=log, stderr=log
                )
                deadline = time.monotonic() + 30
                while True:
                    output = (directory / "server.log").read_text()
                    match = re.search(
                        r"CP2K server listening on 127.0.0.1:(\d+)", output
                    )
                    if match:
                        break
                    if process.poll() is not None or time.monotonic() > deadline:
                        pytest.fail(f"CP2K server failed to start:\n{output}")
                    time.sleep(0.02)
                yield dict(
                    host="127.0.0.1", port=int(match[1]), token=token, timeout=30
                ), process
            finally:
                if process is not None:
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.terminate()
                        try:
                            process.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            process.kill()
                            process.wait(timeout=5)

    return launch


@pytest.mark.parametrize("ranks", [1, 2])
def test_remote_native(
    real_runtime, lj_input, server_factory, ranks, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    with real_runtime.create_force_env(lj_input, output_file="reference.out") as direct:
        expected = direct.calculate(stress=True)
    with server_factory(lj_input, ranks=ranks) as (config, process):
        with SocketEnvironment(**config) as remote:
            assert remote.server_ranks == ranks
            actual = remote.calculate(stress=True)
            assert actual.energy == pytest.approx(expected.energy, abs=1e-12)
            np.testing.assert_allclose(actual.forces, expected.forces, atol=1e-12)
            np.testing.assert_allclose(actual.virial, expected.virial, atol=1e-12)
            actual.forces[:] = 0
            np.testing.assert_allclose(remote.forces, expected.forces, atol=1e-12)
            remote.positions = remote.positions + 0.01
            with pytest.raises(RuntimeError, match="calculate"):
                _ = remote.potential_energy
            assert np.isfinite(remote.calculate(forces=False).energy)
            with pytest.raises(RuntimeError, match="calculate"):
                _ = remote.forces
            with pytest.raises(ValueError, match="right-handed"):
                remote.cell = -np.eye(3)
        assert process.wait(timeout=10) == 0
        remote.close()
        with pytest.raises(RuntimeError, match="closed"):
            remote.calculate()


def test_remote_authentication(lj_input, server_factory):
    with server_factory(lj_input) as (config, process):
        config["token"] = "wrong" * 10
        with pytest.raises(RuntimeError, match="authentication"):
            SocketEnvironment(**config)
        assert process.wait(timeout=10) != 0


def test_remote_abort(lj_input, server_factory):
    with server_factory(lj_input) as (config, process):
        remote = SocketEnvironment(**config)
        process.terminate()
        process.wait(timeout=10)
        with pytest.raises(RuntimeError):
            remote.calculate()
        remote.close()
        with pytest.raises(RuntimeError, match="closed"):
            remote.calculate()


def test_remote_server_error(lj_input, server_factory):
    lj_input["FORCE_EVAL"]["STRESS_TENSOR"] = "NONE"
    with server_factory(lj_input) as (config, process):
        remote = SocketEnvironment(**config)
        with pytest.raises(RuntimeError, match="STRESS_TENSOR"):
            remote.calculate(stress=True)
        assert process.wait(timeout=10) == 0
        with pytest.raises(RuntimeError, match="closed"):
            remote.calculate()


def test_lammps_remote_mpi(lj_input, server_factory, tmp_path):
    if not os.environ.get("CP2K_TEST_LAMMPS") or not os.environ.get(
        "CP2K_TEST_REMOTE_MPI"
    ):
        pytest.skip("Set CP2K_TEST_LAMMPS and CP2K_TEST_REMOTE_MPI")
    # Two LAMMPS client ranks, independent single-rank CP2K server. In particular,
    # no CP2K communicator is passed into the client and only rank 0 uses TCP.
    with server_factory(lj_input) as (config, process):
        path = tmp_path / "connection.json"
        path.write_text(json.dumps(config))
        path.chmod(0o600)
        script = Path(__file__).with_name("lammps_mpi_smoke.py")
        run = subprocess.run(
            [
                os.environ.get("MPIEXEC", "mpiexec"),
                "-n",
                "2",
                sys.executable,
                str(script),
                "--socket-config",
                str(path),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            timeout=45,
        )
        assert run.returncode == 0, run.stdout + run.stderr
        assert process.wait(timeout=10) == 0


@pytest.mark.parametrize("units", ["metal", "real"])
def test_lammps_remote(lj_input, server_factory, units):
    if not os.environ.get("CP2K_TEST_LAMMPS"):
        pytest.skip("Set CP2K_TEST_LAMMPS")
    from lammps import lammps
    from cp2k.lammps import ExternalForce
    from cp2k._units import BOHR_TO_ANGSTROM, HARTREE_TO_EV, HARTREE_TO_KCALMOL

    factor = HARTREE_TO_EV if units == "metal" else HARTREE_TO_KCALMOL
    with server_factory(lj_input) as (config, process):
        lmp = lammps(cmdargs=["-log", "none", "-screen", "none"])
        try:
            with SocketEnvironment(**config, comm=lmp.get_mpi_comm()) as remote:
                expected = remote.calculate(stress=True)
                lmp.commands_string(f"""
units {units}
atom_style atomic
boundary p p p
region box prism 0 20 0 21 0 22 1 2 3
create_box 1 box
create_atoms 1 single 7.2 4.8 4.5
create_atoms 1 single 4 4 4
mass 1 39.948
pair_style zero 8
pair_coeff * *
compute cp_pressure all pressure NULL virial
thermo_style custom step pe c_cp_pressure[1] c_cp_pressure[2] c_cp_pressure[3] c_cp_pressure[4] c_cp_pressure[5] c_cp_pressure[6]
thermo_modify norm no
""")
                with ExternalForce(lmp, remote, atom_ids=[2, 1]) as coupling:
                    coupling.run(0)
                    assert lmp.get_thermo("pe") == pytest.approx(
                        expected.energy * factor, abs=1e-10
                    )
                    tags = lmp.numpy.extract_atom("id")[:2]
                    indices = [1 if tag == 1 else 0 for tag in tags]
                    np.testing.assert_allclose(
                        lmp.numpy.extract_atom("f")[:2],
                        expected.forces[indices] * factor / BOHR_TO_ANGSTROM,
                        atol=1e-10,
                    )
                    pressure = (
                        expected.virial.flat[[0, 4, 8, 1, 2, 5]]
                        * factor
                        / np.linalg.det(remote.cell * BOHR_TO_ANGSTROM)
                        * lmp.extract_global("nktv2p")
                    )
                    np.testing.assert_allclose(
                        lmp.numpy.extract_compute("cp_pressure", 0, 1),
                        pressure,
                        atol=1e-10,
                    )
                    coupling.command(
                        "velocity all create 10 731 mom yes rot no dist gaussian"
                    )
                    coupling.command(
                        "fix thermostat all npt temp 10 10 100 iso 0 0 1000"
                    )
                    coupling.command(
                        "timestep " + ("0.0001" if units == "metal" else "0.1")
                    )
                    coupling.run(3)
                    assert np.isfinite(lmp.get_thermo("pe"))
        finally:
            lmp.close()
        assert process.wait(timeout=10) == 0


@pytest.mark.parametrize("ranks", [1, 2])
@pytest.mark.parametrize("platform", ["Reference", "CPU"])
def test_openmm_remote(lj_input, server_factory, ranks, platform):
    openmm = pytest.importorskip("openmm", minversion="8.6.1")
    from openmm import unit
    from cp2k.openmm import create_force
    from cp2k._units import BOHR_TO_NM, HARTREE_TO_KJMOL

    with server_factory(lj_input, ranks=ranks) as (config, process):
        with SocketEnvironment(**config) as remote:
            expected = remote.calculate()
            system = openmm.System()
            for _ in range(remote.natom):
                system.addParticle(39.948)
            system.setDefaultPeriodicBoxVectors(*remote.cell * BOHR_TO_NM)
            system.addForce(create_force(remote, periodic=True))
            integrator = openmm.VerletIntegrator(0.0001)
            context = openmm.Context(
                system, integrator, openmm.Platform.getPlatformByName(platform)
            )
            context.setPositions(remote.positions * BOHR_TO_NM)
            state = context.getState(getEnergy=True, getForces=True)
            assert state.getPotentialEnergy().value_in_unit(
                unit.kilojoules_per_mole
            ) == pytest.approx(expected.energy * HARTREE_TO_KJMOL, abs=2e-6)
            np.testing.assert_allclose(
                state.getForces(asNumpy=True)._value,
                expected.forces * HARTREE_TO_KJMOL / BOHR_TO_NM,
                rtol=2e-5,
                atol=2e-5,
            )
            integrator.step(3)
            cell = remote.cell * 1.01
            context.setPeriodicBoxVectors(*cell * BOHR_TO_NM)
            assert np.isfinite(
                context.getState(getEnergy=True).getPotentialEnergy()._value
            )
            np.testing.assert_allclose(remote.cell, cell)
            del context, integrator
        assert process.wait(timeout=10) == 0
