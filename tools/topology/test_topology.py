"""Regression tests of native Fortran algebra and the CP2K/Z2Pack adapter."""

import os
from pathlib import Path
import shlex
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest
import z2pack
from z2pack._utils import _gapfind
from cp2k_z2pack import CP2KSystem, nnkp_text, read_loop_mmn

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="session")
def native(tmp_path_factory):
    directory = tmp_path_factory.mktemp("wilson-build")
    compiler = shutil.which(os.environ.get("FC", "gfortran"))
    if compiler is None:
        pytest.fail("Fortran compiler unavailable; set FC to a GNU Fortran compiler")
    libraries = shlex.split(os.environ.get("WILSON_LAPACK_FLAGS", "-llapack -lblas"))
    exe = directory / "wilson_driver"
    subprocess.run(
        [
            compiler,
            "-cpp",
            "-ffree-form",
            "-fcheck=all",
            "-O0",
            "-g",
            str(ROOT / "src/base/kinds.F"),
            str(ROOT / "src/topology_wilson.F"),
            str(Path(__file__).with_name("wilson_driver.f90")),
            *libraries,
            "-o",
            str(exe),
        ],
        cwd=directory,
        check=True,
        capture_output=True,
        text=True,
    )
    return exe


def native_result(native, loops, check=True):
    n = len(loops[0][0])
    lines = [f"{n} {len(loops[0])} {len(loops)}"]
    for loop in loops:
        for m in loop:
            lines.extend(f"{z.real:.17g} {z.imag:.17g}" for z in m.T.ravel())
    result = subprocess.run(
        [str(native)],
        input="\n".join(lines) + "\n",
        capture_output=True,
        text=True,
        check=check,
        env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
    )
    wcc = [
        np.array(line.split()[1:], float)
        for line in result.stdout.splitlines()
        if line.startswith("WCC")
    ]
    z2 = [
        tuple(map(int, line.split()[1:]))
        for line in result.stdout.splitlines()
        if line.strip().startswith("Z2")
    ]
    return result, wcc, z2


def bhz(kx, ky, mass):
    def spin_block(x, y):
        d = mass + np.cos(x) + np.cos(y)
        off = np.sin(x) - 1j * np.sin(y)
        return np.array([[d, off], [off.conjugate(), -d]])

    h = np.zeros((4, 4), complex)
    h[:2, :2] = spin_block(kx, ky)
    h[2:, 2:] = spin_block(-kx, -ky).conj()
    return h


def model_loops(mass, npoints=41, nlines=31, gauge=False):
    rng = np.random.default_rng(191)
    loops = []
    for ky in np.linspace(0, np.pi, nlines):
        states = []
        for kx in np.linspace(0, 2 * np.pi, npoints, endpoint=False):
            _, c = np.linalg.eigh(bhz(kx, ky, mass))
            c = c[:, :2]
            if gauge:
                u, _ = np.linalg.qr(
                    rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
                )
                c = c @ u
            states.append(c)
        loops.append([a.conj().T @ b for a, b in zip(states, states[1:] + states[:1])])
    return loops


@pytest.mark.parametrize("mass,expected", [(-1, 1), (-3, 0), (1, 1), (3, 0)])
def test_native_z2_against_known_model_and_z2pack(native, mass, expected):
    loops = model_loops(mass)
    _, wcc, parity = native_result(native, loops)
    assert parity == [(expected, 0)]
    ref = []
    for loop in loops:
        polar = []
        for m in loop:
            u, _, vh = np.linalg.svd(m)
            polar.append(u @ vh)
        ref.append(z2pack.line.OverlapLineData(polar).wcc)
    # Compare eigenvalues on the unit circle, avoiding the 0/1 branch cut.
    for a, b in zip(wcc, ref):
        assert np.allclose(
            sorted(np.exp(2j * np.pi * a), key=np.angle),
            sorted(np.exp(2j * np.pi * np.array(b)), key=np.angle),
            atol=1e-10,
        )
    result = SimpleNamespace(wcc=ref, gap_pos=[_gapfind(x)[0] for x in ref])
    assert z2pack.invariant.z2(result) == expected


def test_random_nonabelian_gauge_and_reversal(native):
    loops = model_loops(-1, gauge=True)
    _, wcc, parity = native_result(native, loops)
    _, plain, plain_parity = native_result(native, model_loops(-1))
    assert parity == plain_parity == [(1, 0)]
    for a, b in zip(wcc, plain):
        assert abs(sum(np.exp(2j * np.pi * a)) - sum(np.exp(2j * np.pi * b))) < 1e-10
    reverse = [[m.conj().T for m in loop[::-1]] for loop in loops]
    _, reversed_wcc, reverse_parity = native_result(native, reverse)
    assert reverse_parity == parity
    for a, b in zip(wcc, reversed_wcc):
        assert (
            abs(sum(np.exp(2j * np.pi * a)).conjugate() - sum(np.exp(2j * np.pi * b)))
            < 1e-10
        )


def test_singular_link_rejected(native):
    result, _, _ = native_result(native, [[np.zeros((2, 2), complex)]], check=False)
    assert result.returncode != 0 and "-2" in result.stdout


def test_boundary_kramers_pairs_required(native):
    loop = [np.diag(np.exp(2j * np.pi * np.array([0.1, 0.3])))]
    _, _, parity = native_result(native, [loop, loop])
    assert parity == [(-1, -2)]


def test_unresolved_surface_is_not_accepted(native):
    loops = [
        [np.diag(np.exp(2j * np.pi * np.array(w)))]
        for w in ([0.1, 0.1], [0.6, 0.6], [0.1, 0.1])
    ]
    result, _, _ = native_result(native, loops)
    assert "RESOLVED F" in result.stdout


def test_nnkp_nonorthogonal_cell_and_winding():
    cell = [[2, 0, 0], [1, 3, 0], [0, 0, 12]]
    text = nnkp_text([[0.2, 0, 0], [0.6, 0, 0], [1.2, 0, 0]], cell)
    assert "2 1 1 0 0" in text
    assert "begin recip_lattice" in text
    with pytest.raises(ValueError, match="endpoint"):
        nnkp_text([[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0]], cell)
    with pytest.raises(ValueError, match="Nonfinite"):
        nnkp_text([[0, 0, 0], [float("nan"), 0, 0], [0, 0, 0]], cell)


def test_mmn_strict_closure_and_column_order(tmp_path):
    path = tmp_path / "test.mmn"
    points = [[0, 0, 0], [0.5, 0, 0], [1, 0, 0]]
    path.write_text(
        "test\n2 2 1\n1 2 0 0 0\n1 0\n2 1\n3 2\n4 0\n" "2 1 1 0 0\n1 0\n0 0\n0 0\n1 0\n"
    )
    m = read_loop_mmn(path, points, 2)
    assert m[0][0, 1] == 3 + 2j and m[0][1, 0] == 2 + 1j
    with pytest.raises(ValueError, match="bands"):
        read_loop_mmn(path, points, 1)
    path.write_text(path.read_text().replace("2 1 1 0 0", "2 1 0 0 0"))
    with pytest.raises(ValueError, match="connection"):
        read_loop_mmn(path, points, 2)


def test_adapter_retains_requests_and_uses_argument_list(tmp_path, monkeypatch):
    template = tmp_path / "template.inp"
    template.write_text("test input\n")
    restart = tmp_path / "restart.kp"
    restart.write_text("restart\n")
    calls = []

    def execute(argv, **kwargs):
        directory = kwargs["cwd"]
        assert argv == ["cp2k.psmp", "-i", "input.inp"]
        assert kwargs.get("shell", False) is False
        assert (directory / "input.inp").read_text() == "test input\n"
        assert (directory / "restart.kp").read_text() == "restart\n"
        assert "2 1 1 0 0" in (directory / "loop.nnkp").read_text()
        (directory / "loop.mmn").write_text(
            "test\n1 2 1\n1 2 0 0 0\n0.8 0\n2 1 1 0 0\n0 0.9\n"
        )
        calls.append(directory)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(subprocess, "run", execute)
    system = CP2KSystem(
        input_file=template,
        lattice=np.eye(3),
        command=["cp2k.psmp"],
        workdir=tmp_path / "runs",
        input_files=[restart],
        num_bands=1,
    )
    for _ in range(2):
        matrices = system.get_mmn([[0, 0, 0], [0.5, 0, 0], [1, 0, 0]])
        assert np.allclose(matrices, [[[1]], [[1j]]])
        assert (system.last_run / "diagnostics.json").is_file()
        assert (system.last_run / "request.json").is_file()
    assert calls[0] != calls[1]


def test_adapter_reports_failed_calculations(tmp_path, monkeypatch):
    template = tmp_path / "template.inp"
    template.write_text("test input\n")

    def execute(argv, **kwargs):
        kwargs["stdout"].write("SCF failed\n")
        return SimpleNamespace(returncode=7)

    monkeypatch.setattr(subprocess, "run", execute)
    system = CP2KSystem(
        input_file=template,
        lattice=np.eye(3),
        command=["cp2k.psmp"],
        workdir=tmp_path / "runs",
    )
    with pytest.raises(RuntimeError, match=r"CP2K failed \(7\)"):
        system.get_mmn([[0, 0, 0], [0.5, 0, 0], [1, 0, 0]])
    assert (system.last_run / "run.log").read_text() == "SCF failed\n"
