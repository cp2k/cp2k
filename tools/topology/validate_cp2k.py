"""Run CP2K/Z2Pack/native comparisons; retain all inputs, outputs and a JSON report."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import z2pack
from cp2k_z2pack import CP2KSystem, nnkp_text

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/regtesting"))
from compare_wannier90_mmn import read_mmn, compare_mmn


def run(binary, directory, text):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "input.inp").write_text(text)
    with (directory / "run.log").open("w") as log:
        subprocess.run(
            [str(binary), "-i", "input.inp"],
            cwd=directory,
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            env={
                **os.environ,
                "CP2K_DATA_DIR": str(ROOT / "data"),
                "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "1"),
                "OPENBLAS_NUM_THREADS": "1",
            },
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("workdir", type=Path)
    parser.add_argument(
        "--soc", action="store_true", help="Also run the SOC surface validation"
    )
    args = parser.parse_args()
    binary = args.binary.resolve()
    workdir = args.workdir.resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    template = Path(__file__).with_name("examples").joinpath("helium.inp")
    base = template.read_text()
    system = CP2KSystem(
        input_file=template,
        lattice=4 * np.eye(3),
        command=[str(binary)],
        workdir=workdir / "z2pack",
        num_bands=1,
        env={"CP2K_DATA_DIR": str(ROOT / "data")},
    )
    result = z2pack.line.run(
        system=system, line=lambda t: [t, 0, 0], iterator=[8, 16, 32], pos_tol=1e-5
    )
    native = np.loadtxt(system.last_run / "loop.wilson", ndmin=2)[0, 2:]
    difference = np.max(
        np.abs(np.exp(2j * np.pi * native) - np.exp(2j * np.pi * np.array(result.wcc)))
    )
    assert difference < 1e-10
    # Finite GPW grid accuracy, separate from native/Z2Pack algebra agreement.
    assert abs(np.exp(2j * np.pi * native[0]) + 1) < 1e-3

    # Export the same complete mesh twice: automatic neighbours versus imported connections.
    regular = workdir / "regular"
    regular_text = base.replace("KPOINTS_SOURCE NNKP", "KPOINTS_SOURCE SCF").replace(
        "WILSON_LOOP T", "WILSON_LOOP F"
    )
    run(binary, regular, regular_text)
    mesh = read_mmn(str(regular / "loop.mmn"))
    win = (regular / "loop.win").read_text()
    coords = win.split("begin kpoints")[1].split("end kpoints")[0].strip().splitlines()
    # The lattice is copied from the independently generated loop helper.
    prefix = nnkp_text([[0, 0, 0], [0.5, 0, 0], [1, 0, 0]], 4 * np.eye(3)).split(
        "begin kpoints"
    )[0]
    lines = [
        prefix,
        "begin kpoints",
        str(mesh[1]),
        *coords,
        "end kpoints",
        "begin nnkpts",
        str(mesh[2]),
    ]
    lines += [" ".join(map(str, header)) for header, _ in mesh[3]]
    lines += ["end nnkpts", ""]
    explicit = workdir / "explicit"
    explicit.mkdir(exist_ok=True)
    (explicit / "loop.nnkp").write_text("\n".join(lines))
    run(binary, explicit, base.replace("WILSON_LOOP T", "WILSON_LOOP F"))
    raw, singular, *_ = compare_mmn(mesh, read_mmn(str(explicit / "loop.mmn")))
    # Diagnostic only: the legacy path assumes a Hermitian AO operator, whereas
    # directed cross-k overlaps do not have that symmetry. Physical adjoint and
    # zero-step checks are in validate_overlap_geometry.py.

    # Internal mesh refinement, with no Z2Pack call.
    native_surface = workdir / "native-surface"
    run(
        binary,
        native_surface,
        base.replace(
            "KPOINTS_SOURCE NNKP",
            "KPOINTS_SOURCE WILSON\n"
            "        WILSON_MESH 4 3\n        WILSON_MAX_REFINEMENT 2",
        ),
    )
    assert (
        "Wilson surface sampling converged" in (native_surface / "run.log").read_text()
    )
    report = {
        "native_z2pack_circle_error": float(difference),
        "atomic_wcc": native.tolist(),
        "regular_explicit_raw_error": raw,
        "regular_explicit_singular_value_error": singular,
        "native_surface_converged": True,
        "z2pack_convergence": str(result.convergence_report),
        "soc_validation": "not requested",
    }
    if args.soc:
        report.update(validate_soc(binary, workdir, template))
    (workdir / "validation.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


def validate_soc(binary, workdir, template):

    # A closed-shell SOC insulator, including the real DFT -> Z2Pack surface path.
    soc_base = template.with_name("neon-soc.inp").read_text()
    soc_native = workdir / "soc-native"
    run(binary, soc_native, soc_base)
    assert "Converged Z2 invariant: 0" in (soc_native / "run.log").read_text()
    soc_template = workdir / "soc-reference.inp"
    soc_template.write_text(
        soc_base.replace(
            "KPOINTS_SOURCE WILSON",
            "KPOINTS_SOURCE NNKP\n"
            "        NNKP_FILE loop.nnkp\n        WILSON_LOOP T",
        )
        .replace("Z2 T", "Z2 F")
        .replace("SEED_NAME neon", "SEED_NAME loop")
    )
    soc_system = CP2KSystem(
        input_file=soc_template,
        lattice=5 * np.eye(3),
        command=[str(binary)],
        workdir=workdir / "soc-z2pack",
        num_bands=8,
        env={"CP2K_DATA_DIR": str(ROOT / "data")},
    )
    soc_result = z2pack.surface.run(
        system=soc_system,
        surface=lambda s, t: [t, s / 2, 0],
        num_lines=5,
        min_neighbour_dist=0.01,
        iterator=[8, 16, 32],
        pos_tol=1e-3,
    )
    assert z2pack.invariant.z2(soc_result) == 0
    assert not soc_result.convergence_report["line"]["PosCheck"]["FAILED"]
    for check in ("MoveCheck", "GapCheck"):
        assert not soc_result.convergence_report["surface"][check]["FAILED"]
    return {
        "soc_validation": "passed",
        "soc_native_z2": 0,
        "soc_z2pack_z2": z2pack.invariant.z2(soc_result),
        "soc_z2pack_convergence": str(soc_result.convergence_report),
    }


if __name__ == "__main__":
    main()
