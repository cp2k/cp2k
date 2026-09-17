"""DFT checks of cross-k overlap adjoints, normalization, and cell/origin phases."""

import argparse
import json
from pathlib import Path

import numpy as np
from cp2k_z2pack import nnkp_text, read_loop_mmn
from validate_cp2k import run


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("workdir", type=Path)
    args = parser.parse_args()
    root = args.workdir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    template = Path(__file__).with_name("examples").joinpath("helium.inp").read_text()
    cell = np.array([[5.0, 0.0, 0.0], [1.1, 5.0, 0.0], [0.3, 0.7, 5.0]])
    base = template.replace(
        "ABC 4.0 4.0 4.0", "A 5 0 0\n      B 1.1 5 0\n      C .3 .7 5"
    )
    base = base.replace("He 0.5 0.5 0.5", "He 0.23 0.31 0.47")
    points = [[0.13, 0.17, 0.09], [0.39, 0.24, -0.08], [0.13, 0.17, 0.09]]
    pair = root / "adjoint"
    pair.mkdir(exist_ok=True)
    (pair / "loop.nnkp").write_text(nnkp_text(points, cell))
    run(args.binary.resolve(), pair, base)
    matrices = read_loop_mmn(pair / "loop.mmn", points, 1)
    adjoint_error = float(np.max(np.abs(matrices[0] - matrices[1].conj().T)))
    assert adjoint_error < 1e-10

    # At b=0 the cross-k operator must reduce to the AO metric S(k).
    zero = root / "zero-step"
    zero.mkdir(exist_ok=True)
    text = nnkp_text(points, cell)
    coords = text.split("begin kpoints")[1].split("end kpoints")[0].strip().splitlines()
    text = text.replace(coords[2], coords[1])
    (zero / "loop.nnkp").write_text(text)
    run(args.binary.resolve(), zero, base.replace("WILSON_LOOP T", "WILSON_LOOP F"))
    zero_m = read_loop_mmn(zero / "loop.mmn", points, 1)
    metric_error = float(
        max(abs(np.linalg.svd(m, compute_uv=False)[0] - 1) for m in zero_m)
    )
    assert metric_error < 1e-9

    phases = []
    for direction in (1, -1):
        winding = root / f"winding-{direction}"
        winding.mkdir(exist_ok=True)
        k = np.column_stack([np.linspace(0, direction, 17), np.zeros((17, 2))])
        (winding / "loop.nnkp").write_text(nnkp_text(k, cell))
        run(args.binary.resolve(), winding, base)
        phases.append(float(np.loadtxt(winding / "loop.wilson", ndmin=2)[0, 2]))
    phase_error = float(abs(np.exp(2j * np.pi * sum(phases)) - 1))
    assert phase_error < 1e-9
    # Berry phase of an atomic orbital at tau is -2*pi*tau for these ordered links.
    centre_error = float(abs(np.exp(2j * np.pi * (phases[0] + 0.23)) - 1))
    assert (
        centre_error < 0.002
    )  # deliberately modest real-space grid in this smoke test
    report = {
        "adjoint_error": adjoint_error,
        "zero_step_metric_error": metric_error,
        "opposite_winding_error": phase_error,
        "wcc": phases,
        "atomic_phase_error": centre_error,
    }
    (root / "geometry.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
