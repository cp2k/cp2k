"""Check full-basis SOC overlap identities in a bonded multi-atom periodic cell."""

import argparse
import json
from pathlib import Path
import re
import numpy as np
from cp2k_z2pack import nnkp_text, read_loop_mmn
from validate_cp2k import run


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("workdir", type=Path)
    args = parser.parse_args()
    root = args.workdir.resolve()
    base = Path(__file__).with_name("examples").joinpath("stanene-soc.inp").read_text()
    base = re.sub(r"^\s*EXCLUDE_BANDS.*\n", "", base, flags=re.M)
    base = base.replace(
        "KPOINTS_SOURCE WILSON", "KPOINTS_SOURCE NNKP\n        NNKP_FILE loop.nnkp"
    )
    base = base.replace("Z2 T", "Z2 F").replace("SEED_NAME stanene", "SEED_NAME loop")
    cell = [[4.674, 0, 0], [2.337, 4.047802737288466, 0], [0, 0, 20]]
    points = [[0.13, 0.17, 0], [0.39, 0.24, 0], [0.13, 0.17, 0]]
    report = {}
    for name in ("adjoint", "zero-step"):
        directory = root / name
        directory.mkdir(parents=True, exist_ok=True)
        text = nnkp_text(points, cell)
        if name == "zero-step":
            coords = (
                text.split("begin kpoints")[1]
                .split("end kpoints")[0]
                .strip()
                .splitlines()
            )
            text = text.replace(coords[2], coords[1])
        (directory / "loop.nnkp").write_text(text)
        run(args.binary.resolve(), directory, base)
        matrices = read_loop_mmn(directory / "loop.mmn", points, 52)
        if name == "adjoint":
            error = float(np.max(np.abs(matrices[0] - matrices[1].conj().T)))
        else:
            error = float(max(np.max(np.abs(m - np.eye(52))) for m in matrices))
        report[name] = error
        print(name, error, flush=True)
    assert max(report.values()) < 1e-8
    (root / "identities.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
