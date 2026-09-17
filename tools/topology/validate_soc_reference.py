"""Compare the full-basis topology spinors to CP2K's existing SOC band structure."""

import argparse
import json
from pathlib import Path
import re

import numpy as np
from cp2k_z2pack import nnkp_text
from validate_cp2k import run


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("workdir", type=Path)
    parser.add_argument(
        "--stanene", action="store_true", help="Also exercise interatomic SOC blocks"
    )
    args = parser.parse_args()
    root = args.workdir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    name = "stanene-soc.inp" if args.stanene else "neon-soc.inp"
    base = Path(__file__).with_name("examples").joinpath(name).read_text()
    base = base.replace("ADDED_MOS 1", "ADDED_MOS -1")
    base = base.replace("MONKHORST-PACK 2 2 2", "MONKHORST-PACK 4 4 4")
    base = base.replace("MONKHORST-PACK 8 8 1", "MONKHORST-PACK 4 4 4")
    base, removed = re.subn(
        r"^    &PRINT\n.*?^    &END(?: PRINT)?[ \t]*\n", "", base, flags=re.S | re.M
    )
    assert removed == 1, "Expected exactly one DFT print section in the template"
    # Minimal DOS mesh plus a short path: the reference prints the path endpoints.
    properties = """  &PROPERTIES
    &BANDSTRUCTURE
      &BANDSTRUCTURE_PATH
        NPOINTS 2
        SPECIAL_POINT Gamma 0 0 0
        SPECIAL_POINT X 0.5 0 0
      &END
      &DOS
        KPOINTS 1 1 1
      &END
      &SOC
      &END
    &END
  &END
"""
    reference = root / "reference"
    run(
        args.binary.resolve(),
        reference,
        base.replace("  &SUBSYS", properties + "  &SUBSYS"),
    )
    output = (reference / "bandstructure_SCF_and_G0W0_plus_SOC").read_text()
    blocks = output.split("kpoint:")[1:]
    points, values = [], []
    for block in blocks:
        points.append(
            [float(x) for x in block.split("coordinate:")[1].splitlines()[0].split()]
        )
        values.append(
            [
                float(line.split()[-1])
                for line in block.splitlines()
                if "(occ)" in line or "(vir)" in line
            ]
        )
    assert points and len(set(map(len, values))) == 1
    nspinor = len(values[0])
    native = root / "native"
    native.mkdir(exist_ok=True)
    # A directed loop containing all reference points. No topology is inferred from this path.
    loop_points = points + [np.asarray(points[0]) + [1, 0, 0]]
    cell = (
        [[4.674, 0, 0], [2.337, 4.047802737288466, 0], [0, 0, 20]]
        if args.stanene
        else 5 * np.eye(3)
    )
    (native / "loop.nnkp").write_text(nnkp_text(loop_points, cell))
    native_print = """    &PRINT
      &WANNIER90
        KPOINTS_SOURCE NNKP
        NNKP_FILE loop.nnkp
        SEED_NAME loop
        SOC T
      &END
    &END
"""
    run(
        args.binary.resolve(),
        native,
        base.replace("    &XC\n", native_print + "    &XC\n"),
    )
    actual = np.loadtxt(native / "loop.eig")[:, 2].reshape(len(points), nspinor)
    error = float(np.max(np.abs(actual - np.asarray(values))))
    report = {
        "reference_points": points,
        "spinor_bands": nspinor,
        "maximum_eigenvalue_error_eV": error,
        "reference_output_resolution_eV": 0.001,
        "reference_gamma_occupied_eV": values[0][:8],
        "native_gamma_occupied_eV": actual[0, :8].tolist(),
    }
    assert (
        error < 0.00051
    ), "SOC spectra disagree beyond the reference text output precision"
    (root / "comparison.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
