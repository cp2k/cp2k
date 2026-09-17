"""Independent adaptive Z2Pack/CP2K run for the nontrivial DFT SOC benchmark."""

import argparse
import json
from pathlib import Path

import z2pack

from cp2k_z2pack import CP2KSystem
from validate_cp2k import ROOT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("restart", type=Path)
    parser.add_argument("workdir", type=Path)
    args = parser.parse_args()
    root = args.workdir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    text = Path(__file__).with_name("examples").joinpath("stanene-soc.inp").read_text()
    text = text.replace("SCF_GUESS ATOMIC", "SCF_GUESS RESTART")
    text = text.replace(
        "KPOINTS_SOURCE WILSON",
        "KPOINTS_SOURCE NNKP\n        NNKP_FILE loop.nnkp\n        WILSON_LOOP T",
    )
    text = text.replace("Z2 T", "Z2 F").replace("SEED_NAME stanene", "SEED_NAME loop")
    template = root / "stanene-reference.inp"
    template.write_text(text)
    system = CP2KSystem(
        input_file=template,
        input_files=[args.restart],
        lattice=[[4.674, 0, 0], [2.337, 4.047802737288466, 0], [0, 0, 20]],
        command=[str(args.binary.resolve())],
        workdir=root / "lines",
        num_bands=8,
        env={"CP2K_DATA_DIR": str(ROOT / "data")},
    )
    result = z2pack.surface.run(
        system=system,
        surface=lambda s, t: [t, s / 2, 0],
        num_lines=11,
        min_neighbour_dist=1e-4,
        iterator=[16, 32, 64, 128, 256, 512],
        pos_tol=1e-3,
        save_file=str(root / "surface.json"),
        serializer=json,
    )
    # Also save synchronously so checkpoint failures cannot hide in a worker thread.
    z2pack.io.save(result, str(root / "surface.json"), serializer=json)
    assert (root / "surface.json").is_file()
    assert not result.convergence_report["line"]["PosCheck"]["FAILED"]
    for check in ("MoveCheck", "GapCheck"):
        assert not result.convergence_report["surface"][check]["FAILED"]
    invariant = z2pack.invariant.z2(result)
    assert invariant == 1
    assert z2pack.invariant.z2(z2pack.io.load(str(root / "surface.json"))) == invariant
    report = {
        "z2pack_z2": invariant,
        "all_convergence_checks_passed": True,
        "convergence": str(result.convergence_report),
        "cp2k_line_runs": len(list((root / "lines").glob("loop-*"))),
    }
    (root / "validation.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
