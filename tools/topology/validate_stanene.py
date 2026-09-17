"""Nontrivial DFT SOC benchmark: native surface versus Z2Pack overlap analysis."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import z2pack
from z2pack._utils import _gapfind

from validate_cp2k import read_mmn, run


def compare(directory):
    log = (directory / "run.log").read_text()
    assert "Converged Z2 invariant: 1" in log
    native = np.loadtxt(directory / "stanene.wilson", ndmin=2)
    _, npoints, nneighbours, entries = read_mmn(str(directory / "stanene.mmn"))
    assert nneighbours == 1
    nlines = len(native)
    assert npoints % nlines == 0
    per_line = npoints // nlines
    reference, raw_reference = [], []
    maximum_error = 0.0
    for iline in range(nlines):
        loops = [m for _, m in entries[iline * per_line : (iline + 1) * per_line]]
        polar = []
        for overlap in loops:
            u, _, vh = np.linalg.svd(overlap)
            polar.append(u @ vh)
        centres = np.asarray(z2pack.line.OverlapLineData(polar).wcc)
        reference.append(centres)
        raw_reference.append(z2pack.line.OverlapLineData(loops).wcc)
        # Cyclic matching avoids a spurious discrepancy at the 0/1 branch cut.
        a = np.sort(native[iline, 2:])
        error = min(
            np.max(
                np.abs(
                    np.exp(2j * np.pi * a)
                    - np.exp(2j * np.pi * np.roll(centres, shift))
                )
            )
            for shift in range(len(a))
        )
        maximum_error = max(maximum_error, float(error))
    assert maximum_error < 1e-9
    for centres in (reference, raw_reference):
        surface = SimpleNamespace(
            wcc=centres, gap_pos=[_gapfind(c)[0] for c in centres]
        )
        assert z2pack.invariant.z2(surface) == 1
    report = {
        "native_z2": 1,
        "z2pack_polar_z2": 1,
        "z2pack_raw_z2": 1,
        "native_surface_converged": True,
        "lines": nlines,
        "points_per_line": per_line,
        "maximum_wilson_eigenvalue_error": maximum_error,
        "minimum_link_singular_value": float(native[:, 1].min()),
        "scope": "Second-variational DFT SOC benchmark, not a basis/cutoff-converged material prediction",
    }
    (directory / "comparison.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("workdir", type=Path)
    parser.add_argument(
        "--reuse", action="store_true", help="Analyze an existing successful run"
    )
    args = parser.parse_args()
    directory = args.workdir.resolve()
    if not args.reuse:
        template = Path(__file__).with_name("examples") / "stanene-soc.inp"
        run(args.binary.resolve(), directory, template.read_text())
    compare(directory)


if __name__ == "__main__":
    main()
