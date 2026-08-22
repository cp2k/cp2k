#!/usr/bin/env python3
"""Regenerate data/Symmetrix/toy_mace-1-8.json, the toy model this directory tests against.

The model is untrained: its weights come from a fixed seed, so it carries no
physics and the reference energies prove only that the backend is deterministic
and that the atom decomposition is exact. That is what they are for.

Note this does not go through `python -m symmetrix.cli.extract_mace`: that CLI
has no option for the spline resolution and writes indented JSON, whereas the
committed model uses 128 spline points and a compact encoding to stay small.

Needs `mace-torch`, `torch`, `e3nn` and a `symmetrix` checkout on PYTHONPATH.
Bit-identical output requires the same versions (developed against mace-torch
0.3.16); a regenerated model that differs will change the reference values in
TEST_FILES.toml, which then have to be updated together.

Usage:
    python make_toy_model.py [-o ../../../data/Symmetrix/toy_mace-1-8.json]
"""

import argparse
import json

import numpy as np
import torch
from e3nn import o3
from mace import modules
from symmetrix.extract_mace_data import extract_mace_data

NUM_SPLINE_POINTS = 128


def build_model() -> torch.nn.Module:
    """Two-element (H, O) MACE with a fixed seed."""
    torch.manual_seed(42)
    torch.set_default_dtype(torch.float64)
    return modules.ScaleShiftMACE(
        r_max=4.0,
        num_bessel=4,
        num_polynomial_cutoff=5,
        max_ell=3,
        interaction_cls=modules.interaction_classes[
            "RealAgnosticResidualInteractionBlock"
        ],
        # A non-residual first block is what symmetrix supports.
        interaction_cls_first=modules.interaction_classes[
            "RealAgnosticInteractionBlock"
        ],
        num_interactions=2,
        num_elements=2,
        hidden_irreps=o3.Irreps("2x0e"),
        MLP_irreps=o3.Irreps("16x0e"),
        atomic_energies=np.array([-0.5, -75.0]),
        avg_num_neighbors=8.0,
        atomic_numbers=[1, 8],
        correlation=3,
        gate=torch.nn.functional.silu,
        atomic_inter_scale=1.0,
        atomic_inter_shift=0.0,
    ).double()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", default="toy_mace-1-8.json")
    args = parser.parse_args()

    model = build_model()
    data = extract_mace_data(model, ["H", "O"], num_spline_points=NUM_SPLINE_POINTS)
    with open(args.output, "w") as f:
        json.dump(data, f, separators=(",", ":"))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
