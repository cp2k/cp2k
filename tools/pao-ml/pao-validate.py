#!/usr/bin/env python3

# author: Ole Schuett

import torch
import argparse
import numpy as np
from pathlib import Path

from pao.dataset import PaoDataset
from pao.training import loss_function


# ======================================================================================
def main() -> None:
    description = "Validates a given equivariant PAO-ML model against test data."
    parser = argparse.ArgumentParser(description=description)
    # parser.add_argument("--batch", type=int, default=64)
    parser.add_argument("--model", required=True)
    parser.add_argument("--threshold", type=float, default=1.0)

    # Test data files are passed as positional arguments.
    parser.add_argument("test_data", type=Path, nargs="+")
    args = parser.parse_args()

    # Load model.
    model = torch.jit.load(args.model)  # type: ignore
    assert model.pao_model_version == 2
    print(f"Loaded model from file: {args.model}")

    # Load the test data.
    dataset = PaoDataset(
        kind_name=model.kind_name,
        num_layers=model.num_layers,
        cutoff=model.cutoff,
        files=args.test_data,
    )

    # Check compatability between model and test data.
    assert dataset.kind.atomic_number == model.atomic_number
    assert dataset.kind.pao_basis_size == model.pao_basis_size
    assert dataset.kind.prim_basis_name == model.prim_basis_name
    assert dataset.kind.prim_basis_size == model.prim_basis_size
    assert dataset.all_kind_names == model.all_kind_names

    # Compute losses.
    losses = []
    for example in dataset:
        outputs = model(example)
        loss = loss_function(outputs["xblock"], example["xblock"])
        losses.append(loss.item())

    print("minimum loss: {:.8e}".format(np.amin(losses).item()))
    print("median  loss: {:.8e}".format(np.median(losses).item()))
    print("maximum loss: {:.8e}".format(np.amax(losses).item()))

    assert np.median(losses).item() < args.threshold


# ======================================================================================
main()

# EOF
