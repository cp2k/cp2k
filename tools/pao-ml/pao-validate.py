#!/usr/bin/env python3

# author: Ole Schuett

import torch
import argparse
from pathlib import Path

from pao.dataset import PaoDataset
from pao.training import validate_model


# ======================================================================================
def main() -> None:
    description = "Validates a given equivariant PAO-ML model against test data."
    parser = argparse.ArgumentParser(description=description)
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
    stats = validate_model(model, dataset)
    print("validation minimum loss: {:.8e}".format(stats.min_loss))
    print("validation median  loss: {:.8e}".format(stats.median_loss))
    print("validation maximum loss: {:.8e}".format(stats.max_loss))
    assert stats.median_loss < args.threshold


# ======================================================================================
main()

# EOF
