#!/usr/bin/env python3

# author: Ole Schuett

import argparse
import torch
from pathlib import Path

from pao.dataset import PaoDataset, collate_fn
from pao.training import train_model


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="Re-trains an exiting PAO-ML model.")
    parser.add_argument("--epochs", type=int, default=5000)
    parser.add_argument("--batch", type=int, default=64)
    parser.add_argument("--model", type=Path, required=True)

    # Training data files are passed as positional arguments.
    parser.add_argument("training_data", type=Path, nargs="+")
    args = parser.parse_args()

    # Load existing model and ignore most cmd arguments.
    model = torch.jit.load(args.model)  # type: ignore
    assert model.pao_model_version == 2
    print(f"Loaded pre-trained model from file: {args.model}")

    # Load the training data.
    dataset = PaoDataset(
        kind_name=model.kind_name,
        num_layers=model.num_layers,
        cutoff=model.cutoff,
        files=args.training_data,
    )

    # Check compatability between model and training data.
    assert dataset.kind.atomic_number == model.atomic_number
    assert dataset.kind.pao_basis_size == model.pao_basis_size
    assert dataset.kind.prim_basis_name == model.prim_basis_name
    assert dataset.kind.prim_basis_size == model.prim_basis_size
    assert dataset.all_kind_names == model.all_kind_names

    # Train the model.
    train_model(model, dataset, epochs=args.epochs, batch_size=args.batch)

    # Save the model.
    torch.jit.save(model, args.model)  # type: ignore
    print(f"Saved model to file: {args.model}")


# ======================================================================================
main()

# EOF
