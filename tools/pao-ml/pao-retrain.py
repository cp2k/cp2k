#!/usr/bin/env python3

# author: Ole Schuett

import argparse
import torch

from e3nn import o3  # type: ignore
from pathlib import Path
from torch.utils.data import DataLoader

from pao.model import PaoModel
from pao.dataset import PaoDataset
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
    assert model.pao_model_version >= 1
    print(f"Loaded pre-trained model from file: {args.model}")

    # Load the training data.
    dataset = PaoDataset(
        kind_name=model.kind_name,
        num_neighbors=model.num_neighbors,
        files=args.training_data,
    )
    dataloader = DataLoader(dataset, batch_size=args.batch, shuffle=True)
    print(f"Found {len(dataset)} training samples of kind '{model.kind_name}'.")

    # Check compatability between model and training data.
    assert dataset.kind.atomic_number == model.atomic_number
    assert dataset.kind.pao_basis_size == model.pao_basis_size
    assert dataset.kind.prim_basis_name == model.prim_basis_name
    assert dataset.kind.prim_basis_size == model.prim_basis_size
    assert dataset.feature_kind_names == model.feature_kind_names

    # Train the model.
    train_model(model, dataloader, args.epochs)

    # Save the model.
    model.save(args.model)
    print(f"Saved model to file: {args.model}")


# ======================================================================================
main()

# EOF
