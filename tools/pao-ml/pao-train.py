#!/usr/bin/env python3

# author: Ole Schuett

import torch
import argparse
from typing import cast
from pathlib import Path
from torch.utils.data import DataLoader

from pao.model import PaoModel
from pao.dataset import PaoDataset, PaoRecord
from pao.training import train_model


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="Trains an equivariant PAO-ML model.")
    parser.add_argument("--kind", required=True)
    parser.add_argument("--epochs", type=int, default=5000)
    parser.add_argument("--batch", type=int, default=64)
    parser.add_argument("--model", type=Path, default=None)

    # Hyper-parameters - TODO tune default values
    parser.add_argument("--neighbors", type=int, default=5)
    parser.add_argument("--distances", type=int, default=10)
    parser.add_argument("--layers", type=int, default=16)
    parser.add_argument("--cutoff", type=float, default=6.0)

    # Training data files are passed as positional arguments.
    parser.add_argument("training_data", type=Path, nargs="+")
    args = parser.parse_args()

    # Load the training data.
    dataset = PaoDataset(
        kind_name=args.kind, num_neighbors=args.neighbors, files=args.training_data
    )
    dataloader = DataLoader(dataset, batch_size=args.batch, shuffle=True)
    print(f"Found {len(dataset)} training samples of kind '{args.kind}'.")

    # Construct the model.
    model_py = PaoModel(
        kind_name=dataset.kind.name,
        atomic_number=dataset.kind.atomic_number,
        prim_basis_name=dataset.kind.prim_basis_name,
        prim_basis_size=dataset.kind.prim_basis_size,
        pao_basis_size=dataset.kind.pao_basis_size,
        feature_kind_names=dataset.feature_kind_names,
        num_neighbors=args.neighbors,
        num_distances=args.distances,
        num_layers=args.layers,
        cutoff=args.cutoff,
    )

    # Compile the model to TorchScript.
    model = cast(PaoModel, torch.jit.script(model_py))

    num_model_params = sum(p.numel() for p in model.parameters())
    print(f"PAO-ML model will have {num_model_params} parameters.")

    # Train the model.
    train_model(model, dataloader, args.epochs)

    # Save the model.
    default_fn = f"{dataset.kind.prim_basis_name}-PAO{dataset.kind.pao_basis_size}-{args.kind}.pt"
    fn = args.model or default_fn
    model.save(fn)
    print(f"Saved model to file: {fn}")


# ======================================================================================
main()

# EOF
