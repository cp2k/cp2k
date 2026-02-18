#!/usr/bin/env python3

# author: Ole Schuett

import torch
import argparse
import e3nn.util.jit  # type: ignore
from pathlib import Path

from pao.model import PaoModel
from pao.dataset import PaoDataset
from pao.training import train_model


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="Trains an equivariant PAO-ML model.")
    parser.add_argument("--kind", required=True)
    parser.add_argument("--epochs", type=int, default=5000)
    parser.add_argument("--batch", type=int, default=64)
    parser.add_argument("--model", type=Path, default=None)

    # Hyper-parameters - TODO tune default values
    parser.add_argument(
        "--neighbors",
        type=int,
        default=5,
        help="Average number of neighbors used to normalize the edge sums.",
    )
    parser.add_argument(
        "--layers",
        type=int,
        default=1,
        help="Equals the number of message passing hops between neighboring atoms.",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=3.0,
        help="Cutoff distance in Ångström beyond which atoms are no longer neighbors, aka. r_max",
    )

    # Training data files are passed as positional arguments.
    parser.add_argument("training_data", type=Path, nargs="+")
    args = parser.parse_args()

    # Load the training data.
    dataset = PaoDataset(
        kind_name=args.kind,
        num_layers=args.layers,
        cutoff=args.cutoff,
        files=args.training_data,
    )

    # Construct the model.
    model_py = PaoModel(
        kind_name=dataset.kind.name,
        atomic_number=dataset.kind.atomic_number,
        prim_basis_name=dataset.kind.prim_basis_name,
        prim_basis_size=dataset.kind.prim_basis_size,
        pao_basis_size=dataset.kind.pao_basis_size,
        all_kind_names=dataset.all_kind_names,
        num_layers=args.layers,
        cutoff=args.cutoff,
        avg_num_neighbors=args.neighbors,
        # TODO expose other hyper-parameters, e.g. num_features, radial_mlp, etc.
    )

    # Compile the model to TorchScript.
    # TODO: Try PyTorch's Ahead-of-Time Inductor compiler.
    model = e3nn.util.jit.script(model_py)

    num_model_params = sum(p.numel() for p in model.parameters())
    print(f"PAO-ML model will have {num_model_params} parameters.")

    # Train the model.
    train_model(model, dataset, epochs=args.epochs, batch_size=args.batch)

    # Save the model.
    default_fn = f"{dataset.kind.prim_basis_name}-PAO{dataset.kind.pao_basis_size}-{args.kind}.pt"
    fn = args.model or default_fn
    torch.jit.save(model, fn)
    print(f"Saved model to file: {fn}")


# ======================================================================================
main()

# EOF
