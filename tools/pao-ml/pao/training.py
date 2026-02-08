# author: Ole Schuett


from typing import Any, Union, Iterable
from dataclasses import dataclass
import numpy as np
import torch
from torch.utils.data import DataLoader, random_split

from .model import PaoModel
from .dataset import PaoDataset, collate_fn


# ======================================================================================
def loss_function(prediction: torch.Tensor, label: torch.Tensor) -> torch.Tensor:
    # This assumes the columns of prediction and label are orthonormal.
    p1 = prediction.transpose(-2, -1) @ prediction
    p2 = label.transpose(-2, -1) @ label
    return (p1 - p2).pow(2).mean()


# ======================================================================================
def train_model(
    model: Union[PaoModel, torch.jit.ScriptModule],
    dataset: PaoDataset,
    epochs: int,
    batch_size: int,
) -> None:

    # Split dataset into training and validation set.
    generator = torch.Generator().manual_seed(42)
    [training_set, validation_set] = random_split(dataset, [0.8, 0.2], generator)
    validation_percent = 100 * len(validation_set) / len(dataset)
    print(f"Reserving {validation_percent:.0f}% of dataset for validation.")
    assert validation_percent > 0, "Need at least 5 samples."

    # Train the model.
    training_loader = DataLoader(
        training_set, batch_size=batch_size, shuffle=True, collate_fn=collate_fn
    )
    optim = torch.optim.Adam(model.parameters())
    for epoch in range(epochs + 1):
        optim.zero_grad()
        for example in training_loader:
            outputs = model(example)
            loss = loss_function(outputs["xblock"], example["xblock"])
            loss.backward()  # type: ignore
            optim.step()
        if epoch % 100 == 0:
            validation_loss = validate_model(model, validation_set).mean_loss
            print(
                f"epoch: {epoch:5d} | training loss: {loss:.8e} | validation loss: {validation_loss:.8e}"
            )

    validation_loss = validate_model(model, validation_set).mean_loss
    print(
        f"Training complete, final training loss: {loss:.8e}, final validation loss: {validation_loss:.8e}"
    )


# ======================================================================================
@dataclass
class ValidationStats:
    min_loss: float
    median_loss: float
    mean_loss: float
    max_loss: float


# ======================================================================================
def validate_model(
    model: Union[PaoModel, torch.jit.ScriptModule], dataset: Any
) -> ValidationStats:
    losses = []
    with torch.no_grad():
        for example in dataset:
            outputs = model(example)
            loss = loss_function(outputs["xblock"], example["xblock"])
            losses.append(loss.item())

    return ValidationStats(
        min_loss=np.amin(losses).item(),
        median_loss=np.median(losses).item(),
        mean_loss=np.mean(losses).item(),
        max_loss=np.amax(losses).item(),
    )


# EOF
