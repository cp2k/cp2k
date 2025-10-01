# author: Ole Schuett

import torch
from torch.utils.data import DataLoader
from typing import Union, Iterable
from nequip.data import AtomicDataDict  # type: ignore

from .model import PaoModel


# ======================================================================================
def loss_function(prediction: torch.Tensor, label: torch.Tensor) -> torch.Tensor:
    # This assumes the columns of prediction and label are orthonormal.
    p1 = prediction.transpose(-2, -1) @ prediction
    p2 = label.transpose(-2, -1) @ label
    return (p1 - p2).pow(2).mean()


# ======================================================================================
def train_model(
    model: Union[PaoModel, torch.jit.ScriptModule],
    data_source: Iterable[AtomicDataDict],
    epochs: int,
) -> None:
    # Train the model.
    optim = torch.optim.Adam(model.parameters())
    for epoch in range(epochs + 1):
        optim.zero_grad()
        for example in data_source:
            outputs = model(example)
            loss = loss_function(outputs["xblock"], example["xblock"])
            loss.backward()  # type: ignore
            optim.step()
        if epoch % 100 == 0:
            print(f"epoch: {epoch:5d} | loss: {loss:.8e}")

    print(f"Training complete, final loss: {loss:.8e}")


# EOF
