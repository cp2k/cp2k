# author: Ole Schuett

import torch
import numpy as np
import numpy.typing as npt
import scipy.spatial  # type: ignore
from pathlib import Path
from typing import Any, List, Tuple, Iterator
from torch.utils.data import Dataset

from .io import parse_pao_file

PaoRecord = Tuple[torch.Tensor, torch.Tensor, torch.Tensor]


# ======================================================================================
def _as_tensor(x: List[Any] | npt.NDArray[np.float64]) -> torch.Tensor:
    return torch.tensor(np.array(x, dtype=np.float32))


# ======================================================================================
class PaoDataset(Dataset[PaoRecord]):
    def __init__(self, kind_name: str, num_neighbors: int, files: List[Path]):
        self.neighbors_relpos: List[torch.Tensor] = []
        self.neighbors_features: List[torch.Tensor] = []
        self.labels: List[torch.Tensor] = []

        # Load kinds from the first training data file.
        kinds = parse_pao_file(files[0]).kinds
        self.kind = kinds[kind_name]
        self.feature_kind_names = sorted(kinds.keys())

        # Load all training data files.
        for fn in files:
            f = parse_pao_file(fn)

            # Build  k-d tree of atom positions.
            assert 0 < num_neighbors < f.coords.shape[0]
            assert np.all(f.cell == np.diag(np.diagonal(f.cell)))
            boxsize = np.diagonal(f.cell)
            kdtree = scipy.spatial.KDTree(np.mod(f.coords, boxsize), boxsize=boxsize)
            # alternative: https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html

            for i, k in enumerate(f.atom2kind):
                if k == kind_name:
                    # Find indicies of neighbor atoms.
                    nearest = kdtree.query(f.coords[i], num_neighbors + 1)[1]
                    neighbors = [j for j in nearest if j != i]  # filter central atom

                    # Compute relative positions of neighbor atoms.
                    relpos = [f.coords[j] - f.coords[i] for j in neighbors]
                    self.neighbors_relpos.append(_as_tensor(relpos))

                    # Features of neighbor atoms is the one-hot encoding of their kind.
                    features = [
                        f.atom2kind[j] == np.array(self.feature_kind_names)
                        for j in neighbors
                    ]
                    self.neighbors_features.append(_as_tensor(features))

                    # Orthonormalize labels as it's required for the loss_functions.
                    label = np.linalg.svd(f.xblocks[i], full_matrices=False)[2]
                    self.labels.append(_as_tensor(label))

    def __len__(self) -> int:
        return len(self.labels)

    def __getitem__(self, i: int) -> PaoRecord:
        return self.neighbors_relpos[i], self.neighbors_features[i], self.labels[i]

    # To make mypy happy.
    def __iter__(self) -> Iterator[PaoRecord]:
        for i in range(len(self)):
            yield self[i]


# EOF
