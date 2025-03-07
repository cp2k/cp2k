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
            assert 0 < num_neighbors
            assert np.all(f.cell == np.diag(np.diagonal(f.cell)))
            boxsize = np.diagonal(f.cell)
            kdtree = scipy.spatial.KDTree(np.mod(f.coords, boxsize), boxsize=boxsize)
            # alternative: https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html

            # Small systems might contain less atoms than the model's num_neighbors.
            num_available_neighbors = min(f.coords.shape[0] - 1, num_neighbors)
            for iatom, ikind in enumerate(f.atom2kind):
                if ikind != kind_name:
                    continue

                # Find indicies of neighbor atoms.
                nearest = kdtree.query(f.coords[iatom], num_available_neighbors + 1)[1]
                assert nearest[0] == iatom  # the central atom should be the nearest

                # Compute relative positions and features of neighbor atoms.
                relpos = np.zeros([num_neighbors, 3])
                features = np.zeros([num_neighbors, len(self.feature_kind_names)])
                for jneighbor in range(num_available_neighbors):
                    jatom = nearest[jneighbor + 1]  # +1 to skip over central atom
                    relpos[jneighbor, :] = f.coords[jatom] - f.coords[iatom]
                    # Features of neighbor atoms is the one-hot encoding of their kind.
                    one_hot = f.atom2kind[jatom] == np.array(self.feature_kind_names)
                    features[jneighbor, :] = one_hot
                self.neighbors_relpos.append(_as_tensor(relpos))
                self.neighbors_features.append(_as_tensor(features))

                # Orthonormalize labels as it's required for the loss_functions.
                label = np.linalg.svd(f.xblocks[iatom], full_matrices=False)[2]
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
