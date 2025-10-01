# author: Ole Schuett

import torch
import numpy as np
import numpy.typing as npt
from numpy.linalg import norm
from pathlib import Path
from typing import Any, List, Tuple, Iterator
from torch.utils.data import Dataset
from nequip.data import AtomicDataDict, register_fields  # type: ignore

from .io import parse_pao_file

NDArray = npt.NDArray[np.float64]

# ======================================================================================
register_fields(graph_fields=["xblock"])  # no need to register "central_edge_index"

collate_fn = AtomicDataDict.batched_from_list


# ======================================================================================
def _as_float_tensor(x: List[Any] | npt.NDArray[np.float64]) -> torch.Tensor:
    return torch.tensor(np.array(x, dtype=np.float32))


# ======================================================================================
def _as_int_tensor(x: List[Any] | npt.NDArray[np.int64]) -> torch.Tensor:
    return torch.tensor(np.array(x, dtype=np.int64))


# ======================================================================================
class PaoDataset(Dataset[AtomicDataDict]):
    def __init__(
        self,
        kind_name: str,
        num_layers: int,
        cutoff: float,
        files: List[Path],
    ):
        self.examples: List[AtomicDataDict] = []

        # Load kinds from the first training data file.
        kinds = parse_pao_file(files[0]).kinds
        self.kind = kinds[kind_name]
        self.all_kind_names = sorted(kinds.keys())

        # Load all training data files.
        for fn in files:
            f = parse_pao_file(fn)
            natoms = f.coords.shape[0]
            all_cell_shifts = [
                i * f.cell[0, :] + j * f.cell[1, :] + k * f.cell[2, :]
                for i in (-1, 0, +1)
                for j in (-1, 0, +1)
                for k in (-1, 0, +1)
            ]

            for iatom, ikind in enumerate(f.atom2kind):
                if ikind != kind_name:
                    continue

                # Find atoms that are reachable within num_layers * cutoff.
                neighbor_pos: List[NDArray] = [f.coords[iatom]]
                neighbor_atom_types: List[int] = [self.all_kind_names.index(ikind)]
                for jatom, jkind in enumerate(f.atom2kind):
                    for cell_shift in all_cell_shifts:
                        if jatom == iatom and all(cell_shift == 0.0):
                            continue  # central atom is always the first neighbor
                        jpos = f.coords[jatom] + cell_shift
                        if norm(jpos - f.coords[iatom]) < num_layers * cutoff:
                            neighbor_pos.append(jpos)
                            neighbor_atom_types.append(self.all_kind_names.index(jkind))

                # Build connectivity graph of neighbors.
                edge_index: List[Tuple[int, int]] = []
                edge_vectors: List[NDArray] = []
                for ineighbor, ipos in enumerate(neighbor_pos):
                    for jneighbor, jpos in enumerate(neighbor_pos):
                        if ineighbor == jneighbor:
                            continue  # ConvNetLayer already includes self connections
                        edge_vector = jpos - ipos
                        if norm(edge_vector) < cutoff:
                            edge_index.append((ineighbor, jneighbor))
                            edge_vectors.append(edge_vector)

                # TODO remove edges that are more than num_layers hops await from iatom.

                # Orthonormalize labels as it's required for the loss_functions.
                xblock = np.linalg.svd(f.xblocks[iatom], full_matrices=False)[2]

                # Collect all tensors into an AtomicDataDict for NequIP.
                self.examples.append(
                    {
                        "atom_types": _as_int_tensor(neighbor_atom_types),
                        "edge_index": _as_int_tensor(edge_index).T,
                        "edge_vectors": _as_float_tensor(edge_vectors),
                        # Adding an extra leading dimension to simplify the batching.
                        "xblock": _as_float_tensor(xblock).unsqueeze(0),
                        # The "pos" key is used by batched_from_list() to compute edge_index offsets.
                        "pos": torch.zeros(len(neighbor_pos)),
                        # Within an example the central atom is always the first atom.
                        # However, during batching multiple examples are concatenated.
                        # To keep track of the central atom we're using a key that
                        # contains "edge_index" because in batched_from_list() the
                        # correct offsets are added to those fields during batching.
                        "central_edge_index": _as_int_tensor([[0, 0]]).T,
                    }
                )

        # Print some stats.
        avg_neighbors = np.mean([e["atom_types"].shape[0] for e in self.examples])
        avg_edges = np.mean([e["edge_index"].shape[1] for e in self.examples])
        print(
            f"Found {len(self.examples)} examples of kind '{kind_name}' with "
            f"on average {avg_neighbors:.1f} neighbors and {avg_edges:.1f} edges."
        )

    # ----------------------------------------------------------------------------------
    def __len__(self) -> int:
        return len(self.examples)

    def __getitem__(self, i: int) -> AtomicDataDict:
        return self.examples[i]

    # To make mypy happy.
    def __iter__(self) -> Iterator[AtomicDataDict]:
        for i in range(len(self)):
            yield self[i]


# EOF
