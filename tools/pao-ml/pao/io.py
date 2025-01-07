# author: Ole Schuett

import re
from pathlib import Path
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple
import numpy as np
import numpy.typing as npt

KindName = str
NDArray = npt.NDArray[np.float64]


# ======================================================================================
@dataclass
class AtomicKind:
    name: str
    atomic_number: int
    nparams: int = -1
    prim_basis_size: int = -1
    prim_basis_name: str = ""
    pao_basis_size: int = -1


# ======================================================================================
@dataclass
class PaoFile:
    kinds: Dict[KindName, AtomicKind]
    atom2kind: List[KindName]
    cell: NDArray
    coords: NDArray
    xblocks: List[NDArray]


# ======================================================================================
def parse_pao_file(path: Path) -> PaoFile:
    assert str(path).endswith(".pao")

    ikind2name = {}  # maps kind index to kind name
    atom2kind: List[KindName] = []  # maps atom index to kind name
    kinds: Dict[KindName, AtomicKind] = {}
    cell: NDArray = np.zeros((3, 3))
    coords_list: List[List[str]] = []
    xblocks: List[NDArray] = []

    for line in path.read_text().strip().split("\n"):
        parts = line.split()
        if parts[0] == "Parametrization":
            assert parts[1] == "EQUIVARIANT"

        elif parts[0] == "Kind":
            ikind = int(parts[1])
            kind = AtomicKind(name=parts[2], atomic_number=int(parts[3]))
            kinds[kind.name] = kind
            ikind2name[ikind] = kind.name

        elif parts[0] == "NParams":
            ikind = int(parts[1])
            kinds[ikind2name[ikind]].nparams = int(parts[2])

        elif parts[0] == "PrimBasis":
            ikind = int(parts[1])
            kinds[ikind2name[ikind]].prim_basis_size = int(parts[2])
            kinds[ikind2name[ikind]].prim_basis_name = parts[3]

        elif parts[0] == "PaoBasis":
            ikind = int(parts[1])
            kinds[ikind2name[ikind]].pao_basis_size = int(parts[2])

        elif parts[0] == "Cell":
            cell = np.array(parts[1:], float).reshape(3, 3)

        elif parts[0] == "Atom":
            atom2kind.append(parts[2])
            coords_list.append(parts[3:])

        elif parts[0] == "Xblock":
            xblocks.append(np.array(parts[2:], float))

    # Convert coordinates to numpy array.
    coords = np.array(coords_list, float)

    # Reshape xblocks.
    for iatom, kind_name in enumerate(atom2kind):
        n = kinds[kind_name].prim_basis_size
        m = kinds[kind_name].pao_basis_size
        xblocks[iatom] = xblocks[iatom].reshape(m, n)

    return PaoFile(kinds, atom2kind, cell, coords, xblocks)


# ======================================================================================
def write_pao_file(
    path: Path,
    kinds: Dict[KindName, AtomicKind],
    atom2kind: List[KindName],
    coords: NDArray,
    xblocks: List[NDArray],
) -> None:
    natoms = coords.shape[0]
    assert coords.shape[1] == 3
    assert len(xblocks) == natoms

    output = []
    output.append("Version 4")
    output.append("Parametrization EQUIVARIANT")
    output.append(f"Nkinds {len(kinds)}")
    for ikind, (kind_name, kind) in enumerate(kinds.items()):
        i = ikind + 1
        output.append(f"Kind {i} {kind_name} {kind.atomic_number}")
        output.append(f"NParams {i} {kind.nparams}")
        output.append(f"PrimBasis {i} {kind.prim_basis_size} {kind.prim_basis_name}")
        output.append(f"PaoBasis {i} {kind.pao_basis_size}")
        output.append(f"NPaoPotentials {i} 0")
    output.append("Cell 8.0 0.0 0.0   0.0 8.0 0.0   0.0 0.0 8.0")
    output.append(f"Natoms {natoms}")

    for iatom in range(natoms):
        c = coords[iatom, :]
        output.append(f"Atom {iatom + 1} {atom2kind[iatom]} {c[0]} {c[1]} {c[2]}")

    for iatom in range(natoms):
        kind = kinds[atom2kind[iatom]]
        assert len(xblocks[iatom].shape) == 2
        assert xblocks[iatom].shape[0] == kind.pao_basis_size
        assert xblocks[iatom].shape[1] == kind.prim_basis_size
        x = xblocks[iatom].flatten()
        y = " ".join(["%f" % i for i in x])
        output.append(f"Xblock {iatom + 1} {y}")

    output.append("THE_END")
    path.write_text("\n".join(output))


# ======================================================================================
def read_cp2k_energy(path: Path) -> float:
    try:
        content = path.read_text()
        m = re.search(r"ENERGY\|(.*)", content)
        assert m
        return float(m.group(1).split()[-1])
    except:
        print(f"error with: {path}")
    return float("NaN")


# EOF
