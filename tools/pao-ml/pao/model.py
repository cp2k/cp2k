# author: Ole Schuett

import torch
from typing import Any, List, Dict

from e3nn import o3  # type: ignore
from e3nn.o3._linear import Linear  # type: ignore
from nequip.data import AtomicDataDict  # type: ignore
from nequip.nn import GraphModel, SequentialGraphNetwork, ConvNetLayer, GraphModuleMixin  # type: ignore
from nequip.nn.embedding import (  # type: ignore
    NodeTypeEmbed,
    PolynomialCutoff,
    EdgeLengthNormalizer,
    BesselEdgeLengthEncoding,
    SphericalHarmonicEdgeAttrs,
)

# from nequip.model import model_builder


# ======================================================================================
# Based on NequIPGNNEnergyModel:
# https://github.com/mir-group/nequip/blob/main/nequip/model/nequip_models.py
# @model_builder  # TODO: Probably needed for PyTorch's Ahead-of-Time Inductor compiler.
class PaoModel(SequentialGraphNetwork):  # type: ignore

    # Ensure these get picked up as attributes by TorchScript.
    # https://pytorch.org/docs/stable/generated/torch.jit.Attribute.html
    pao_model_version: int
    kind_name: str
    atomic_number: int
    prim_basis_name: str
    prim_basis_size: int
    pao_basis_size: int
    all_kind_names: List[str]
    num_layers: int
    cutoff: float

    def __init__(
        self,
        kind_name: str,
        atomic_number: int,
        prim_basis_name: str,
        prim_basis_size: int,
        pao_basis_size: int,
        all_kind_names: List[str],
        avg_num_neighbors: float,
        num_layers: int,
        cutoff: float,
        # convnet params
        parity: bool = True,  # TODO really needed?
        num_features: int = 32,
        radial_mlp_depth: int = 2,
        radial_mlp_width: int = 64,
        # edge length encoding
        num_bessels: int = 8,
        bessel_trainable: bool = False,
        polynomial_cutoff_p: int = 6,
    ):
        assert num_layers > 0
        assert all(tn.isalnum() for tn in all_kind_names)

        self.pao_model_version = 2
        self.kind_name = kind_name
        self.atomic_number = atomic_number
        self.prim_basis_name = prim_basis_name
        self.prim_basis_size = prim_basis_size
        self.pao_basis_size = pao_basis_size
        self.all_kind_names = all_kind_names
        self.num_layers = num_layers
        self.cutoff = cutoff

        # Irreps of primary basis
        # TODO: Export the specs directly from cp2k as part of the .pao files.
        prim_basis_specs = {
            "DZVP-MOLOPT-GTH/H": "2x0e + 1x1o",  # two s-shells, one p-shell
            "DZVP-MOLOPT-GTH/O": "2x0e + 2x1o + 1x2e",  # two s, two p, one d-shell
            "TZV2P-MOLOPT-GGA-GTH-q1/H": "3x0e + 2x1o + 1x2e",
            "TZV2P-MOLOPT-GGA-GTH-q6/O": "3x0e + 3x1o + 2x2e + 1x3o",
        }
        basis_specs_key = f"{prim_basis_name}/{kind_name}"
        irreps_prim_basis = o3.Irreps(prim_basis_specs[basis_specs_key])
        assert prim_basis_size == irreps_prim_basis.dim

        # Spherical harmonics
        lmax = 2 * irreps_prim_basis.lmax
        irreps_edge_sh = repr(
            o3.Irreps.spherical_harmonics(lmax=lmax, p=-1 if parity else 1)
        )

        # Convert a single set of parameters uniformly for every layer
        feature_irreps_hidden = repr(
            o3.Irreps(
                [
                    (num_features, (l, p))
                    for p in ((1, -1) if parity else (1,))
                    for l in range(lmax + 1)
                ]
            )
        )

        # Assemble modules for SequentialGraphNetwork.
        modules: Dict[str, GraphModuleMixin] = dict()

        # Edge tensor embedding
        modules["spharm"] = SphericalHarmonicEdgeAttrs(irreps_edge_sh=irreps_edge_sh)

        # Edge scalar embedding
        modules["edge_norm"] = EdgeLengthNormalizer(
            r_max=cutoff,
            type_names=all_kind_names,
            irreps_in=modules["spharm"].irreps_out,
        )

        # Edge length encoding
        modules["bessel_encode"] = BesselEdgeLengthEncoding(
            num_bessels=num_bessels,
            trainable=bessel_trainable,
            cutoff=PolynomialCutoff(polynomial_cutoff_p),
            irreps_in=modules["edge_norm"].irreps_out,
        )

        # Node scalar embedding
        modules["type_embed"] = NodeTypeEmbed(
            type_names=all_kind_names,
            num_features=num_features,
            irreps_in=modules["bessel_encode"].irreps_out,
        )

        # ConvNetLayer layers
        prev_irreps_out = modules["type_embed"].irreps_out
        for ilayer in range(num_layers):
            modules[f"convnet_layer_{ilayer}"] = ConvNetLayer(
                irreps_in=prev_irreps_out,
                feature_irreps_hidden=feature_irreps_hidden,
                convolution_kwargs={
                    "radial_mlp_depth": radial_mlp_depth,
                    "radial_mlp_width": radial_mlp_width,
                    "avg_num_neighbors": avg_num_neighbors,
                    # to ensure isolated atom limit
                    # "use_sc": layer_i != 0,  # TODO user self connection
                },
            )
            prev_irreps_out = modules[f"convnet_layer_{ilayer}"].irreps_out

        # Readout into PAO basis
        modules["pao_readout"] = PaoReadout(
            irreps_in=prev_irreps_out,
            irreps_prim_basis=irreps_prim_basis,
            pao_basis_size=pao_basis_size,
        )

        super().__init__(modules)


# ======================================================================================
class PaoReadout(GraphModuleMixin, torch.nn.Module):  # type: ignore
    def __init__(
        self, irreps_in: Any, irreps_prim_basis: Any, pao_basis_size: int
    ) -> None:
        super().__init__()
        self._init_irreps(  # type: ignore
            irreps_in=irreps_in,
            required_irreps_in=["node_features"],
            irreps_out={},  # no irreps to output
        )
        self.pao_basis_size = pao_basis_size
        self.matrix = SymmetricMatrix(irreps_prim_basis)  # auxiliary Hamiltonian
        self.linear = Linear(
            irreps_in=self.irreps_in["node_features"], irreps_out=self.matrix.irreps_in  # type: ignore
        )
        # CP2K uses the yzx convention, while e3nn uses xyz.
        # https://docs.e3nn.org/en/stable/guide/change_of_basis.html
        yzx_to_xyz = torch.tensor([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        self.D_yzx_to_xyz = irreps_prim_basis.D_from_matrix(yzx_to_xyz)

    def forward(self, data: AtomicDataDict.Type) -> AtomicDataDict.Type:
        # Due to batching there can be multiple central atoms.
        central_atoms = data["central_edge_index"][0]
        central_atom_features = data["node_features"].index_select(0, central_atoms)
        h_aux_vec = self.linear(central_atom_features)
        h_aux_matrix = self.matrix(h_aux_vec)
        u_matrix = torch.linalg.eigh(h_aux_matrix)[1]
        xblock = u_matrix[..., : self.pao_basis_size].transpose(-2, -1)
        outputs = {"xblock": xblock @ self.D_yzx_to_xyz}
        return outputs


# ======================================================================================
def flatten_irreps(irreps: o3.Irreps) -> List[o3.Irrep]:
    "Helper function to turn an Irreps object into a list of individual Irrep objects."
    result = []
    for mul, ir in irreps:
        result += mul * [ir]
    return result


# ======================================================================================
def dim(l: int) -> int:
    return 2 * l + 1


# ======================================================================================
class SymmetricMatrix(torch.nn.Module):
    def __init__(self, irreps_prim_basis: Any):
        super().__init__()
        self.irreps_prim_basis = irreps_prim_basis
        self.irreps_prim_basis_ls: List[int] = irreps_prim_basis.ls

        # Compute irreps required to represent a matrix
        self.irreps_in = o3.Irreps()
        self.wigner_blocks = []
        for i, a in enumerate(flatten_irreps(irreps_prim_basis)):
            for j, b in enumerate(flatten_irreps(irreps_prim_basis)):
                if j > i:
                    continue  # skip upper triangle
                self.irreps_in += a * b

                # Pre-compute wigner blocks
                for lk in range(abs(a.l - b.l), a.l + b.l + 1):
                    self.wigner_blocks.append(o3.wigner_3j(a.l, b.l, lk))

        self.irreps_in_ls = self.irreps_in.ls

    # ----------------------------------------------------------------------------------
    def forward(self, vector: torch.Tensor) -> torch.Tensor:
        assert vector.shape[-1] == sum(dim(l) for l in self.irreps_in_ls)
        basis_size = sum(dim(l) for l in self.irreps_prim_basis_ls)
        matrix = torch.zeros(vector.shape[:-1] + (basis_size, basis_size))
        matrix[..., :, :] = torch.eye(basis_size)  # ensure matrix is diagonalizable
        c = 0  # position in vector
        z = 0  # position in self.wigner_blocks
        for i, li in enumerate(self.irreps_prim_basis_ls):
            a = sum(dim(l) for l in self.irreps_prim_basis_ls[:i])  # first matrix row
            for j, lj in enumerate(self.irreps_prim_basis_ls):
                if j > i:
                    continue  # skip upper triangle
                b = sum(dim(l) for l in self.irreps_prim_basis_ls[:j])  # 1st matrix col
                for lk in range(abs(li - lj), li + lj + 1):
                    # TODO the wigner blocks are mostly zeros - not sure pytorch takes advantage of that.
                    wigner_block = self.wigner_blocks[z]
                    coeffs = vector[..., c : c + dim(lk)]
                    block = torch.tensordot(coeffs, wigner_block, dims=[[-1], [-1]])
                    matrix[..., a : a + dim(li), b : b + dim(lj)] += block
                    c += dim(lk)
                    z += 1

        return matrix


# EOF
