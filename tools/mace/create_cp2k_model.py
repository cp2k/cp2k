# pylint: disable=wrong-import-position
"""
Export a trained MACE model to a TorchScript file that can be loaded by CP2K.

The export script must be executed in an environment where MACE is installed.

Required packages:
  * torch: required for loading and converting the model.
  * e3nn: required by MACE for JIT compilation.
  * mace-torch: required because the trained MACE model is loaded using
    torch and its original model classes must be available.

The MACE version used for export should be compatible with the version
used for training the model.

After exporting, the CP2K-compatible `.pth` file can be loaded by CP2K
through the LibTorch interface. No MACE, e3nn, or Python installation is
required during CP2K simulations.

Usage:
    python create_cp2k_model.py my_mace.model --dtype float64 --head default
    -> writes  my_mace.model-cp2k.pth
"""

import argparse
import os
from typing import Any, Dict, Optional

os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"

import torch
from e3nn.util import jit  # type: ignore
from e3nn.util.jit import compile_mode  # type: ignore

# Z -> chemical symbol table (index == atomic number). Covers the dummy X plus
# all 118 known elements, matching CP2K's src/common/periodic_table.F.
chemical_symbols = ["X"] + """
    H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co
    Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb
    Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os
    Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm
    Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
""".split()


@compile_mode("script")
class CP2K_MACE(torch.nn.Module):
    """MACE model wrapped for CP2K's single-dictionary libtorch interface."""

    head: torch.Tensor

    def __init__(self, model: Any, head: Optional[str] = None) -> None:
        super().__init__()
        self.model = model
        self.register_buffer("atomic_numbers", model.atomic_numbers)
        self.register_buffer("r_max", model.r_max)
        self.register_buffer("num_interactions", model.num_interactions)
        self.num_types: int = int(model.atomic_numbers.shape[0])

        if not hasattr(model, "heads"):
            model.heads = [None]
        if head is not None and head in model.heads:
            head_idx = model.heads.index(head)
        else:
            head_idx = len(model.heads) - 1
        self.register_buffer(
            "head", torch.tensor(head_idx, dtype=torch.long).unsqueeze(0)
        )

        for param in self.model.parameters():
            param.requires_grad = False

    def forward(
        self, data: Dict[str, torch.Tensor]
    ) -> Dict[str, Optional[torch.Tensor]]:
        pos = data["pos"]
        edge_index = data["edge_index"]
        unit_shifts = data["edge_cell_shift"]
        cell = data["cell"].view(3, 3)
        atom_types = data["atom_types"].to(torch.long)
        n_atoms = pos.shape[0]

        # one-hot node attributes expected by MACE
        node_attrs = torch.nn.functional.one_hot(
            atom_types, num_classes=self.num_types
        ).to(pos.dtype)

        # real-space edge shift vectors: unit_shifts @ cell  (Angstrom)
        shifts = torch.matmul(unit_shifts, cell)

        batch = torch.zeros(n_atoms, dtype=torch.long, device=pos.device)
        ptr = torch.tensor([0, n_atoms], dtype=torch.long, device=pos.device)

        mace_in: Dict[str, torch.Tensor] = {
            "positions": pos,
            "node_attrs": node_attrs,
            "edge_index": edge_index,
            "shifts": shifts,
            "unit_shifts": unit_shifts,
            "cell": cell,
            "batch": batch,
            "ptr": ptr,
            "head": self.head,
        }

        out = self.model(
            mace_in,
            training=False,
            compute_force=True,
            compute_virials=True,
            compute_stress=False,
            compute_displacement=True,
        )

        node_energy = out["node_energy"]
        forces = out["forces"]
        virials = out["virials"]

        atomic_energy: Optional[torch.Tensor] = None
        if node_energy is not None:
            atomic_energy = node_energy.unsqueeze(-1)

        virial: Optional[torch.Tensor] = None
        if virials is not None:
            virial = virials.view(3, 3)

        return {
            "atomic_energy": atomic_energy,
            "forces": forces,
            "virial": virial,
        }


def build_metadata(model: Any, dtype: str) -> Dict[str, str]:
    z = model.atomic_numbers.to(torch.long).tolist()
    type_names = " ".join(chemical_symbols[int(zi)] for zi in z)
    return {
        "num_types": str(len(z)),
        "r_max": repr(float(model.r_max)),
        "type_names": type_names,
        "model_dtype": dtype,
        "allow_tf32": "0",
        # MACE uses a single global cutoff; leave per-edge-type cutoffs empty so
        # CP2K's reader falls back to filling the cutoff matrix with r_max**2.
        "per_edge_type_cutoff": "",
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("model_path", help="Path to the trained MACE .model file")
    p.add_argument(
        "--head", default=None, help="Model head to export (multi-head models)"
    )
    p.add_argument("--dtype", choices=["float64", "float32"], default="float64")
    p.add_argument(
        "--output", default=None, help="Output path (default: <model>-cp2k.pth)"
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    model = torch.load(args.model_path, map_location="cpu")
    if args.dtype == "float64":
        model = model.double()
    else:
        print("Converting model to float32, this may cause loss of precision.")
        model = model.float()
    model = model.to("cpu")

    wrapper = CP2K_MACE(model, head=args.head)
    wrapper.eval()

    scripted = jit.compile(wrapper)

    metadata = build_metadata(model, args.dtype)
    extra_files = {k: v.encode("utf-8") for k, v in metadata.items()}

    out_path = args.output or (args.model_path + "-cp2k.pth")
    scripted.save(out_path, _extra_files=extra_files)

    print(f"Wrote CP2K MACE model to: {out_path}")
    print("Metadata:")
    for k, v in metadata.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
