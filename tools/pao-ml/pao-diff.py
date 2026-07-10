#!/usr/bin/env python3

# author: Ole Schuett

import numpy as np
import argparse
from pathlib import Path

from pao.io import parse_pao_file
from pao.dataset import prepare_xblock
from pao.training import loss_function


# ======================================================================================
def main() -> None:
    description = "Compares two PAO files by computing the loss between their xblocks."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("pao_file_a", type=Path)
    parser.add_argument("pao_file_b", type=Path)
    args = parser.parse_args()

    a = parse_pao_file(args.pao_file_a)
    b = parse_pao_file(args.pao_file_b)

    assert a.kinds == b.kinds
    assert a.atom2kind == b.atom2kind
    assert np.all(a.cell == b.cell)
    assert np.all(a.coords == b.coords)

    for iatom in range(len(a.atom2kind)):
        xa = prepare_xblock(a.xblocks[iatom])
        xb = prepare_xblock(b.xblocks[iatom])
        loss = loss_function(xa, xb)
        print(f"iatom: {iatom+1}  kind: {a.atom2kind[iatom]}  loss: {loss:.8e}")


# ======================================================================================
main()

# EOF
