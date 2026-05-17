#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compare Wannier90 .mmn files block by block.

The raw overlap matrices depend on the band gauge at each k-point. Their
singular values are invariant under unitary rotations within equally sized band
subspaces and are therefore useful when validating reconstructed k-point MOs
against a full-mesh reference.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence, TypeAlias

import numpy as np
from numpy.typing import NDArray

Header: TypeAlias = tuple[int, int, int, int, int]
MmnMatrix: TypeAlias = NDArray[np.complex128]
MmnBlock: TypeAlias = tuple[Header, MmnMatrix]
MmnData: TypeAlias = tuple[int, int, int, list[MmnBlock]]
WorstBlock: TypeAlias = tuple[int, Header] | None


def parse_header(line: str, filename: Path) -> Header:
    fields = line.split()
    if len(fields) < 5:
        raise ValueError(f"{filename}: incomplete block header")
    try:
        return (
            int(fields[0]),
            int(fields[1]),
            int(fields[2]),
            int(fields[3]),
            int(fields[4]),
        )
    except ValueError as exc:
        raise ValueError(f"{filename}: could not parse block header") from exc


def parse_complex(line: str, filename: Path) -> complex:
    fields = line.split()
    if len(fields) < 2:
        raise ValueError(f"{filename}: incomplete block value")
    try:
        real, imag = (float(x) for x in fields[:2])
    except ValueError as exc:
        raise ValueError(f"{filename}: could not parse block value") from exc
    return complex(real, imag)


def read_mmn(filename: str) -> MmnData:
    path = Path(filename)
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"{filename}: not a Wannier90 .mmn file")
    try:
        nbands, nkpts, nneigh = [int(x) for x in lines[1].split()[:3]]
    except ValueError as exc:
        raise ValueError(f"{filename}: could not parse .mmn dimensions") from exc

    blocks: list[MmnBlock] = []
    idx = 2
    for _ in range(nkpts * nneigh):
        if idx >= len(lines):
            raise ValueError(f"{filename}: missing block header")
        header = parse_header(lines[idx], path)
        idx += 1
        matrix: MmnMatrix = np.zeros((nbands, nbands), dtype=np.complex128)
        for iband in range(nbands):
            for jband in range(nbands):
                if idx >= len(lines):
                    raise ValueError(f"{filename}: missing block values")
                value = parse_complex(lines[idx], path)
                idx += 1
                matrix[jband, iband] = value
        blocks.append((header, matrix))
    return nbands, nkpts, nneigh, blocks


def compare_mmn(
    left: MmnData, right: MmnData
) -> tuple[float, float, WorstBlock, WorstBlock]:
    left_dims = left[:3]
    right_dims = right[:3]
    if left_dims != right_dims:
        raise ValueError(f"incompatible .mmn dimensions: {left_dims} != {right_dims}")

    raw_max = 0.0
    sv_max = 0.0
    worst_raw: WorstBlock = None
    worst_sv: WorstBlock = None
    for iblock, ((left_header, left_matrix), (right_header, right_matrix)) in enumerate(
        zip(left[3], right[3]), start=1
    ):
        if left_header != right_header:
            raise ValueError(
                f"block {iblock}: incompatible headers {left_header} != {right_header}"
            )
        raw = float(np.max(np.abs(left_matrix - right_matrix)))
        left_sv = np.linalg.svd(left_matrix, compute_uv=False)
        right_sv = np.linalg.svd(right_matrix, compute_uv=False)
        sv = float(np.max(np.abs(left_sv - right_sv)))
        if raw > raw_max:
            raw_max = raw
            worst_raw = (iblock, left_header)
        if sv > sv_max:
            sv_max = sv
            worst_sv = (iblock, left_header)
    return raw_max, sv_max, worst_raw, worst_sv


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left", help="first Wannier90 .mmn file")
    parser.add_argument("right", help="second Wannier90 .mmn file")
    parser.add_argument("--raw-tol", type=float, default=None)
    parser.add_argument("--sv-tol", type=float, default=None)
    args = parser.parse_args(argv)

    raw_max, sv_max, worst_raw, worst_sv = compare_mmn(
        read_mmn(args.left), read_mmn(args.right)
    )
    print(f"raw_max = {raw_max:.16e}")
    print(f"sv_max  = {sv_max:.16e}")
    if worst_raw is not None:
        print(f"worst_raw_block = {worst_raw[0]} header={worst_raw[1]}")
    if worst_sv is not None:
        print(f"worst_sv_block  = {worst_sv[0]} header={worst_sv[1]}")

    failed = False
    if args.raw_tol is not None and raw_max > args.raw_tol:
        failed = True
    if args.sv_tol is not None and sv_max > args.sv_tol:
        failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
