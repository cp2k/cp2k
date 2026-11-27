#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
#   CP2K: A general program to perform molecular dynamics simulations
#   Copyright 2000-2026 CP2K developers group <https://cp2k.org>
#
#   SPDX-License-Identifier: GPL-2.0-or-later
# --------------------------------------------------------------------------------------------------

"""Unit tests for tools/band_structure/plot_band_structure.py."""

from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np

from plot_band_structure import (
    EV_PER_HARTREE,
    compute_kpath_geometry,
    compute_spectral_map,
    export_csv,
    main,
    parse_bs_file,
    parse_cp2k_out,
)

SAMPLE_STD_BS = """# Set 1: 2 special points, 3 k-points, 2 bands
#  Special point 1   0.00000000   0.00000000   0.00000000  GAMMA
#  Special point 2   0.50000000   0.00000000   0.00000000  X
#  Point 1  Spin 1:   0.00000000   0.00000000   0.00000000  0.33333333
#   Band    Energy [eV]     Occupation
       1    -5.00000000     2.00000000
       2    -1.00000000     0.00000000
#  Point 2  Spin 1:   0.25000000   0.00000000   0.00000000  0.33333333
#   Band    Energy [eV]     Occupation
       1    -4.50000000     2.00000000
       2    -1.20000000     0.00000000
#  Point 3  Spin 1:   0.50000000   0.00000000   0.00000000  0.33333333
#   Band    Energy [eV]     Occupation
       1    -4.00000000     2.00000000
       2    -1.50000000     0.00000000
# Set 2: 2 special points, 2 k-points, 2 bands
#  Special point 1   0.33333333   0.33333333   0.00000000  K
#  Special point 2   0.00000000   0.00000000   0.00000000  GAMMA
#  Point 1  Spin 1:   0.33333333   0.33333333   0.00000000  0.50000000
#   Band    Energy [eV]     Occupation
       1    -4.20000000     2.00000000
       2    -1.40000000     0.00000000
#  Point 2  Spin 1:   0.00000000   0.00000000   0.00000000  0.50000000
#   Band    Energy [eV]     Occupation
       1    -5.00000000     2.00000000
       2    -1.00000000     0.00000000
"""

SAMPLE_CP2K_OUT = """
 CELL| Vector a [angstrom]:       3.155     0.000     0.000   |a| =     3.155000
 CELL| Vector b [angstrom]:      -1.577     2.732     0.000   |b| =     3.155000
 CELL| Vector c [angstrom]:       0.000     0.000    18.000   |c| =    18.000000
  Fermi energy:                                                -0.15000000000000
  Fermi energy:                                                -0.14700000000000
"""

SAMPLE_UNFOLD_BS = """# UNFOLD_BAND| primitive cell vector 1 [Bohr] =    5.962086    0.000000    0.000000
# UNFOLD_BAND| primitive cell vector 2 [Bohr] =   -2.981043    5.163318    0.000000
# UNFOLD_BAND| primitive cell vector 3 [Bohr] =    0.000000    0.000000   34.015070
# UNFOLD_BAND| atom groups: 2 (columns 5..6 = per-group weights)
# Set 1: 2 special points, 3 k-points, 2 bands
#  Special point 1   0.50000000   0.00000000   0.00000000  M
#  Special point 2   0.00000000   0.00000000   0.00000000  GAMMA
#  Point 1  Spin 1:   0.50000000   0.00000000   0.00000000  0.33333333
#   Band    Energy [eV]     Occupation    Unfold_W      W_grp1      W_grp2
       1    -2.00000000     1.00000000    0.80000000  0.60000000  0.20000000
       2     1.00000000     0.00000000    0.50000000  0.40000000  0.10000000
#  Point 1  Spin 2:   0.50000000   0.00000000   0.00000000  0.33333333
#   Band    Energy [eV]     Occupation    Unfold_W      W_grp1      W_grp2
       1    -1.90000000     1.00000000    0.75000000  0.55000000  0.20000000
       2     1.10000000     0.00000000    0.45000000  0.35000000  0.10000000
#  Point 2  Spin 1:   0.25000000   0.00000000   0.00000000  0.33333333
#   Band    Energy [eV]     Occupation    Unfold_W      W_grp1      W_grp2
       1    -1.50000000     1.00000000    0.85000000  0.65000000  0.20000000
       2     0.80000000     0.00000000    0.55000000  0.45000000  0.10000000
#  Point 2  Spin 2:   0.25000000   0.00000000   0.00000000  0.33333333
#   Band    Energy [eV]     Occupation    Unfold_W      W_grp1      W_grp2
       1    -1.45000000     1.00000000    0.80000000  0.60000000  0.20000000
       2     0.85000000     0.00000000    0.50000000  0.40000000  0.10000000
#  Point 3  Spin 1:   0.00000000   0.00000000   0.00000000  0.33333333
#   Band    Energy [eV]     Occupation    Unfold_W      W_grp1      W_grp2
       1    -1.00000000     1.00000000    0.90000000  0.70000000  0.20000000
       2     0.50000000     0.00000000    0.60000000  0.50000000  0.10000000
#  Point 3  Spin 2:   0.00000000   0.00000000   0.00000000  0.33333333
#   Band    Energy [eV]     Occupation    Unfold_W      W_grp1      W_grp2
       1    -0.95000000     1.00000000    0.85000000  0.65000000  0.20000000
       2     0.55000000     0.00000000    0.55000000  0.45000000  0.10000000
"""


class TestBandStructureTool(unittest.TestCase):
    def test_standard_multiset_and_cp2k_out(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bs_path = root / "std.bs"
            out_path = root / "std.out"
            csv_path = root / "std.csv"
            bs_path.write_text(SAMPLE_STD_BS, encoding="utf-8")
            out_path.write_text(SAMPLE_CP2K_OUT, encoding="utf-8")

            bs = parse_bs_file(bs_path)
            self.assertFalse(bs.is_unfolded)
            self.assertEqual(len(bs.sets), 2)
            self.assertEqual(bs.spins, [1])

            out_info = parse_cp2k_out(out_path)
            self.assertIsNotNone(out_info.fermi_ev)
            assert out_info.fermi_ev is not None
            self.assertAlmostEqual(out_info.fermi_ev, -0.147 * EV_PER_HARTREE, places=6)

            geom = compute_kpath_geometry(bs, out_info.cell_angstrom)
            self.assertEqual(geom.tick_labels, ["\u0393", "X|K", "\u0393"])
            self.assertEqual(len(geom.break_positions), 1)

            export_csv(bs, geom, energy_zero=out_info.fermi_ev, out_path=csv_path)
            with csv_path.open("r", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            self.assertEqual(
                rows[0],
                [
                    "set",
                    "k_index",
                    "k_dist",
                    "kx",
                    "ky",
                    "kz",
                    "spin",
                    "band",
                    "energy_ev",
                    "occupation",
                ],
            )
            self.assertEqual(len(rows) - 1, 10)

    def test_unfolded_uks_multigroup_and_spectrum(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bs_path = root / "unfold.bs"
            csv_path = root / "unfold.csv"
            bs_path.write_text(SAMPLE_UNFOLD_BS, encoding="utf-8")

            bs = parse_bs_file(bs_path)
            self.assertTrue(bs.is_unfolded)
            self.assertEqual(bs.n_groups, 2)
            self.assertEqual(bs.spins, [1, 2])
            self.assertIsNotNone(bs.prim_cell_angstrom)

            geom = compute_kpath_geometry(bs, None)
            self.assertEqual(geom.tick_labels, ["M", "\u0393"])

            e_grid = np.linspace(-3.0, 2.0, 201, dtype=np.float64)
            pts_s1 = bs.sets[0].points_by_spin[1]
            spec_tot = compute_spectral_map(pts_s1, 0.0, e_grid, 0.05, None)
            spec_g1 = compute_spectral_map(pts_s1, 0.0, e_grid, 0.05, 1)
            spec_g2 = compute_spectral_map(pts_s1, 0.0, e_grid, 0.05, 2)
            self.assertTrue(np.allclose(spec_tot, spec_g1 + spec_g2, atol=1.0e-12))

            export_csv(bs, geom, energy_zero=0.0, out_path=csv_path)
            with csv_path.open("r", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            self.assertEqual(rows[0][-3:], ["weight", "w_grp1", "w_grp2"])
            self.assertEqual(len(rows) - 1, 12)

    def test_appended_bs_uses_last_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bs_path = Path(tmp) / "appended.bs"
            bs_path.write_text(
                SAMPLE_STD_BS + "\n" + SAMPLE_UNFOLD_BS, encoding="utf-8"
            )
            bs = parse_bs_file(bs_path)
            self.assertTrue(bs.is_unfolded)
            self.assertEqual(len(bs.sets), 1)
            self.assertEqual(bs.n_groups, 2)

    def test_plotting_modes_if_matplotlib_installed(self) -> None:
        if importlib.util.find_spec("matplotlib") is None:
            self.skipTest("matplotlib is not installed")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            std_bs = root / "std.bs"
            unf_bs = root / "unf.bs"
            out_f = root / "run.out"
            std_bs.write_text(SAMPLE_STD_BS, encoding="utf-8")
            unf_bs.write_text(SAMPLE_UNFOLD_BS, encoding="utf-8")
            out_f.write_text(SAMPLE_CP2K_OUT, encoding="utf-8")

            for argv, out_img in [
                (
                    [
                        str(std_bs),
                        "--cp2k-out",
                        str(out_f),
                        "-o",
                        str(root / "std.png"),
                    ],
                    root / "std.png",
                ),
                (
                    [
                        str(unf_bs),
                        "--mode",
                        "spectrum",
                        "--spin",
                        "1",
                        "-o",
                        str(root / "spec.png"),
                    ],
                    root / "spec.png",
                ),
                (
                    [
                        str(unf_bs),
                        "--mode",
                        "lines",
                        "--spin",
                        "1",
                        "--group",
                        "2",
                        "--dot-stride",
                        "1",
                        "--cmap",
                        "RdYlBu_r",
                        "-o",
                        str(root / "lines.png"),
                    ],
                    root / "lines.png",
                ),
            ]:
                rc = main(argv)
                self.assertEqual(rc, 0)
                self.assertTrue(out_img.is_file() and out_img.stat().st_size > 0)


if __name__ == "__main__":
    unittest.main()
