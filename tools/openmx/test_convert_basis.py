# SPDX-License-Identifier: GPL-2.0-or-later
"""Offline tests; no OpenMX installation or network required."""

import contextlib
import io
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.special import lpmv

from convert_basis import (
    ANGULAR_PERMUTATIONS,
    BasisSpec,
    PaoFile,
    fit_channel,
    gaussian_overlap,
    main,
    primitives,
    quadrature,
    render_basis,
)


def fixture(path, lmax=3, multiplicity=2):
    r = np.geomspace(5e-4, 12.0, 220)
    lines = [
        "AtomSpecies 6",
        "PAO.Lmax " + str(lmax),
        "PAO.Mul " + str(multiplicity),
        "grid.num.output " + str(len(r)),
        "radial.cutoff.pao 10.0",
        "valence.electron 4",
    ]
    for l in range(lmax + 1):
        g, _ = primitives(r, np.array([0.4, 1.2])[:multiplicity], l)
        lines.append(f"<pseudo.atomic.orbitals.L={l}")
        for rr, row in zip(r, g):
            lines.append(" ".join(f"{v:.16e}" for v in [np.log(rr), rr, *row]))
        lines.append(f"pseudo.atomic.orbitals.L={l}>")
    path.write_text("\n".join(lines) + "\n")


def cp2k_harmonic(l, m, xyz):
    # Independent port of CP2K's Cartesian-to-spherical formula, not the
    # permutation table being tested. Truncation of Fortran integer division
    # matters for negative odd phase exponents.
    def choose(n, k):
        return math.comb(n, k) if 0 <= k <= n else 0

    def dfac(n):
        return math.prod(range(n, 0, -2))

    value = np.zeros(len(xyz))
    ma = abs(m)
    for lx in range(l + 1):
        for ly in range(l - lx + 1):
            lz = l - lx - ly
            j = lx + ly - ma
            if j < 0 or j % 2:
                continue
            j //= 2
            s1 = 0.0
            for i in range((l - ma) // 2 + 1):
                s2 = 0.0
                for k in range(j + 1):
                    if (m < 0 and abs(ma - lx) % 2 == 1) or (
                        m > 0 and abs(ma - lx) % 2 == 0
                    ):
                        s = (-1.0) ** int((ma - lx + 2 * k) / 2) * np.sqrt(2)
                    elif m == 0 and lx % 2 == 0:
                        s = (-1.0) ** (k - lx // 2)
                    else:
                        s = 0.0
                    s2 += choose(j, k) * choose(ma, lx - 2 * k) * s
                s1 += (
                    choose(l, i)
                    * choose(i, j)
                    * (-1.0) ** i
                    * math.factorial(2 * l - 2 * i)
                    / math.factorial(l - ma - 2 * i)
                    * s2
                )
            fac = math.factorial
            c2s = (
                np.sqrt(
                    fac(2 * lx)
                    * fac(2 * ly)
                    * fac(2 * lz)
                    * fac(l)
                    * fac(l - ma)
                    / (fac(lx) * fac(ly) * fac(lz) * fac(2 * l) * fac(l + ma))
                )
                * s1
                / (2**l * fac(l))
            )
            slm = (
                np.sqrt(dfac(2 * l + 1) / (4 * np.pi))
                * c2s
                / np.sqrt(dfac(2 * lx - 1) * dfac(2 * ly - 1) * dfac(2 * lz - 1))
            )
            value += slm * xyz[:, 0] ** lx * xyz[:, 1] ** ly * xyz[:, 2] ** lz
    return value


def openmx_harmonics(l, xyz):
    # Get_Orbitals.c evaluated on the unit sphere; coefficients written in
    # exact analytic form to avoid its rounded decimal constants.
    x, y, z = xyz.T
    pi = np.pi
    if l == 0:
        return np.ones((len(x), 1)) / np.sqrt(4 * pi)
    if l == 1:
        return np.array([x, y, z]).T * np.sqrt(3 / (4 * pi))
    if l == 2:
        return np.array(
            [
                np.sqrt(5 / (16 * pi)) * (3 * z * z - 1),
                np.sqrt(15 / (16 * pi)) * (x * x - y * y),
                np.sqrt(15 / (4 * pi)) * x * y,
                np.sqrt(15 / (4 * pi)) * x * z,
                np.sqrt(15 / (4 * pi)) * y * z,
            ]
        ).T
    if l >= 4:
        # Independent ComplexSH + Set_Comp2Real construction from OpenMX 3.9.
        phi = np.arctan2(y, x)
        harmonics = np.zeros((len(x), 2 * l + 1), dtype=complex)
        for m in range(l + 1):
            positive = (
                np.sqrt(
                    (2 * l + 1)
                    / (4 * pi)
                    * math.factorial(l - m)
                    / math.factorial(l + m)
                )
                * lpmv(m, l, z)
                * np.exp(1j * m * phi)
            )
            harmonics[:, l + m] = positive
            harmonics[:, l - m] = (-1) ** m * positive.conj()
        transform = np.zeros((2 * l + 1, 2 * l + 1), dtype=complex)
        transform[0, l] = 1
        for m in range(1, l + 1):
            transform[2 * m - 1, l - m] = 1 / np.sqrt(2)
            transform[2 * m - 1, l + m] = (-1) ** m / np.sqrt(2)
            transform[2 * m, l - m] = 1j / np.sqrt(2)
            transform[2 * m, l + m] = -1j * (-1) ** m / np.sqrt(2)
        result = harmonics @ transform.T
        np.testing.assert_allclose(result.imag, 0, atol=1e-15)
        return result.real
    return np.array(
        [
            np.sqrt(7 / (16 * pi)) * z * (5 * z * z - 3),
            np.sqrt(21 / (32 * pi)) * x * (5 * z * z - 1),
            np.sqrt(21 / (32 * pi)) * y * (5 * z * z - 1),
            np.sqrt(105 / (16 * pi)) * z * (x * x - y * y),
            np.sqrt(105 / (4 * pi)) * x * y * z,
            np.sqrt(35 / (32 * pi)) * x * (x * x - 3 * y * y),
            np.sqrt(35 / (32 * pi)) * y * (3 * x * x - y * y),
        ]
    ).T


class ConverterTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="openmx-test-")
        self.root = Path(self.directory.name)
        self.path = self.root / "C6.0.pao"
        fixture(self.path)

    def tearDown(self):
        self.directory.cleanup()

    def test_selection(self):
        s = BasisSpec.parse("C6.0-s2p2d1")
        self.assertEqual(s.nao, 13)
        self.assertEqual(s.counts, (2, 2, 1))
        self.assertEqual(sorted(s.permutation()), list(range(13)))
        self.assertEqual(BasisSpec.parse("C6.0-s1f1").counts, (1, 0, 0, 1))

    def test_bad_selection(self):
        for text in [
            "../C6.0-s1",
            "C6.0-s0",
            "C6.0-p1s1",
            "C6.0-s1s1",
            "C6.0-h1",
            "C6.0",
        ]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                BasisSpec.parse(text)

    def test_harmonic_signs_and_order(self):
        xyz = np.random.default_rng(714).normal(size=(40, 3))
        xyz /= np.linalg.norm(xyz, axis=1)[:, None]
        for l in range(5):
            cp = np.array([cp2k_harmonic(l, m, xyz) for m in range(-l, l + 1)]).T
            omx = openmx_harmonics(l, xyz)[:, ANGULAR_PERMUTATIONS[l]]
            np.testing.assert_allclose(cp, omx, atol=2e-15)

    def test_g_channel_file(self):
        fixture(self.path, lmax=4)
        pao = PaoFile.read(self.path)
        self.assertEqual(pao.orbitals[4].shape, (220, 2))
        self.assertEqual(BasisSpec.parse("C6.0-g1").nao, 9)

    def test_read_and_derivatives(self):
        pao = PaoFile.read(self.path)
        self.assertEqual(pao.element, "C")
        self.assertEqual(len(pao.sha256), 64)
        value, derivative = pao.radial(1, 2, np.array([1e-4, 0.5, 13.0]))
        self.assertTrue(np.isfinite(derivative).all())
        np.testing.assert_array_equal(value[-1], 0)
        np.testing.assert_array_equal(derivative[-1], 0)
        with self.assertRaises(ValueError):
            pao.radial(1, 3, np.array([1.0]))

    def test_fortran_exponents(self):
        self.path.write_text(self.path.read_text().replace("e-", "D-"))
        self.assertEqual(PaoFile.read(self.path).radius.shape, (220,))

    def test_reject_invalid_grid_and_blocks(self):
        original = self.path.read_text()
        for text in [
            original.replace("PAO.Mul 2", "PAO.Mul 3"),
            original.replace("pseudo.atomic.orbitals.L=2>", "broken>"),
            original.replace("AtomSpecies 6", "AtomSpecies 6\nAtomSpecies 8"),
            original.replace("valence.electron 4", "valence.electron NaN"),
        ]:
            self.path.write_text(text)
            with self.assertRaises(ValueError):
                PaoFile.read(self.path)

    def test_gaussian_norm_and_derivative(self):
        r, w = quadrature(np.linspace(0, 20, 400))
        alpha = np.array([0.05, 0.4, 2.0])
        for l in range(4):
            g, dg = primitives(r, alpha, l)
            np.testing.assert_allclose(
                g.T @ ((w * r * r)[:, None] * g), gaussian_overlap(alpha, l), atol=2e-13
            )
            hi, _ = primitives(r + 1e-6, alpha, l)
            lo, _ = primitives(r - 1e-6, alpha, l)
            np.testing.assert_allclose((hi - lo) / 2e-6, dg, atol=5e-10)

    def test_fit_and_cp2k_contraction(self):
        pao = PaoFile.read(self.path)
        fit = fit_channel(pao, 1, 2, nprimitive=24)
        self.assertLess(max(fit.diagnostics["normalized_radial_l2_errors"]), 2e-4)
        self.assertLess(fit.diagnostics["analytic_vs_quadrature_norm_max_error"], 1e-7)
        # CP2K first multiplies coefficients by alpha**((2*l+3)/4),
        # then normalizes each contraction. This must equal the exported radial fit.
        alpha, c, l = fit.exponents, fit.coefficients, fit.l
        raw = c * (2**l * (2 / np.pi) ** 0.75 * alpha ** ((2 * l + 3) / 4))[:, None]
        norm2 = np.einsum(
            "ij,ik,kj->j",
            raw,
            0.5 * gamma_local(l + 1.5) / (alpha[:, None] + alpha[None, :]) ** (l + 1.5),
            raw,
        )
        r = np.linspace(0.001, 12, 400)
        cp2k_radial = (r[:, None] ** l * np.exp(-r[:, None] ** 2 * alpha)) @ (
            raw / np.sqrt(norm2)
        )
        g, _ = primitives(r, alpha, l)
        np.testing.assert_allclose(cp2k_radial, g @ c, atol=5e-11)
        output = render_basis(BasisSpec.parse("C6.0-p2"), pao, [fit])
        self.assertIn("C OMX-FIT-C6.0-p2", output)
        self.assertIn("2 1 1 24 2", output)

    def test_cli_rejects_overwrite_and_source_collision(self):
        args = ["--pao-dir", str(self.root), "--output", str(self.path), "C6.0-s1"]
        before = self.path.read_bytes()
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main(args)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main(args + ["--overwrite"])
        self.assertEqual(before, self.path.read_bytes())

    def test_cli_rejects_element_prefix_mismatch(self):
        # H is a prefix of He, but they must not be accepted as the same element.
        source = self.root / "He6.0.pao"
        source.write_text(
            self.path.read_text().replace("AtomSpecies 6", "AtomSpecies 1")
        )
        args = [
            "--pao-dir",
            str(self.root),
            "--output",
            str(self.root / "basis"),
            "He6.0-s1",
        ]
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main(args)
        self.assertFalse((self.root / "basis").exists())


def gamma_local(x):
    return math.gamma(x)


if __name__ == "__main__":
    unittest.main()
