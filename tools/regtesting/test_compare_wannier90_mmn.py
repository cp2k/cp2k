#!/usr/bin/env python3
"""Unit tests for the gauge-independent overlap reciprocity check."""

import unittest

import numpy as np

from compare_wannier90_mmn import check_reciprocity


class ReciprocityTests(unittest.TestCase):
    def setUp(self):
        self.forward = (1, 2, 1, 0, -1)
        self.reverse = (2, 1, -1, 0, 1)
        self.matrix = np.array([[0.9, 0.2j], [0.1 - 0.2j, 0.8j]])

    def data(self, reverse=None):
        if reverse is None:
            reverse = self.matrix.conj().T
        return (2, 2, 1, [(self.forward, self.matrix), (self.reverse, reverse)])

    def test_nonhermitian_block_and_reverse(self):
        self.assertEqual(check_reciprocity(self.data()), 0.0)

    def test_complex_unitary_gauge(self):
        left = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)
        right = np.diag([1j, -1])
        forward = left.conj().T @ self.matrix @ right
        reverse = right.conj().T @ self.matrix.conj().T @ left
        data = (2, 2, 1, [(self.forward, forward), (self.reverse, reverse)])
        self.assertLess(check_reciprocity(data), 1e-14)

    def test_wrong_sign(self):
        self.assertGreater(check_reciprocity(self.data(-self.matrix.conj().T)), 1.0)

    def test_missing_reverse(self):
        data = (2, 2, 1, [(self.forward, self.matrix), ((2, 1, 0, 0, 0), self.matrix)])
        with self.assertRaisesRegex(ValueError, "Missing reverse"):
            check_reciprocity(data)

    def test_duplicate_header(self):
        data = (2, 2, 1, [(self.forward, self.matrix), (self.forward, self.matrix)])
        with self.assertRaisesRegex(ValueError, "Repeated or missing"):
            check_reciprocity(data)

    def test_nonfinite_overlap(self):
        reverse = self.matrix.conj().T.copy()
        reverse[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "Invalid overlap"):
            check_reciprocity(self.data(reverse))


if __name__ == "__main__":
    unittest.main()
