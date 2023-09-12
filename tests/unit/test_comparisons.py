from unittest import TestCase

import numpy as np

from python_pdb.comparisons import rmsd


class TestRmsd(TestCase):
    def test(self):
        arr1 = np.array([[1.00, 0.00, 0.00],
                         [3.00, 0.00, 1.00],
                         [2.00, 0.00, 1.00]])
        arr2 = np.array([[6.00, 0.00, 0.00],
                         [1.00, 0.70, 1.00],
                         [4.00, 0.00, 1.00]])

        self.assertAlmostEqual(rmsd(arr1, arr2), 3.3411575)

    def test_zeros(self):
        arr1 = np.zeros((10, 3))
        arr2 = np.zeros((10, 3))

        self.assertAlmostEqual(rmsd(arr1, arr2), 0.00)
