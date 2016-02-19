import random
import unittest

import numpy as np
from numpy.testing import assert_array_equal

from quiver2 import BandedMatrix


class TestBlockMatrix(unittest.TestCase):

    def test_inband(self):
        m = BandedMatrix((13, 11), 5)
        self.assertTrue(m.inband(0, 0))
        self.assertTrue(m.inband(7, 0))
        self.assertTrue(not m.inband(8, 0))
        self.assertTrue(m.inband(0, 5))
        self.assertTrue(not m.inband(0, 6))

    def test_sym(self):
        m = BandedMatrix((3, 3), 0)
        m.data[:] = 1
        assert_array_equal(m.full(), np.diag(np.ones(3)))

    def test_sym_band(self):
        m = BandedMatrix((3, 3), 1)
        m.data[:] = 1
        expected = np.ones((3, 3))
        expected[-1, 0] = 0
        expected[0, -1] = 0
        assert_array_equal(m.full(), expected)

    def test_wide(self):
        m = BandedMatrix((3, 4), 0)
        m.data[:] = 1
        a = np.ones(4)
        b = np.ones(3)
        expected = (np.diag(a, 0) + np.diag(b, 1))[:3, :]
        assert_array_equal(m.full(), expected)

    def test_wide_band(self):
        m = BandedMatrix((3, 5), 1)
        m.data[:] = 1
        expected = np.ones((3, 5))
        expected[-1, 0] = expected[0, -1] = 0
        assert_array_equal(m.full(), expected)

    def test_wide_col(self):
        m = BandedMatrix((3, 5), 1)
        m.data[:] = 1
        first = np.ones(2)
        middle = np.ones(3)
        last = first
        assert_array_equal(m.get_col(0), first)
        assert_array_equal(m.get_col(1), middle)
        assert_array_equal(m.get_col(2), middle)
        assert_array_equal(m.get_col(3), middle)
        assert_array_equal(m.get_col(4), last)

    def test_tall(self):
        m = BandedMatrix((4, 3), 0)
        m.data[:] = 1
        a = np.ones(4)
        b = np.ones(3)
        expected = (np.diag(a, 0) + np.diag(b, -1))[:, :3]
        assert_array_equal(m.full(), expected)

    def test_tall_band(self):
        m = BandedMatrix((5, 3), 1)
        m.data[:] = 1
        expected = np.ones((5, 3))
        expected[-1, 0] = expected[0, -1] = 0
        assert_array_equal(m.full(), expected)


if __name__ == '__main__':
    unittest.main()
