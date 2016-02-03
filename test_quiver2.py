import unittest

import numpy as np
from numpy.testing import assert_array_equal

from sample import phred
from quiver2 import seq_to_array
from quiver2 import forward
from quiver2 import backward
from quiver2 import substitution
from quiver2 import insertion
from quiver2 import deletion
from quiver2 import score_mutation


class TestQuiver2(unittest.TestCase):

    def test_perfect_forward(self):
        log_del = -10
        log_ins = -5
        template = seq_to_array('AA')
        seq_array = seq_to_array('AA')
        phreds = np.repeat(30, len(seq_array))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        expected = np.array([[  0, -10, -20],
                             [ -5,   0, -10],
                             [-10,  -5,   0]])
        assert_array_equal(A, expected)

    def test_imperfect_forward(self):
        log_del = -10
        log_ins = -5
        template = seq_to_array('AA')
        seq_array = seq_to_array('AT')
        phreds = np.repeat(3, len(seq_array))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        expected = np.array([[  0, -10, -20],
                             [ -5,  0,  -10],
                             [-10, -5,  -3]])
        assert_array_equal(A, expected)

    def test_substitution(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.01)
        mutation = ('substitution', 1, 0)
        template = seq_to_array('ATAA')
        seq_array = seq_to_array('AAAA')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        Acols, Bcol = substitution(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        expected = forward(seq_array, phreds, seq_array, log_ins, log_del)[:, [2, 3]]
        assert_array_equal(Acols, expected)
        assert_array_equal(Bcol, B[:, 3])

    def test_insertion(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.01)
        mutation = ('insertion', 1, 3)
        template = seq_to_array('AAA')
        seq_array = seq_to_array('ATAA')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        Acols, Bcol = insertion(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        expected = forward(seq_array, phreds, seq_array, log_ins, log_del)[:, [2, 3]]
        assert_array_equal(Acols, expected)
        assert_array_equal(Bcol, B[:, 2])

    def test_deletion(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.01)
        mutation = ('deletion', 1, None)
        template = seq_to_array('ATAA')
        seq_array = seq_to_array('AAA')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        Acols, Bcol = deletion(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        expected = forward(seq_array, phreds, seq_array, log_ins, log_del)[:, [2, 3]]
        assert_array_equal(Acols, expected)
        assert_array_equal(Bcol, B[:, 3])

    def test_score_mutation(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.01)
        mutation = ('substitution', 1, 0)
        template = seq_to_array('ATAA')
        seq_array = seq_to_array('AAAA')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        result = score_mutation(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        expected = forward(seq_array, phreds, seq_array, log_ins, log_del)[-1, -1]
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
