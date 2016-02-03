import unittest

import numpy as np

from sample import phred
from quiver2 import seq_to_array
from quiver2 import forward
from quiver2 import backward
from quiver2 import substitution
from quiver2 import score_mutation


class TestQuiver2(unittest.TestCase):

    def test_perfect_forward(self):
        log_ins = -10
        log_del = -10
        template = seq_to_array('AA')
        seq_array = seq_to_array('AA')
        phreds = np.repeat(30, len(seq_array))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        expected = np.array([[0,   -10, -20],
                             [-10,   0, -10],
                             [-20, -10,   0]])
        self.assertTrue(np.all(A == expected))

    def test_imperfect_forward(self):
        log_ins = -10
        log_del = -10
        template = seq_to_array('AA')
        seq_array = seq_to_array('AT')
        phreds = np.repeat(5, len(seq_array))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        expected = np.array([[0,   -10, -20],
                             [-10,   0, -10],
                             [-20, -10,   -5]])
        self.assertTrue(np.all(A == expected))


    def test_substitution(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.001)
        mutation = ('substitution', 1, 0)
        template = seq_to_array('ATAA')
        seq_array = seq_to_array('AAAA')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        Acols, Bcol = substitution(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        expected = forward(seq_array, phreds, seq_array, log_ins, log_del)[:, 1:4]
        self.assertTrue(np.all(Acols == expected))
        self.assertTrue(np.all(Bcol == B[3]))


    def test_score_mutation(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.001)
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
