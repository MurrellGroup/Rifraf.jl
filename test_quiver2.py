import random
import unittest

import numpy as np
from numpy.testing import assert_array_equal

from sample import phred
from sample import sample_from_template
from sample import sample
from quiver2 import update_template
from quiver2 import forward
from quiver2 import backward
from quiver2 import score_mutation
from quiver2 import updated_col
from quiver2 import quiver2


def random_seq(length):
    return ''.join(random.choice('ACGT') for _ in range(length))


def score_slow(template, sequence, log_p, log_ins, log_del, bandwidth):
    A = forward(sequence, log_p, template, log_ins, log_del, bandwidth)
    nrows, ncols = A.shape
    return A.get_elt(nrows - 1, ncols - 1)


class TestQuiver2(unittest.TestCase):

    # TODO: test asymmetric matrices
    def test_perfect_forward(self):
        log_del = -10
        log_ins = -5
        template = 'AA'
        seq = 'AA'
        log_p = np.repeat(-3, len(seq))
        A = forward(seq, log_p, template, log_ins, log_del, bandwidth=1).full()
        expected = np.array([[  0, -10, 0],
                             [ -5,   0, -10],
                             [0,  -5,   0]])
        assert_array_equal(A, expected)

    def test_backward(self):
        log_del = -10
        log_ins = -5
        template = 'AA'
        seq = 'AT'
        log_p = np.repeat(-3, len(seq))
        B = backward(seq, log_p, template, log_ins, log_del, bandwidth=1).full()
        expected = np.array([[-3, -5, 0],
                             [-13, -3, -5],
                             [0, -10, 0]])
        assert_array_equal(B, expected)

    def test_imperfect_forward(self):
        log_del = -10
        log_ins = -5
        template = 'AA'
        seq = 'AT'
        log_p = np.repeat(-3, len(seq))
        A = forward(seq, log_p, template, log_ins, log_del, bandwidth=1).full()
        expected = np.array([[  0, -10, 0],
                             [ -5,  0,  -10],
                             [0, -5,  -3]])
        assert_array_equal(A, expected)

    def test_updated_col_substitution(self):
        log_del = -10
        log_ins = -5
        template = 'TAAG'
        seq = 'ATAG'
        pos = 1
        log_p = np.repeat(-3, len(seq))
        A = forward(seq, log_p, template, log_ins, log_del, bandwidth=1)
        result = updated_col(pos, pos + 1, 'T', template, seq, log_p, A, log_ins, log_del)
        expected = forward(seq, log_p, 'TTAG', log_ins, log_del, bandwidth=1).get_col(pos + 1)
        assert_array_equal(result, expected)

    def test_random_mutations(self):
        point_rate = 0.1
        insertion_rate = 0.01
        deletion_rate = 0.01
        log_ins = np.log10(insertion_rate)
        log_del = np.log10(deletion_rate)
        for _ in range(1000):
            template_len = random.randint(5, 20)
            template_seq = random_seq(template_len)
            template = template_seq
            seq = sample_from_template(template_seq, point_rate, insertion_rate, deletion_rate)
            bandwidth = max(2 * np.abs(len(template) - len(seq)), 5)
            f = random.choice(['substitution', 'insertion', 'deletion'])
            maxpos = len(template) if f == 'insertion' else len(template) - 1
            pos = random.randint(0, maxpos)
            base = random.choice('ACGT')
            mutation = (f,
                        pos,
                        base)
            new_template = update_template(template, mutation)
            phreds = phred(np.repeat(point_rate + insertion_rate + deletion_rate, len(seq)))
            log_p = -phreds / 10
            A = forward(seq, log_p, template, log_ins, log_del, bandwidth)
            B = backward(seq, log_p, template, log_ins, log_del, bandwidth)
            M = forward(seq, log_p, new_template, log_ins, log_del, bandwidth)
            nr, nc = M.shape
            if f == 'substitution':
                # this is the only case where the updated column exactly matches the full result
                col = updated_col(pos, pos + 1, base, template, seq, log_p, A, log_ins, log_del)
                exp_col = M.get_col(pos + 1)
                assert_array_equal(col, exp_col)
            score = score_mutation(mutation, template, seq, log_p, A, B, log_ins, log_del, bandwidth)
            score2 = score_slow(new_template, seq, log_p, log_ins, log_del, bandwidth)
            score3 = M.get_elt(nr - 1, nc - 1)
            score4 = backward(seq, log_p, new_template, log_ins, log_del, bandwidth).get_elt(0, 0)
            self.assertAlmostEqual(score, score2)
            self.assertAlmostEqual(score, score3)
            self.assertAlmostEqual(score, score4)

    def test_quiver2(self):
        for _ in range(10):
            template, reads, phreds = sample(n=10, length=20, error_rate=3/100, insertion_rate=3/100, deletion_rate=3/100)
            result = quiver2(reads[0], reads, phreds, np.log10(3/100), np.log10(3/100), verbose=False, bandwidth=3, min_dist=9)
            self.assertEqual(template, result)


if __name__ == '__main__':
    unittest.main()
