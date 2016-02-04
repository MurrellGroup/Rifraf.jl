import random
import unittest

import numpy as np
from numpy.testing import assert_array_equal

from sample import phred
from sample import sample_from_template
from quiver2 import seq_to_array
from quiver2 import update_template
from quiver2 import forward
from quiver2 import backward
from quiver2 import substitution
from quiver2 import insertion
from quiver2 import deletion
from quiver2 import score_mutation
from quiver2 import score_slow


def random_seq(length):
    return ''.join(random.choice('ACGT') for _ in range(length))


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

    def test_backward(self):
        log_del = -10
        log_ins = -5
        template = seq_to_array('AA')
        seq_array = seq_to_array('AT')
        phreds = np.repeat(30, len(seq_array))
        B = backward(seq_array, phreds, template, log_ins, log_del)
        expected = np.array([[-3, -5, -10],
                             [-13, -3, -5],
                             [-20, -10, 0]])
        assert_array_equal(B, expected)

    def test_imperfect_forward(self):
        log_del = -10
        log_ins = -5
        template = seq_to_array('AA')
        seq_array = seq_to_array('AT')
        phreds = np.repeat(30, len(seq_array))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        expected = np.array([[  0, -10, -20],
                             [ -5,  0,  -10],
                             [-10, -5,  -3]])
        assert_array_equal(A, expected)

    def test_substitution(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.01)
        mutation = (substitution, 1, 0)
        template = seq_to_array('ATTA')
        new_template = seq_to_array('AATA')
        assert_array_equal(update_template(template, mutation), new_template)
        seq_array = seq_to_array('AAAA')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        Acol, Bcol = substitution(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        Aup = forward(seq_array, phreds, new_template, log_ins, log_del)
        score = score_mutation(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        # score2 = sum(scores_slow(new_template, [seq_array], [phreds], log_ins, log_del))
        # self.assertEqual(score, score2)
        self.assertEqual(score, Aup[-1, -1])

    def test_substitution_last(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.01)
        mutation = (substitution, 2, 0)
        template = seq_to_array('ATT')
        new_template = seq_to_array('ATA')
        assert_array_equal(update_template(template, mutation), new_template)
        seq_array = seq_to_array('AAA')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        Acol, Bcol = substitution(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        Aup = forward(seq_array, phreds, new_template, log_ins, log_del)[:, [2, 3]]
        score = score_mutation(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        # score2 = sum(scores_slow(new_template, [seq_array], [phreds], log_ins, log_del))
        # self.assertEqual(score, score2)
        self.assertEqual(score, Aup[-1, -1])

    def test_insertion(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.01)
        mutation = (insertion, 1, 3)
        template = seq_to_array('AAA')
        new_template = seq_to_array('ATAA')
        assert_array_equal(update_template(template, mutation), new_template)
        seq_array = seq_to_array('ATTAA')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        Acol, Bcol = insertion(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        Aup = forward(seq_array, phreds, new_template, log_ins, log_del)
        assert_array_equal(Bcol, B[:, 1])
        score = score_mutation(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        # score2 = sum(scores_slow(new_template, [seq_array], [phreds], log_ins, log_del))
        # self.assertEqual(score, score2)
        self.assertEqual(score, Aup[-1, -1])

    def test_insertion_last(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.01)
        mutation = (insertion, 3, 3)
        template = seq_to_array('AAA')
        new_template = seq_to_array('AAAT')
        assert_array_equal(update_template(template, mutation), new_template)
        seq_array = seq_to_array('AAATT')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        Acol, Bcol = insertion(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        Aup = forward(seq_array, phreds, new_template, log_ins, log_del)
        Bup = backward(seq_array, phreds, new_template, log_ins, log_del)
        score = score_mutation(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        # score2 = sum(scores_slow(new_template, [seq_array], [phreds], log_ins, log_del))
        # self.assertEqual(score, score2)
        self.assertEqual(score, Aup[-1, -1])

    def test_deletion(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.01)
        mutation = (deletion, 1, None)
        template = seq_to_array('ATTAA')
        new_template = seq_to_array('ATAA')
        assert_array_equal(update_template(template, mutation), new_template)
        seq_array = seq_to_array('AAA')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        Acol, Bcol = deletion(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        Aup = forward(seq_array, phreds, new_template, log_ins, log_del)
        assert_array_equal(Bcol, B[:, 2])
        score = score_mutation(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        # score2 = sum(scores_slow(new_template, [seq_array], [phreds], log_ins, log_del))
        # self.assertEqual(score, score2)
        self.assertEqual(score, Aup[-1, -1])

    def test_deletion_last(self):
        log_ins = np.log10(0.001)
        log_del = np.log10(0.01)
        mutation = (deletion, 3, None)
        template = seq_to_array('AAAT')
        new_template = seq_to_array('AAA')
        assert_array_equal(update_template(template, mutation), new_template)
        seq_array = seq_to_array('AAT')
        phreds = phred(np.repeat(0.1, len(seq_array)))
        A = forward(seq_array, phreds, template, log_ins, log_del)
        B = backward(seq_array, phreds, template, log_ins, log_del)
        Acol, Bcol = deletion(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        Aup = forward(seq_array, phreds, new_template, log_ins, log_del)
        score = score_mutation(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
        # score2 = sum(scores_slow(new_template, [seq_array], [phreds], log_ins, log_del))
        # self.assertEqual(score, score2)
        self.assertEqual(score, Aup[-1, -1])

    def test_random_mutations(self):
        point_rate = 0.1
        insertion_rate = 0.01
        deletion_rate = 0.01
        log_ins = np.log10(insertion_rate)
        log_del = np.log10(deletion_rate)
        for _ in range(10000):
            template_len = random.randint(3, 20)
            template_seq = random_seq(template_len)
            template = seq_to_array(template_seq)
            seq = sample_from_template(template_seq, point_rate, insertion_rate, deletion_rate)
            seq_array = seq_to_array(seq)
            f = random.choice([substitution, insertion, deletion])
            mutation = (f,
                        random.randint(0, len(template) - 2),
                        random.choice([0, 1, 2, 3]))
            new_template = update_template(template, mutation)
            phreds = phred(np.repeat(point_rate + insertion_rate + deletion_rate, len(seq_array)))
            A = forward(seq_array, phreds, template, log_ins, log_del)
            B = backward(seq_array, phreds, template, log_ins, log_del)
            score = score_mutation(mutation, template, seq_array, phreds, A, B, log_ins, log_del)
            score2 = score_slow(new_template, seq_array, phreds, log_ins, log_del)
            score3 = forward(seq_array, phreds, new_template, log_ins, log_del)[-1, -1]
            self.assertAlmostEqual(score, score2)
            self.assertAlmostEqual(score, score3)


if __name__ == '__main__':
    unittest.main()
