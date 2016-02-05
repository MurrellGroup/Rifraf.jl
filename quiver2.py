import random
import numpy as np
import scipy.sparse as sp

# TODO: try column-major order
# TODO: rewrite BandedMatrix as Cython extension type
# TODO: either write slow code in Cython, or rewrite in Julia
# TODO: use pbdagcon for initial template
# TODO: do not completely recompute As and Bs each iteration. Only recompute from (first) Acol on.
# TODO: compress duplicate reads and account for their scores


class OutsideBandException(Exception):
    pass


class BandedMatrix(object):
    """A sparse banded matrix.

    Currently only supports limited slicing and no arithmetic
    operations.

    """

    def __init__(self, shape, bandwidth, offset=0):
        nrows, ncols = shape
        if not nrows > 0 and ncols > 0:
            raise Exception()
        self.shape = shape
        self.bandwidth = bandwidth
        self.data = np.zeros((2 * bandwidth + 1, ncols))
        self.offset = offset

    def _convert_row(self, i, j, open=False):
        x = 1 if open else 0
        y = 0 if open else 1
        if i < max(j - self.bandwidth - self.offset, 0):
            raise OutsideBandException()
        if i > min(j + self.bandwidth - self.offset + x, self.shape[0] - y):
            raise OutsideBandException()
        return i - (j - self.bandwidth - self.offset)

    def range(self, j):
        start = max(j - self.bandwidth - self.offset, 0)
        stop = min(j - self.offset + self.bandwidth + 1, self.shape[0])
        return start, stop

    def data_range(self, j):
        start, stop = self.range(j)
        dstart = self._convert_row(start, j)
        dstop = self._convert_row(stop, j, open=True)
        return dstart, dstop

    def __getitem__(self, key):
        # TODO: only supports limited slicing
        i, j = key
        if not isinstance(j, int):
            raise Exception()
        if j < 0:
            j = self.shape[1] + j
        if isinstance(i, int):
            if i < 0:
                i = self.shape[0] + i
            try:
                row = self._convert_row(i, j)
                return self.data[row, j]
            except OutsideBandException:
                return 0
        elif i == slice(None, None, None):
            start, stop = self.data_range(j)
            # TODO: return BandedMatrix
            return self.data[start:stop, j]
        else:
            raise Exception()

    def __setitem__(self, key, val):
        # TODO: implement multidim slicing and negative indices
        i, j = key
        row = self._convert_row(i, j)
        self.data[row, j] = val

    def todense(self):
        result = np.zeros(self.shape)
        for j in range(self.shape[1]):
            start, stop = self.range(j)
            dstart, dstop = self.data_range(j)
            result[start:stop, j] = self.data[dstart:dstop, j]
        return result


def update(arr, i, j, actual_j, s_base, t_base, log_p, log_ins, log_del, bandwidth):
    # TODO: log(1-p) may not be negligible
    if i == actual_j - bandwidth:
        return max([arr[i, j - 1] + log_del,  # deletion
                    arr[i - 1, j - 1] + (0 if s_base == t_base else log_p)])
    elif i == actual_j + bandwidth:
        return max([arr[i - 1, j] + log_ins,  # insertion
                    arr[i - 1, j - 1] + (0 if s_base == t_base else log_p)])
    elif actual_j - bandwidth < i < actual_j + bandwidth:
        return max([arr[i - 1, j] + log_ins,  # insertion
                    arr[i, j - 1] + log_del,  # deletion
                    arr[i - 1, j - 1] + (0 if s_base == t_base else log_p)])
    else:
        raise Exception('update called outside bandwidth')


def forward(s, log_p, t, log_ins, log_del, bandwidth):
    result = BandedMatrix((len(s) + 1, len(t) + 1), bandwidth)
    for i in range(min(bandwidth + 1, len(s) + 1)):
        result[i, 0] = log_ins * i
    for j in range(min(bandwidth + 1, len(t) + 1)):
        result[0, j] = log_del * j
    for i in range(1, len(s) + 1):
        for j in range(max(1, i - bandwidth), min(len(t) + 1, i + bandwidth + 1)):
            result[i, j] = update(result, i, j, j, s[i-1], t[j-1], log_p[i-1], log_ins, log_del, bandwidth)
    return result


def backward(s, log_p, t, log_ins, log_del, bandwidth):
    s = list(reversed(s))
    log_p = list(reversed(log_p))
    t = list(reversed(t))
    result = forward(s, log_p, t, log_ins, log_del, bandwidth)
    result.data = np.flipud(np.fliplr(result.data))
    result.offset = len(t) - len(s)
    return result


def mutations(template):
    """Returns (name, position, base)"""
    for j in range(len(template)):
        for base in 'ACGT':
            if template[j] == base:
                continue
            yield ['substitution', j, base]
        yield ['deletion', j, None]
        for base in 'ACGT':
            yield ['insertion', j, base]
    # insertion after last
    for base in 'ACGT':
        yield ['insertion', len(template), base]


def updated_col(pos, base, template, seq, log_p, A, B, log_ins, log_del, bandwidth):
    Acols = BandedMatrix((A.shape[0], 2), bandwidth, offset=-pos)
    Acols.data[:, 0] = A.data[:, pos]
    if pos < bandwidth:
        Acols[0, 1] = Acols[0, 0] + log_del
    j = pos + 1
    for i in range(max(1, j - bandwidth), min(A.shape[0], j + bandwidth + 1)):
        Acols[i, 1] = update(Acols, i, 1, j, seq[i-1], base, log_p[i-1], log_ins, log_del, bandwidth)
    return Acols[:, 1]


b_offset = {'substitution': 1,
            'insertion': 0,
            'deletion': 1}


def score_mutation(mutation, template, seq, log_p, A, B, log_ins, log_del, bandwidth):
    """Score a mutation using the forward-backward trick."""
    f, pos, base = mutation
    if f == 'deletion':
        aj = pos
        Acol = A[:, aj]
    else:
        aj = pos + 1
        Acol = updated_col(pos, base, template, seq, log_p, A, B, log_ins, log_del, bandwidth)
    bj = pos + b_offset[f]
    Bcol = B[:, bj]

    # need to chop off beginning of Acol and end of Bcol and align
    a_start, a_stop = A.range(aj)
    b_start, b_stop = B.range(bj)

    amin = max(b_start - a_start, 0)
    amax = len(Acol) - max(a_stop - b_stop, 0)

    bmin = max(a_start - b_start, 0)
    bmax = len(Bcol) - max(b_stop - a_stop, 0)

    assert amax - amin == bmax - bmin
    assert Acol[amin:amax].shape == Bcol[bmin:bmax].shape

    return (Acol[amin:amax] + Bcol[bmin:bmax]).max()


def update_template(template, mutation):
    f, pos, base = mutation
    if f == 'substitution':
        return ''.join([template[:pos], base, template[pos + 1:]])
    elif f == 'insertion':
        return ''.join([template[:pos], base, template[pos:]])
    elif f == 'deletion':
        return ''.join([template[:pos], template[pos + 1:]])
    else:
        raise Exception('unknown mutation: {}'.format(mtype))


def apply_mutations(template, mutations):
    for i in range(len(mutations)):
        score, mutation = mutations[i]
        template = update_template(template, mutation)
        mtype, pos, _ = mutation
        if mtype == 'insertion':
            for ii in range(i + 1, len(mutations)):
                if mutations[ii][1][1] > pos:
                    mutations[ii][1][1] += 1
        elif mtype == 'deletion':
            for ii in range(i + 1, len(mutations)):
                if mutations[ii][1][1] > pos:
                    mutations[ii][1][1] -= 1
    return template


def choose_candidates(candidates, min_dist, i, max_multi_iters):
    final_cands = []
    if i < max_multi_iters:
        posns = set()
        for c in reversed(sorted(candidates, key=lambda c: c[0])):
            _, (_, posn, _) = c
            if any(abs(posn - p) < min_dist for p in posns):
                continue
            posns.add(posn)
            final_cands.append(c)
    else:
        final_cands = [list(sorted(candidates, key=lambda c: c[0]))[-1]]
    return final_cands


def quiver2(sequences, phreds, log_ins, log_del, bandwidth=10,
            min_dist=9, max_iters=100, max_multi_iters=50, seed=None,
            verbose=False):
    """
    sequences: list of dna sequences

    phreds: list of numpy array

    """
    log_ps = list(-phred / 10 for phred in phreds)
    state = np.random.RandomState(seed)
    template = state.choice(sequences)

    lens = list(len(s) for s in sequences)
    min_bandwidth = 2 * (max(lens) - min(lens))
    bandwidth = max(min_bandwidth, bandwidth)
    if bandwidth < 1:
        raise Exception('{} bandwidth is too small'.format(bandwidth))

    As = list(forward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
    Bs = list(backward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
    cur_score = sum(A[-1, -1] for A in As)
    for i in range(max_iters):
        if verbose:
            print("iteration {}".format(i))
        candidates = []
        for mutation in mutations(template):
            score = sum(score_mutation(mutation, template, seq, log_p, A, B, log_ins, log_del, bandwidth)
                        for seq, log_p, A, B in zip(sequences, log_ps, As, Bs))
            if score > cur_score:
                candidates.append([score, mutation])
        if not candidates:
            break
        chosen_cands = choose_candidates(candidates, min_dist, i, max_multi_iters)
        template = apply_mutations(template, chosen_cands)
        As = list(forward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
        Bs = list(backward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
        cur_score = sum(A[-1, -1] for A in As)
        if verbose:
            print('  kept {} mutations'.format(len(chosen_cands)))
            print('  score: {}'.format(cur_score))
    return template
