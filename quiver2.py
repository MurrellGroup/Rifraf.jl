import numpy as np
import scipy.sparse as sp

# FIXME: asymmetrical bandwidth is wrong thing to do. offset so it is symmetric.
# TODO: test BandedMatrix
# TODO: improve slicing and arithmetic operations of BandedMatrix
# TODO: port to Julia and compare speed

# TODO: use faster version of partial order aligner
# TODO: benchmark against actual quiver implementation


class OutsideBandError(Exception):
    pass


class OutOfBoundsError(Exception):
    pass


class BandedMatrix(object):
    """A sparse banded matrix.

    Currently only supports limited slicing and no arithmetic
    operations.

    """

    def __init__(self, shape, bandwidth, offset=0):
        nrows, ncols = shape
        if not nrows > 0 and ncols > 0:
            raise Exception('shape is not positive: {}'.format(shape))
        self.shape = shape
        self.bandwidth = bandwidth
        self.data = np.zeros((2 * bandwidth + 1, ncols))
        self.offset = offset

    @property
    def nrows(self):
        return self.shape[0]

    @property
    def ncols(self):
        return self.shape[1]

    def check_row(self, i):
        if i < 0:
            i = self.nrows + i
        if i >= self.nrows:
            raise OutOfBoundsError()
        return i

    def check_col(self, j, newcol=False):
        """if newcol is true, allow j=ncols, for inserting a new column"""
        if j < 0:
            j = self.ncols + j
        if j >= self.ncols + (1 if newcol else 0):
            raise OutOfBoundsError()
        return j

    def _convert_row(self, i, j, extend=False):
        """convert array row to data row"""
        o = 1 if extend else 0
        if i < max(j - self.bandwidth - self.offset, 0):
            raise OutsideBandError()
        if i > min(j + self.bandwidth - self.offset + o, self.nrows - (1 - o)):
            raise OutsideBandError()
        result = i - (j - self.bandwidth - self.offset)
        if result < 0 or result - o >= self.data.shape[0]:
            raise OutOfBoundsError()
        return result

    def range(self, j, newcol=False):
        """array-indexed row range of band in column j"""
        j = self.check_col(j, newcol)
        start = max(j - self.bandwidth - self.offset, 0)
        stop = min(j - self.offset + self.bandwidth + 1, self.nrows)
        assert start < stop
        if start == stop:
            raise OutOfBoundsError()
        return start, stop

    def data_range(self, j):
        """row range for column j in self.data"""
        j = self.check_col(j)
        start, stop = self.range(j)
        dstart = self._convert_row(start, j)
        dstop = self._convert_row(stop, j, extend=True)
        assert dstart < dstop
        return dstart, dstop

    def inband(self, i, j):
        start, stop = self.range(j)
        return i >= start and i < stop

    def get_elt(self, i, j):
        i = self.check_row(i)
        j = self.check_col(j)
        if self.inband(i, j):
            row = self._convert_row(i, j)
            return self.data[row, j]
        return 0

    def get_col(self, j):
        j = self.check_col(j)
        start, stop = self.data_range(j)
        return self.data[start:stop, j]

    def __getitem__(self, key):
        # TODO: only supports limited slicing
        # TODO: return BandedMatrix
        i, j = key
        if not isinstance(j, int):
            raise Exception()
        if isinstance(i, int):
            return self.get_elt(i, j)
        elif i == slice(None, None, None):
            return self.get_col(j)
        else:
            raise Exception()

    def __setitem__(self, key, val):
        # TODO: implement multidim slicing and negative indices
        i, j = key
        i = self.check_row(i)
        j = self.check_col(j)
        row = self._convert_row(i, j)
        self.data[row, j] = val

    def todense(self):
        result = np.zeros((self.nrows, self.ncols))
        for j in range(self.ncols):
            start, stop = self.range(j)
            dstart, dstop = self.data_range(j)
            result[start:stop, j] = self.data[dstart:dstop, j]
        return result


def update(arr, i, j, actual_j, s_base, t_base,
           log_p, log_ins, log_del, bandwidth):
    # TODO: log(1-p) may not be negligible
    sub = (0 if s_base == t_base else log_p)
    if i == actual_j - bandwidth:
        return max(arr[i, j - 1] + log_del,  # deletion
                   arr[i - 1, j - 1] + sub)
    elif i == actual_j + bandwidth:
        return max(arr[i - 1, j] + log_ins,  # insertion
                   arr[i - 1, j - 1] + sub)
    elif actual_j - bandwidth < i < actual_j + bandwidth:
        return max(max(arr[i - 1, j] + log_ins,  # insertion
                       arr[i, j - 1] + log_del),  # deletion
                   arr[i - 1, j - 1] + sub)
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


def updated_col(pos, base, template, seq, log_p,
                A, B, log_ins, log_del, bandwidth):
    Acols = BandedMatrix((A.nrows, 2), bandwidth, offset=-pos)
    Acols.data[:, 0] = A.data[:, pos]
    if pos < bandwidth:
        Acols[0, 1] = Acols[0, 0] + log_del
    j = pos + 1
    for i in range(max(1, j - bandwidth), min(A.nrows, j + bandwidth + 1)):
        Acols[i, 1] = update(Acols, i, 1, j, seq[i-1], base, log_p[i-1], log_ins, log_del, bandwidth)
    return Acols[:, 1]


def backward(s, log_p, t, log_ins, log_del, bandwidth):
    s = ''.join(reversed(s))
    log_p = np.array(list(reversed(log_p)))
    t = ''.join(reversed(t))
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


def score_mutation(mutation, template, seq, log_p, A, B, log_ins, log_del, bandwidth):
    """Score a mutation using the forward-backward trick."""
    mtype, pos, base = mutation
    if mtype == 'deletion':
        aj = pos
        Acol = A[:, aj]
    else:
        aj = pos + 1
        Acol = updated_col(pos, base, template, seq, log_p, A, B, log_ins, log_del, bandwidth)
    bj = pos + (0 if mtype == 'insertion' else 1)
    Bcol = B[:, bj]

    # need to chop off beginning of Acol and end of Bcol and align
    a_start, a_stop = A.range(aj, newcol=(mtype == "insertion" and aj == A.shape[1]))
    b_start, b_stop = B.range(bj)

    amin = max(b_start - a_start, 0)
    amax = len(Acol) - max(a_stop - b_stop, 0)

    bmin = max(a_start - b_start, 0)
    bmax = len(Bcol) - max(b_stop - a_stop, 0)

    # TODO: do this in cython
    return (Acol[amin:amax] + Bcol[bmin:bmax]).max()


def update_template(template, mutation):
    mtype, pos, base = mutation
    if mtype == 'substitution':
        return ''.join([template[:pos], base, template[pos + 1:]])
    elif mtype == 'insertion':
        return ''.join([template[:pos], base, template[pos:]])
    elif mtype == 'deletion':
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


def quiver2(template, sequences, phreds, log_ins, log_del, bandwidth=None,
            min_dist=9, max_iters=100, max_multi_iters=50, seed=None,
            verbose=False):
    """Generate an alignment-free consensus.

    sequences: list of dna sequences

    phreds: list of arrays of phred scores

    """
    log_ps = list((-phred / 10) for phred in phreds)

    _bandwidth = bandwidth
    lens = list(len(s) for s in sequences)
    maxlen = max(lens)
    minlen = min(lens)
    upper = max(maxlen, len(template))
    lower = min(minlen, len(template))
    min_bandwidth = max(1, 2 * (upper - lower))
    if _bandwidth is None:
        _bandwidth = min_bandwidth
        _bandwidth = max(min_bandwidth, _bandwidth)
    if _bandwidth < min_bandwidth:
        raise Exception('minimum bandwidth is {}, but given {} '.format(min_bandwidth, _bandwidth))

    if verbose:
        print("computing initial alignments")
    As = list(forward(s, p, template, log_ins, log_del, _bandwidth) for s, p in zip(sequences, log_ps))
    Bs = list(backward(s, p, template, log_ins, log_del, _bandwidth) for s, p in zip(sequences, log_ps))
    cur_score = sum(A[-1, -1] for A in As)
    for i in range(max_iters):
        if verbose:
            print("iteration {}".format(i))
            if bandwidth is None:
                print("  bandwidth: {}".format(_bandwidth))
        candidates = []
        for mutation in mutations(template):
            score = sum(score_mutation(mutation, template, seq, log_p, A, B, log_ins, log_del, _bandwidth)
                        for seq, log_p, A, B in zip(sequences, log_ps, As, Bs))
            if score > cur_score:
                candidates.append([score, mutation])
        if not candidates:
            break
        if verbose:
            print("  found {} candidate mutations".format(len(candidates)))
        chosen_cands = choose_candidates(candidates, min_dist, i, max_multi_iters)
        template = apply_mutations(template, chosen_cands)
        As = list(forward(s, p, template, log_ins, log_del, _bandwidth) for s, p in zip(sequences, log_ps))
        Bs = list(backward(s, p, template, log_ins, log_del, _bandwidth) for s, p in zip(sequences, log_ps))
        cur_score = sum(A[-1, -1] for A in As)

        if bandwidth is None:
            upper = max(maxlen, len(template))
            lower = min(minlen, len(template))
            _bandwidth = max(min_bandwidth, upper - lower)

        if verbose:
            print('  kept {} mutations'.format(len(chosen_cands)))
            print('  score: {}'.format(cur_score))
    return template
