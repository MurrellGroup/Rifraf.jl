import numpy as np

# TODO: cache BandedMatrix calls
# TODO: test BandedMatrix
# TODO: improve slicing and arithmetic operations of BandedMatrix
# TODO: port to Julia and compare speed

# TODO: use faster version of partial order aligner
# TODO: benchmark against actual quiver implementation


class OutsideBandError(Exception):
    pass


class OutOfBoundsError(Exception):
    pass


class BandedMatrix:
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
        if not j < self.ncols + (1 if newcol else 0):
            raise OutOfBoundsError()
        return j

    def square_row_range(self, j, newcol=False):
        """array-indexed theoretical row range of band in column j, assuming
        square matrix and band extending above first row and below
        last.

        """
        j = self.check_col(j, newcol)
        start = j - self.bandwidth - self.offset
        stop = j - self.offset + self.bandwidth + 1
        return start, stop

    def row_range(self, j, newcol=False):
        """array-indexed real row range of band in column j"""
        # FIXME: assumes every column contains nonzero entries
        start, stop = self.square_row_range(j, newcol)
        start = max(0, start)
        stop = min(stop, self.nrows)
        assert start < stop
        return start, stop

    def inband(self, i, j):
        start, stop = self.row_range(j)
        return start <= i < stop

    def _convert_row(self, i, j):
        """convert array row to data row"""
        if not self.inband(i, j):
            raise OutsideBandError()
        start, stop = self.row_range(j)
        result = i - (j - self.bandwidth - self.offset)
        assert 0 <= result < self.data.shape[0]
        return result

    def data_row_range(self, j):
        """row range for column j in self.data"""
        j = self.check_col(j)
        start, stop = self.row_range(j)
        dstart = self._convert_row(start, j)
        dstop = self._convert_row(stop - 1, j) + 1
        assert dstart < dstop
        return dstart, dstop

    def get_elt(self, i, j):
        i = self.check_row(i)
        j = self.check_col(j)
        if self.inband(i, j):
            row = self._convert_row(i, j)
            return self.data[row, j]
        return 0

    def get_col(self, j):
        j = self.check_col(j)
        start, stop = self.data_row_range(j)
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
            raise Exception('given slicing not recognized: {}'.format(key))

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
            start, stop = self.row_range(j)
            dstart, dstop = self.data_row_range(j)
            result[start:stop, j] = self.data[dstart:dstop, j]
        return result

    def col_indices(self, i):
        """Generate all band indices"""
        i = self.check_row(i)
        # FIXME: this is inefficient
        for j in range(self.ncols):
            if self.inband(i, j):
                yield j

    def row_indices(self, j):
        """Generate all band indices"""
        j = self.check_col(j)
        start, stop = self.row_range(j)
        return range(start, stop)

    def indices(self):
        """Generate all band indices in [1:, 1:]"""
        for j in range(1, self.ncols):
            for i in range(*self.row_range(j)):
                if i > 0:
                    yield i, j

    def predecessors(self, i, j):
        """Generate all left, up, and left-up band indices"""
        mods = ((-1, 0), (0, -1), (-1, -1))
        for imod, jmod in mods:
            ii = i + imod
            jj = j + jmod
            if ii < 0 or jj < 0:
                continue
            if self.inband(ii, jj):
                yield ii, jj

    def flip(self):
        # FIXME: not symmetric
        self.data = np.flipud(np.fliplr(self.data))
        self.offset += divmod(self.nrows - self.ncols, 2)[1]


class CenteredBandedMatrix(BandedMatrix):
    """A banded matrix guaranteed to have the band touch elements [0, 0] and [-1, -1], plus some optional extra."""
    def __init__(self, shape, bandwidth):
        self.given_bandwidth = bandwidth
        nrows, ncols = shape
        diff = ncols - nrows
        offset, r = divmod(diff, 2)
        self.extra_bandwidth = max(bandwidth, r)
        bandwidth = self.extra_bandwidth + np.abs(offset)
        super().__init__(shape, bandwidth, offset=offset)


def update(arr, i, j, s_base, t_base, log_p, log_ins, log_del):
    # TODO: log(1-p) may not be negligible
    sub = (0 if s_base == t_base else log_p)
    scores = []
    if arr.inband(i - 1, j):
        # insertion
        scores.append(arr[i - 1, j] + log_ins)
    if arr.inband(i, j - 1):
        # deletion
        scores.append(arr[i, j - 1] + log_del)
    if arr.inband(i - 1, j - 1):
        scores.append(arr[i - 1, j - 1] + sub)
    return max(scores)


def forward(s, log_p, t, log_ins, log_del, bandwidth):
    result = CenteredBandedMatrix((len(s) + 1, len(t) + 1), bandwidth)
    for i in result.row_indices(0):
        result[i, 0] = log_ins * i
    for j in result.col_indices(0):
        result[0, j] = log_del * j
    for i, j in result.indices():
        result[i, j] = update(result, i, j, s[i-1], t[j-1], log_p[i-1], log_ins, log_del)
    return result


def updated_col(pos, base, template, seq, log_p,
                A, B, log_ins, log_del, bandwidth):
    Acols = BandedMatrix((A.nrows, 2), A.bandwidth, offset=A.offset - pos)
    Acols.data[:, 0] = A.data[:, pos]
    if Acols.inband(0, 0) and Acols.inband(0, 1):
        Acols[0, 1] = Acols[0, 0] + log_del
    for i in Acols.row_indices(1):
        Acols[i, 1] = update(Acols, i, 1, seq[i-1], base, log_p[i-1], log_ins, log_del)
    return Acols[:, 1]


def backward(s, log_p, t, log_ins, log_del, bandwidth):
    s = ''.join(reversed(s))
    log_p = np.array(list(reversed(log_p)))
    t = ''.join(reversed(t))
    result = forward(s, log_p, t, log_ins, log_del, bandwidth)
    result.flip()
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
    a_start, a_stop = A.row_range(aj, newcol=(mtype == "insertion" and aj == A.shape[1]))
    b_start, b_stop = B.row_range(bj)

    amin = max(b_start - a_start, 0)
    amax = len(Acol) - max(a_stop - b_stop, 0)

    bmin = max(a_start - b_start, 0)
    bmax = len(Bcol) - max(b_stop - a_stop, 0)

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


def choose_candidates(candidates, min_dist)
    final_cands = []
    posns = set()
    for c in reversed(sorted(candidates, key=lambda c: c[0])):
        _, (_, posn, _) = c
        if any(abs(posn - p) < min_dist for p in posns):
            continue
        posns.add(posn)
        final_cands.append(c)
    return final_cands


def quiver2(template, sequences, phreds, log_ins, log_del, bandwidth=10,
            min_dist=9, max_iters=100, verbose=False):
    """Generate an alignment-free consensus.

    sequences: list of dna sequences

    phreds: list of arrays of phred scores

    """
    log_ps = list((-phred / 10) for phred in phreds)

    if bandwidth < 0:
        raise Exception('bandwidth cannot be negative: {} '.format(bandwidth))

    if verbose:
        print("computing initial alignments")
    As = list(forward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
    Bs = list(backward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
    best_score = sum(A[-1, -1] for A in As)
    for i in range(max_iters):
        old_score = best_score
        if verbose:
            print("iteration {}".format(i))
        candidates = []
        for mutation in mutations(template):
            score = sum(score_mutation(mutation, template, seq, log_p, A, B, log_ins, log_del, bandwidth)
                        for seq, log_p, A, B in zip(sequences, log_ps, As, Bs))
            if score > best_score:
                candidates.append([score, mutation])
        if not candidates:
            break
        if verbose:
            print("  found {} candidate mutations".format(len(candidates)))
        chosen_cands = choose_candidates(candidates, min_dist)
        template = apply_mutations(template, chosen_cands)
        As = list(forward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
        Bs = list(backward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
        best_score = sum(A[-1, -1] for A in As)

        # detect if a single mutation is better
        # FIXME: code duplication
        if best_score < chosen_cands[0][0]:
            if verbose:
                print('  rejecting multiple candidates in favor of best')
            chosen_cands = [chosen_cands[0]]
            template = apply_mutations(template, chosen_cands)
            As = list(forward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
            Bs = list(backward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
            best_score = sum(A[-1, -1] for A in As)
            assert best_score == chosen_cands[0][0]
        if verbose:
            print('  kept {} mutations'.format(len(chosen_cands)))
            print('  score: {:.4f}'.format(best_score))
        assert best_score > old_score
    return template
