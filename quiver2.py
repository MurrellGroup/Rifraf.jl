import numpy as np

# FIXME: insertions/delections can cause wrong answers in
# `updated_col()` and functions that rely on it, because they change
# the shape of the matrix.

# TODO: cache BandedMatrix calls, or precompute ranges, if they take too much time
# TODO: bounds checks in BandedMatrix
# TODO: test BandedMatrix
# TODO: improve slicing and arithmetic operations of BandedMatrix
# TODO: port to Julia and compare speed
# TODO: stochastic gradient descent. weight choices by sequence quality.

# TODO: use faster version of partial order aligner
# TODO: benchmark against actual quiver implementation


class BandedMatrix:
    """A sparse banded matrix.

    Band is always wide enough to include [0, 0] and [-1, -1].

    Currently only supports limited slicing and no arithmetic
    operations.

    Warning: does not do bounds checks.

    """

    def __init__(self, shape, bandwidth=0):
        nrows, ncols = shape
        if not nrows > 0 and ncols > 0:
            raise Exception('shape is not positive: {}'.format(shape))
        self.shape = shape
        self.bandwidth = bandwidth
        self.diff = np.abs(nrows - ncols)
        data_nrows = 2 * bandwidth + self.diff + 1
        self.data = np.zeros((data_nrows, ncols))
        self.h_offset = max(ncols - nrows, 0)
        self.v_offset = max(nrows - ncols, 0)

    @property
    def nrows(self):
        return self.shape[0]

    @property
    def ncols(self):
        return self.shape[1]

    def data_row(self, i, j):
        """Returns k where self.data[k, j] stores M[i, j]"""
        return (i - j) + self.h_offset + self.bandwidth

    def __getitem__(self, key):
        i, j = key
        while i < 0:
            i = self.nrows + i
        while j < 0:
            j = self.ncols + j
        row = self.data_row(i, j)
        return self.data[row, j]

    def row_range(self, j):
        """[start, stop) of band in column j"""
        start = max(0, j - self.h_offset - self.bandwidth)
        stop = min(j + self.v_offset + self.bandwidth + 1, self.nrows)
        return start, stop

    def data_row_range(self, j):
        """[start, stop) of band in self.data column j"""
        start = max(0, self.h_offset + self.bandwidth - j)
        a, b = self.row_range(j)
        stop = start + (b - a)
        return start, stop

    def get_col(self, j):
        """banded part only of column j"""
        dstart, dstop = self.data_row_range(j)
        return self.data[dstart:dstop, j]

    def __setitem__(self, key, val):
        i, j = key
        row = self.data_row(i, j)
        self.data[row, j] = val

    def inband(self, i, j):
        """whether position [i, j] is in the band"""
        row = self.data_row(i, j)
        return 0 <= row < self.data.shape[0]

    def full(self):
        """Construct a dense np.array of this matrix."""
        result = np.zeros((self.nrows, self.ncols))
        for j in range(self.ncols):
            start, stop = self.row_range(j)
            dstart, dstop = self.data_row_range(j)
            result[start:stop, j] = self.data[dstart:dstop, j]
        return result

    def flip(self):
        self.data = np.flipud(np.fliplr(self.data))


def update(arr, i, j, s_base, t_base, log_p, log_ins, log_del):
    """compute alignment score of [i, j]"""
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
    """Forward matrix for aligning `s` to `t`."""
    result = BandedMatrix((len(s) + 1, len(t) + 1), bandwidth)
    for i in range(min(result.shape[0], result.v_offset + bandwidth + 1)):
        result[i, 0] = log_ins * i
    for j in range(min(result.shape[1], result.h_offset + bandwidth + 1)):
        result[0, j] = log_del * j
    # TODO: reverse order of iteration?
    for j in range(1, result.ncols):
        start, stop = result.row_range(j)
        start = max(start, 1)
        for i in range(start, stop):
            result[i, j] = update(result, i, j, s[i-1], t[j-1], log_p[i-1], log_ins, log_del)
    return result


def backward(s, log_p, t, log_ins, log_del, bandwidth):
    """Backward matrix for aligning `s` to `t`."""
    s = ''.join(reversed(s))
    log_p = np.array(list(reversed(log_p)))
    t = ''.join(reversed(t))
    result = forward(s, log_p, t, log_ins, log_del, bandwidth)
    result.flip()
    return result


def updated_col(prev_col, new_col, base, template, seq, log_p, A, log_ins, log_del):
    """Compute a new column for `base` in a particular position.

    prev_col: column from which to calculate new result

    new_col: used only for shape. If this is a substitution, it is
             prev_col + 1; if it is a deletion it is the same as
             prev_col.

    """
    prev = A.get_col(prev_col)
    prev_start, prev_stop = A.row_range(prev_col)
    row_start, row_stop = A.row_range(new_col)
    offset = 1 - (row_start - prev_start)
    result = np.empty((row_stop - row_start))
    for i, real_i in enumerate(range(row_start, row_stop)):
        scores = []
        if i > 0:
            # insertion
            scores.append(result[i - 1] + log_ins)
        if real_i > prev_start and real_i <= prev_stop:
            # match / substitution
            scores.append(prev[i - offset] + (0 if base == seq[real_i-1] else log_p[real_i-1]))
        if real_i >= prev_start and real_i < prev_stop:
            # deletion
            scores.append(prev[i - offset + 1] + log_del)
        result[i] = max(scores)
    return result


def score_mutation(mutation, template, seq, log_p, A, B, log_ins, log_del, bandwidth):
    """Score a mutation using the forward-backward trick."""
    mtype, pos, base = mutation
    bj = pos + (0 if mtype == 'insertion' else 1)
    Bcol = B.get_col(bj)
    if mtype == 'deletion':
        aj = pos
        Acol = A.get_col(aj)
        # need to chop off beginning of Acol and end of Bcol and align
        a_start, a_stop = A.row_range(aj)
        b_start, b_stop = B.row_range(bj)
        amin = max(b_start - a_start, 0)
        amax = len(Acol) - max(a_stop - b_stop, 0)
        bmin = max(a_start - b_start, 0)
        bmax = len(Bcol) - max(b_stop - a_stop, 0)
        return (Acol[amin:amax] + Bcol[bmin:bmax]).max()
    aj = pos + (1 if mtype == 'substitution' else 0)
    Acol = updated_col(pos, aj, base, template, seq, log_p, A, log_ins, log_del)
    return (Acol + Bcol).max()


def mutations(template):
    """Yields all mutations as tuple (name, position, base)"""
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


def update_template(template, mutation):
    """Apply a mutation to a template"""
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
    """Apply multiple mutations to a template.

    Keeps track of changes in positions due to insertions and
    deletions.

    """
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


def choose_candidates(candidates, min_dist):
    """Choose best mutations that are at least `min_dist` bases distant."""
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

    log_ins: penalty for insertion

    log_del: penalty for deletion

    bandwidth: extra bandwidth in banded alignment matrices

    min_dist: minimum distance between multiple mutations accepted in
              a single iteration

    max_iters: maximum iterations before termination

    verbose: print progress to the command line

    """
    log_ps = list((-phred / 10) for phred in phreds)

    if bandwidth < 0:
        raise Exception('bandwidth cannot be negative: {} '.format(bandwidth))

    if verbose:
        print("computing initial alignments")
    As = list(forward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
    Bs = list(backward(s, p, template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
    best_score = sum(A[-1, -1] for A in As)
    best_template = template
    current_score = best_score
    current_template = best_template
    if verbose:
        print('initial score: {:.4f}'.format(best_score))
    converged = False
    total_mutations = 0
    for i in range(max_iters):
        old_template = current_template
        old_score = current_score
        if verbose:
            print("iteration {}".format(i))
        candidates = []
        for mutation in mutations(current_template):
            score = sum(score_mutation(mutation, current_template, seq, log_p, A, B, log_ins, log_del, bandwidth)
                        for seq, log_p, A, B in zip(sequences, log_ps, As, Bs))
            if score > current_score:
                candidates.append([score, mutation])
        if not candidates:
            converged = True
            break
        if verbose:
            print("  found {} candidate mutations".format(len(candidates)))
        chosen_cands = choose_candidates(candidates, min_dist)
        if verbose:
            print('  filtered to {} candidate mutations'.format(len(chosen_cands)))
        current_template = apply_mutations(current_template, chosen_cands)
        As = list(forward(s, p, current_template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
        Bs = list(backward(s, p, current_template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
        current_score = sum(A[-1, -1] for A in As)

        # detect if a single mutation is better
        # FIXME: code duplication
        # FIXME: this may not actually work, because score_mutation() is not exact
        if current_score < chosen_cands[0][0]:
            if verbose:
                print('  rejecting multiple candidates in favor of best')
            chosen_cands = [chosen_cands[0]]
            current_template = apply_mutations(old_template, chosen_cands)
            As = list(forward(s, p, current_template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
            Bs = list(backward(s, p, current_template, log_ins, log_del, bandwidth) for s, p in zip(sequences, log_ps))
            current_score = sum(A[-1, -1] for A in As)
        total_mutations += len(chosen_cands)
        if verbose:
            print('  score: {:.4f}'.format(current_score))
        if current_score > best_score:
            best_template = current_template
            best_score = current_score
        if old_score > current_score:
            if verbose:
                # this can happen because score_mutation() is not exact for insertions and deletions
                # FIXME: detect if this keeps happening and return best overall template
                print('Warning: score decreased.')
    info = {
        'converged': converged,
        'iterations': i,
        'mutations': total_mutations,
    }
    return best_template, info
