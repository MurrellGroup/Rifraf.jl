import random

import numpy as np

# TODO: sparse matrices for banding
# TODO: either write slow code in Cython, or rewrite in Julia
# TODO: use pbdagcon for initial template
# TODO: do not completely recompute As and Bs each iteration. Only recompute from (first) Acol on.
# TODO: compress duplicate reads and account for their scores


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
    result = np.zeros((len(s) + 1, len(t) + 1))
    nrows = min(bandwidth + 1, len(s) + 1)
    ncols = min(bandwidth + 1, len(t) + 1)
    result[:nrows, 0] = log_ins * np.arange(nrows)
    result[0, :ncols] = log_del * np.arange(ncols)
    for i in range(1, len(s) + 1):
        for j in range(max(1, i - bandwidth), min(len(t) + 1, i + bandwidth + 1)):
            result[i, j] = update(result, i, j, j, s[i-1], t[j-1], log_p[i-1], log_ins, log_del, bandwidth)
    return result


def backward(s, log_p, t, log_ins, log_del, bandwidth):
    s = list(reversed(s))
    log_p = list(reversed(log_p))
    t = list(reversed(t))
    return np.flipud(np.fliplr(forward(s, log_p, t, log_ins, log_del, bandwidth)))


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
    Acols = np.zeros((A.shape[0], 2))
    Acols[:, 0] = A[:, pos]
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
    amin, amax = aj - bandwidth, aj + bandwidth + 1
    bcenter = bj - (A.shape[1] - A.shape[0])
    bmin, bmax = bcenter - bandwidth, bcenter + bandwidth + 1
    # need to only consider positions within bandwidth of both columns
    index = range(max(0, amin, bmin), min(A.shape[0], amax, bmax))
    return (Acol + Bcol)[index].max()


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
