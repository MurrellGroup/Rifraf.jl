import random

import numpy as np

# TODO: remove sequence <-> array conversion
# TODO: implement banding.
# TODO: either write slow code in Cython, or rewrite in Julia


def seq_to_array(seq):
    convert = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return np.array(list(convert[base] for base in seq))


def array_to_seq(a):
    convert = 'ACGT'
    return ''.join(convert[i] for i in a)


def update(arr, i, j, s_base, t_base, log_p, log_ins, log_del):
    return max([arr[i - 1, j] + log_ins,  # insertion
                arr[i, j - 1] + log_del,  # deletion
                # TODO: log(1-p) may not be negligible
                arr[i - 1, j - 1] + (0 if s_base == t_base else log_p)])


def forward(s, log_p, t, log_ins, log_del):
    result = np.zeros((len(s) + 1, len(t) + 1))
    # invariant: result[i, j] is prob of aligning s[:i] to t[:j]
    result[:, 0] = log_ins * np.arange(len(s) + 1)
    result[0, :] = log_del * np.arange(len(t) + 1)
    for i in range(1, len(s) + 1):
        for j in range(1, len(t) + 1):
            result[i, j] = update(result, i, j, s[i-1], t[j-1], log_p[i-1], log_ins, log_del)
    return result


def backward(s, log_p, t, log_ins, log_del):
    s = list(reversed(s))
    log_p = list(reversed(log_p))
    t = list(reversed(t))
    return np.flipud(np.fliplr(forward(s, log_p, t, log_ins, log_del)))


def mutations(template):
    """Returns (function, position, base)"""
    for j in range(len(template)):
        # mutation
        for base in range(4):
            if template[j] == base:
                continue
            yield (substitution, j, base)
        # deletion
        yield (deletion, j, None)
        # insertion
        for base in range(4):
            yield (insertion, j, base)
    # insertion after last
    for base in range(4):
        yield (insertion, len(template), base)


def substitution(mutation, template, seq_array, log_p, A, B, log_ins, log_del):
    mtype, pos, base = mutation
    Acols = np.copy(A[:, pos:pos+2])
    j = 1
    for i in range(1, A.shape[0]):
        Acols[i, j] = update(Acols, i, j, seq_array[i-1], base, log_p[i-1], log_ins, log_del)
    if pos == len(template) - 1:
        return Acols[:, 1], None
    return Acols[:, 1], B[:, pos + 1]


def deletion(mutation, template, seq_array, log_p, A, B, log_ins, log_del):
    _, pos, _ = mutation
    if pos == len(template) - 1:
        return A[:, -2], None
    Acols = np.zeros((A.shape[0], 2))
    Acols[:, 0] = A[:, pos]
    Acols[0, 1] = Acols[0, 0] + log_del
    mybase = template[pos + 1]
    j = 1
    for i in range(1, A.shape[0]):
        Acols[i, j] = update(Acols, i, j, seq_array[i-1], mybase, log_p[i-1], log_ins, log_del)
    return Acols[:, 0], B[:, pos + 1]


def insertion(mutation, template, seq_array, log_p, A, B, log_ins, log_del):
    _, pos, base = mutation
    Acols = np.zeros((A.shape[0], 2))
    Acols[:, 0] = A[:, pos]
    Acols[0, :] = Acols[0, 0] + np.arange(2) * log_del
    j = 1
    for i in range(1, A.shape[0]):
        Acols[i, j] = update(Acols, i, j, seq_array[i-1], base, log_p[i-1], log_ins, log_del)
    if pos == len(template):
        return Acols[:, 1], None
    return Acols[:, 1], B[:, pos]


def score_mutation(mutation, template, seq_array, log_p, A, B, log_ins, log_del):
    """Score a mutation using the forward-backward trick."""
    f, _, _ = mutation
    Acol, Bcol = f(mutation, template, seq_array, log_p, A, B, log_ins, log_del)
    if Bcol is None:
        return Acol[-1]
    return (Acol + Bcol).max()


def score_slow(template, sequence, log_p, log_ins, log_del):
    return forward(sequence, log_p, template, log_ins, log_del)[-1, -1]


def update_template(template, mutation):
    f, pos, base = mutation
    if f == substitution:
        result = np.copy(template)
        result[pos] = base
        return result
    elif f == insertion:
        return np.insert(template, pos, base)
    elif f == deletion:
        return np.delete(template, pos)
    else:
        raise Exception('unknown mutation: {}'.format(mtype))


def quiver2(sequences, phreds, log_ins, log_del, maxiter=100, seed=None, verbose=False):
    """
    sequences: list of dna sequences

    phreds: list of numpy array

    """
    seq_arrays = list(seq_to_array(s) for s in sequences)
    log_ps = list(-phred / 10 for phred in phreds)
    del sequences
    # TODO: use pbdagcon for initial template
    state = np.random.RandomState(seed)
    template = np.copy(state.choice(seq_arrays))

    As = list(forward(s, p, template, log_ins, log_del) for s, p in zip(seq_arrays, log_ps))
    Bs = list(backward(s, p, template, log_ins, log_del) for s, p in zip(seq_arrays, log_ps))
    best_score = sum(A[-1, -1] for A in As)
    orig_best_score = sum(A[-1, -1] for A in As)
    for i in range(maxiter):
        if verbose:
            print("iteration {}".format(i))
        best_mutation = None
        for mutation in mutations(template):
            score = sum(score_mutation(mutation, template, seq_array, log_p, A, B, log_ins, log_del)
                        for seq_array, log_p, A, B in zip(seq_arrays, log_ps, As, Bs))
            if score > best_score:
                best_mutation = mutation
                best_score = score
        if best_mutation is None:
            break
        if verbose:
            print('score: {}'.format(best_score))
        template = update_template(template, best_mutation)
        As = list(forward(s, p, template, log_ins, log_del) for s, p in zip(seq_arrays, log_ps))
        Bs = list(backward(s, p, template, log_ins, log_del) for s, p in zip(seq_arrays, log_ps))
    return array_to_seq(template)
