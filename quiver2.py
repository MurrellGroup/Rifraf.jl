import random

import numpy as np


def seq_to_array(seq):
    convert = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return np.array(list(convert[base] for base in seq))


def array_to_seq(a):
    convert = 'ACGT'
    return ''.join(convert[i] for i in a)


def update(arr, i, j, s_base, t_base, phred, log_ins, log_del):
    return max([arr[i - 1, j] + log_ins,  # insertion
                arr[i, j - 1] + log_del,  # deletion
                # TODO: phred(1-p) may not be negligible
                arr[i - 1, j - 1] + (0 if s_base == t_base else -phred / 10)])


def forward(s, phred, t, log_ins, log_del):
    result = np.zeros((len(s) + 1, len(t) + 1))
    # invariant: result[i, j] is prob of aligning s[:i] to t[:j]
    result[:, 0] = log_ins * np.arange(len(s) + 1)
    result[0, :] = log_del * np.arange(len(t) + 1)
    for i in range(1, len(s) + 1):
        for j in range(1, len(t) + 1):
            result[i, j] = update(result, i, j, s[i-1], t[j-1], phred[i-1], log_ins, log_del)
    return result


def backward(s, phred, t, log_ins, log_del):
    s = list(reversed(s))
    phred = list(reversed(phred))
    t = list(reversed(t))
    return np.flipud(np.fliplr(forward(s, phred, t, log_ins, log_del)))


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


def substitution(mutation, template, seq_array, phred, A, B, log_ins, log_del):
    mtype, pos, base = mutation
    Acols = np.copy(A[:, pos:pos+2])
    j = 1
    for i in range(1, A.shape[0]):
        Acols[i, j] = update(Acols, i, j, seq_array[i-1], base, phred[i-1], log_ins, log_del)
    if pos == len(template) - 1:
        return Acols[:, 1], None
    return Acols[:, 1], B[:, pos + 1]


def deletion(mutation, template, seq_array, phred, A, B, log_ins, log_del):
    _, pos, _ = mutation
    if pos == len(template) - 1:
        return A[:, -2], None
    Acols = np.zeros((A.shape[0], 2))
    Acols[:, 0] = A[:, pos]
    Acols[0, 1] = Acols[0, 0] + log_del
    mybase = template[pos + 1]
    j = 1
    for i in range(1, A.shape[0]):
        Acols[i, j] = update(Acols, i, j, seq_array[i-1], mybase, phred[i-1], log_ins, log_del)
    return Acols[:, 0], B[:, pos + 1]


def insertion(mutation, template, seq_array, phred, A, B, log_ins, log_del):
    _, pos, base = mutation
    Acols = np.zeros((A.shape[0], 2))
    Acols[:, 0] = A[:, pos]
    Acols[0, :] = Acols[0, 0] + np.arange(2) * log_del
    j = 1
    for i in range(1, A.shape[0]):
        Acols[i, j] = update(Acols, i, j, seq_array[i-1], base, phred[i-1], log_ins, log_del)
    if pos == len(template):
        return Acols[:, 1], None
    return Acols[:, 1], B[:, pos]


def score_mutation(mutation, template, seq_array, phred, A, B, log_ins, log_del):
    """Score a mutation using the forward-backward trick."""
    f, _, _ = mutation
    Acol, Bcol = f(mutation, template, seq_array, phred, A, B, log_ins, log_del)
    if Bcol is None:
        return Acol[-1]
    return (Acol + Bcol).max()


def score_slow(template, sequence, phred, log_ins, log_del):
    return forward(sequence, phred, template, log_ins, log_del)[-1, -1]


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


def quiver2(sequences, phreds, log_ins, log_del, maxiter=100):
    """
    sequences: list of dna sequences

    phreds: list of numpy array

    """
    seq_arrays = list(seq_to_array(s) for s in sequences)
    del sequences
    # choose random sequence as initial template
    # TODO: use pbdagcon for initial template
    template = np.copy(random.choice(seq_arrays))

    As = list(forward(s, p, template, log_ins, log_del) for s, p in zip(seq_arrays, phreds))
    Bs = list(backward(s, p, template, log_ins, log_del) for s, p in zip(seq_arrays, phreds))
    best_score = sum(A[-1, -1] for A in As)
    orig_best_score = sum(A[-1, -1] for A in As)
    # iterate: consider all changes and choose best until convergence
    for i in range(maxiter):
        best_mutation = None
        for mutation in mutations(template):
            scores = list(score_mutation(mutation, template, seq_array, phred, A, B, log_ins, log_del)
                          for seq_array, phred, A, B in zip(seq_arrays, phreds, As, Bs))
            cand_template = update_template(template, mutation)
            scores2 = list(score_slow(cand_template, s, p, log_ins, log_del) for s, p in zip(seq_arrays, phreds))
            assert np.all(np.array(scores) == np.array(scores2))
            score = sum(scores)
            if score > best_score:
                best_mutation = mutation
                best_score = score
        if best_mutation is None:
            # no better template found
            break
        new_template = update_template(template, best_mutation)
        new_As = list(forward(s, p, new_template, log_ins, log_del) for s, p in zip(seq_arrays, phreds))
        new_Bs = list(backward(s, p, new_template, log_ins, log_del) for s, p in zip(seq_arrays, phreds))
        new_score = sum(a[-1, -1] for a in new_As)
        assert(np.any(list(o[-1, -1] != n[-1, -1] for o, n in zip(As, new_As))))
        assert(np.any(list(o[0, 0] != n[0, 0] for o, n in zip(Bs, new_Bs))))
        assert new_score == best_score
        assert new_score > orig_best_score
        template = new_template
        As = new_As
        Bs = new_Bs
    return array_to_seq(template)
