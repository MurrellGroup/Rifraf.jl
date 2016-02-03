import numpy as np


def seq_to_array(seq):
    convert = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return np.array(list(convert[base] for base in seq))


def array_to_seq(a):
    convert = 'ACGT'
    return ''.join(convert[i] for i in a)


def forward(s, phred, t, log_ins, log_del):
    result = np.zeros((len(s) + 1, len(t) + 1))
    # invariant: result[i, j] is prob of aligning s[:i] to t[:j]
    result[:, 0] = log_ins * np.arange(len(s) + 1)
    result[0, :] = log_del * np.arange(len(t) + 1)
    for i in range(1, len(s) + 1):
        for j in range(1, len(t) + 1):
            result[i, j] = max([result[i - 1, j] + log_ins,  # insertion
                                result[i, j - 1] + log_del,  # deletion
                                # TODO: phred(1-p) may not be negligible
                                result[i - 1, j - 1] + (0 if s[i - 1] == t[j - 1] else -phred[i - 1])])
    return result


def backward(s, phred, t, log_ins, log_del):
    s = list(reversed(s))
    phred = list(reversed(phred))
    t = list(reversed(t))
    return np.flipud(np.fliplr(forward(s, phred, t, log_ins, log_del)))


def mutate(template):
    """Returns (type, position, base)"""
    for j in range(len(template)):
        # mutation
        for base in range(4):
            if template[j] == base:
                continue
            yield ('substitution', j, base)
        # deletion
        yield ('deletion', j, None)
        # insertion
        for base in range(4):
            yield ('insertion', j, base)
    # insertion after last
    for base in range(4):
        yield ('insertion', len(template), base)


def substitution(mutation, template, seq_array, phred, A, B, log_ins, log_del):
    mtype, pos, base = mutation
    Acols = np.zeros((A.shape[0], 3))
    Acols[:, 0] = A[:, pos]
    Acols[0, :] = A[0, pos] + np.arange(3) * log_del
    for i in range(1, A.shape[0]):
        for j in (1, 2):
            # only need to update last two columns
            mybase = base if j == 1 else template[pos + 1]
            Acols[i, j] = max([Acols[i - 1, j] + log_ins,  # insertion
                               Acols[i, j - 1] + log_del,  # deletion
                               Acols[i - 1, j - 1] + (0 if seq_array[i - 1] == mybase else -phred[i - 1])])
    return Acols[:, 1:], B[:, pos + 2]


def deletion(mutation, template, seq_array, phred, A, B, log_ins, log_del):
    _, pos, _ = mutation
    Acols = np.zeros((A.shape[0], 2))
    Acols[:, 0] = A[:, pos]
    Acols[0, 1] = Acols[0, 0] + log_del
    mybase = template[pos + 1]
    for i in range(1, A.shape[0]):
        Acols[i, 1] = max([Acols[i - 1, 1] + log_ins,  # insertion
                           Acols[i, 0] + log_del,  # deletion
                           Acols[i - 1, 0] + (0 if seq_array[i - 1] == mybase else -phred[i - 1])])
    return Acols, B[:, pos + 2]


def insertion(mutation, template, seq_array, phred, A, B, log_ins, log_del):
    _, pos, base = mutation
    Acols = np.zeros((A.shape[0], 3))
    Acols[:, 0] = A[:, pos]
    Acols[0, :] = A[0, pos] + np.arange(3) * log_del
    for i in range(1, A.shape[0]):
        for j in (1, 2):
            # only need to update last two columns
            mybase = base if j == 1 else template[pos]
            Acols[i, j] = max([Acols[i - 1, j] + log_ins,  # insertion
                               Acols[i, j - 1] + log_del,  # deletion
                               Acols[i - 1, j - 1] + (0 if seq_array[i - 1] == mybase else -phred[i - 1])])
    return Acols[:, 1:], B[:, pos + 1]


def score_mutation(mutation, template, seq_array, phred, A, B, log_ins, log_del):
    """Score a mutation using the forward-backward trick."""
    mtype, pos, base = mutation
    if mtype == 'substitution':
        if pos == A.shape[1] - 2:
            raise Exception('not implemented yet')
        # Acols contains columns for pos, and pos + 1
        Acols, Bcol = substitution(mutation, template, seq_array, phred, A, B, log_ins, log_del)
    elif mtype == 'insertion':
        pass
    elif mtype == 'deletion':
        pass
    else:
        raise Exception('unknown mutation type: {}'.format(mutation))
    # start with deletion
    result = Acols[0, 1] + Bcol[0]
    for i in range(1, A.shape[0]):
        # all possible ways of combining alignments subalignments
        result = max([result,
                      Acols[i - 1, 1] + Bcol[i],  # insertion
                      Acols[i, 0] + Bcol[i],  # deletion
                      Acols[i - 1, 0] + Bcol[i]])  # match
    return result



def quiver2(sequences, phreds, log_ins, log_del, maxiter=100):
    """
    sequences: list of dna sequences

    phreds: list of numpy array

    """
    seq_arrays = list(seq_to_array(s) for s in sequences)
    # choose first sequence as initial template
    # TODO: choose highest quality sequence as a template
    template = np.copy(seq_arrays[0])

    As = list(forward(s, p, template, log_ins, log_del) for s, p in zip(sequences, phreds))
    Bs = list(backward(s, p, template, log_ins, log_del) for s, p in zip(sequences, phreds))
    best_score = sum(A[-1, -1] for A in As)

    # iterate: consider all changes and choose best until convergence
    for i in range(maxiter):
        best_mutation = None
        for mutation in mutations(template):
            cand_score = sum(score_mutation(mutation, template, seq_array, phred, A, B, log_ins, log_del)
                             for seq_array, phred, A, B in zip(seq_arrays, phreds, As, Bs))
            if cand_score > best_score:
                best_mutation = mutation
                best_score = cand_score
        if best_mutation is None:
            # no better template found
            break
        template = update_template(template, best_mutation)
        As = list(forward(s, p, template, log_ins, log_del) for s, p in zip(sequences, phreds))
        Bs = list(backward(s, p, template, log_ins, log_del) for s, p in zip(sequences, phreds))
        best_score = sum(A[-1, -1] for A in As)
    return array_to_seq(template)
