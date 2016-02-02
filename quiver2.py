import numpy as np


def seq_to_array(seq):
    convert = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return np.array(list(convert[base] for base in seq))


def array_to_seq(a):
    convert = 'ACGT'
    return ''.join(convert[i] for i in a)


def forward(s, qv, t, log_ins, log_del):
    result = np.zeros((len(s) + 1, len(t) + 1))
    # invariant: result[i, j] is prob of aligning s[:i] to t[:j]
    result[:, 0] = log_ins * np.arange(len(s) + 1)
    result[0, :] = log_del * np.arange(len(t) + 1)
    for i in range(1, len(s) + 1):
        for j in range(1, len(t) + 1):
            result[i, j] = max([result[i - 1, j] + log_ins,  # insertion
                                result[i, j - 1] + log_del,  # deletion
                                # TODO: phred(1-p) may not be negligible
                                result[i - 1, j - 1] + (0 if s[i - 1] == t[j - 1] else qv[i - 1])])
    return result


def score(template, sequences, phreds, log_ins, log_del):
    return sum(forward(s, qv, template, log_ins, log_del)[-1, -1]
               for s, qv in zip(sequences, phreds))


# def backward(s, qv, t, log_ins, log_del):
#     return np.flipud(np.fliplr(forward(s[:,:,-1], qv[:,:,-1], t[:,:,-1], log_ins, log_del)))


# def _mutation_score(j, val, s, t, A, B):
#     if val == t[j]:
#         # no need to update
#         return B[0, 0]
#     if j == 0:
#         # update B[:, 0] and use it
#         for i in reversed(range(len(A))):
#             B[i, j] = max([result[i - 1, j] + log_ins,  # insertion
#                            result[i, j - 1] + log_del,  # deletion
#                            # TODO: phred(1-p) may not be negligible
#                            result[i - 1, j - 1] + (0 if s[i] == t[j] else qv[i])])
#         return B[0, 0]

#     if j == len(A) - 1:
#         # update A[:, j] and use it
#         for i in range(len(A)):
#             A[i, j] = max([result[i - 1, j] + log_ins,  # insertion
#                            result[i, j - 1] + log_del,  # deletion
#                            # TODO: phred(1-p) may not be negligible
#                            result[i - 1, j - 1] + (0 if s[i] == val else qv[i])])
#         x, y = A.shape
#         return A[x - 1, y - 1]
#     # need updated A[:, j]
#     # need un-updated B[:, j+1]

def mutate(template):
    # yield all single indel variants of a sequences
    for j in range(len(template)):
        # mutation
        for val in range(4):
            result = np.copy(template)
            result[j] = val
            yield result
        # deletion
        result = np.copy(template)
        yield np.delete(result, j)
        # insertion
        for val in range(4):
            result = np.copy(template)
            result = np.insert(result, j, val)
            yield result
    # insertion after last
    result = np.copy(template)
    result = np.insert(result, len(result), val)
    yield result


def quiver2(sequences, phreds, log_ins, log_del, maxiter=100):
    """
    sequences: list of dna sequences

    phreds: list of numpy array

    """
    seq_arrays = list(seq_to_array(s) for s in sequences)
    log_ps = list(-p for p in phreds)
    # choose first sequence as initial template
    # TODO: choose highest quality sequence as a template
    template = np.copy(seq_arrays[0])

    # # compute forward matrices A and backward matrices B
    # As = list(forward(s, p, template, log_ins, log_del) for s, p in zip(sequences, phreds))
    # Bs = list(backward(s, p, template, log_ins, log_del) for s, p in zip(sequences, phreds))

    # iterate: consider all changes and choose best until convergence
    for i in range(maxiter):
        old = template
        best = template
        best_score = score(template, seq_arrays, log_ps, log_ins, log_del)
        # TODO: score candidates with column trick
        for cand in mutate(template):
            cand_score = score(cand, seq_arrays, log_ps, log_ins, log_del)
            if cand_score > best_score:
                best = cand
                best_score = cand_score
        template = best
        if np.all(old == template):
            break
    return array_to_seq(template)
