import numpy as np


def forward(s, qv, t, log_ins, log_del):
    result = np.zeros((len(s), len(t)))
    for i in range(len(s)):
        for j in range(len(t)):
            result[i, j] = max([result[i - 1, j] + log_ins,  # insertion
                                result[i, j - 1] + log_del,  # deletion
                                # TODO: phred(1-p) may not be negligible
                                result[i - 1, j - 1] + (0 if s[i] == t[j] else qv[i])])
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


def _quiver2(sequences, phreds, log_ins, log_del, maxiter=100):
    """
    sequences: list of numpy arrays

    phreds: list of numpy array

    """
    log_ps = -phreds
    # choose first sequence as initial template
    # TODO: choose highest quality sequence as a template
    template = np.copy(sequences[0])

    # # compute forward matrices A and backward matrices B
    # As = list(forward(s, p, template, log_ins, log_del) for s, p in zip(sequences, phreds))
    # Bs = list(backward(s, p, template, log_ins, log_del) for s, p in zip(sequences, phreds))

    # iterate: consider all changes and choose best until convergence
    for i in range(maxiter):
        old = template
        best = template
        best_score = score(template, sequences, log_ps, log_ins, log_del)
        for cand in mutate(template):
            cand_score = score(cand, sequences, log_ps, log_ins, log_del)
            if cand_score > best_score:
                best = cand
                best_score = cand_score
        template = best
        if np.all(old == template):
            break
    return template
