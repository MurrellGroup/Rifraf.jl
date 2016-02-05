import numpy as np
import scipy.sparse as sp

from poapy.poagraph import POAGraph
from poapy.seqgraphalignment import SeqGraphAlignment

from _quiver2 import BandedMatrix
from _quiver2 import update
from _quiver2 import forward

# FIXME: asymmetrical bandwidth is wrong
# TODO: use faster version of partial order aligner
# TODO: try column-major order
# TODO: rewrite BandedMatrix as Cython extension type
# TODO: either write slow code in Cython, or rewrite in Julia
# TODO: use pbdagcon for initial template


def backward(s, log_p, t, log_ins, log_del, bandwidth):
    s = ''.join(reversed(s))
    log_p = np.array(list(reversed(log_p)), dtype=np.float32)
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


def updated_col(pos, base, template, seq, log_p, A, B, log_ins, log_del, bandwidth):
    Acols = BandedMatrix((A.shape[0], 2), bandwidth, offset=-pos)
    Acols.data[:, 0] = A.data[:, pos]
    if pos < bandwidth:
        Acols[0, 1] = Acols[0, 0] + log_del
    j = pos + 1
    for i in range(max(1, j - bandwidth), min(A.shape[0], j + bandwidth + 1)):
        Acols[i, 1] = update(Acols, i, 1, j, seq[i-1], base, log_p[i-1], log_ins, log_del, bandwidth)
    return Acols[:, 1]


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
    a_start, a_stop = A.range(aj)
    b_start, b_stop = B.range(bj)

    amin = max(b_start - a_start, 0)
    amax = len(Acol) - max(a_stop - b_stop, 0)

    bmin = max(a_start - b_start, 0)
    bmax = len(Bcol) - max(b_stop - a_stop, 0)

    return (np.asarray(Acol[amin:amax]) + np.asarray(Bcol[bmin:bmax])).max()


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


def quiver2(sequences, phreds, log_ins, log_del, bandwidth=10,
            min_dist=9, max_iters=100, max_multi_iters=50, seed=None,
            verbose=False):
    """Generate an alignment-free consensus.

    sequences: list of dna sequences

    phreds: list of arrays of phred scores

    """
    log_ps = list((-phred / 10).astype(np.float32) for phred in phreds)

    graph = POAGraph(sequences[0])
    for sequence in sequences[1:]:
        alignment = SeqGraphAlignment(sequence, graph, globalAlign=True)
        graph.incorporateSeqAlignment(alignment, sequence)
    _, bases, _ = graph.consensus()
    template = ''.join(bases).upper()

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
        raise Exception('minimum bandwidth is {}, but given {} '.format(_bandwidth))

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
