"""Convenience type for passing around a sequence and its per-base
error probabilities.

"""
type RifrafSequence
    seq::DNASeq
    match_scores::Vector{LogProb}
    mismatch_scores::Vector{LogProb}
    ins_scores::Vector{LogProb}
    del_scores::Vector{LogProb}
    codon_ins_scores::Vector{LogProb}
    codon_del_scores::Vector{LogProb}
    bandwidth::Int

end

function RifrafSequence(seq::DNASeq, error_log_p::Vector{LogProb},
                        bandwidth::Int, scores::Scores)
    if bandwidth < 1
        error("bandwidth must be positive")
    end

    if length(seq) != length(error_log_p)
        error("length mismatch")
    end
    if length(seq) == 0
        return RifrafSequence(seq, LogProb[], LogProb[], LogProb[],
                              LogProb[], LogProb[], LogProb[], bandwidth)
    end
    if minimum(error_log_p) == -Inf
        error("a log error probability is negative infinity")
    end
    if maximum(error_log_p) > 0.0
        bad_value = maximum(error_log_p)
        error("a log error probability is > 0: $bad_value")
    end
    # all scores are symmetric on the ends so there are fewer
    # branches in the alignment code
    match_scores = Vector{LogProb}(length(error_log_p))
    mismatch_scores = Vector{LogProb}(length(error_log_p))
    ins_scores = Vector{LogProb}(length(error_log_p))
    del_scores = Vector{LogProb}(length(error_log_p) + 1)

    match_scores[:] = log10(1.0 - exp10(error_log_p))
    mismatch_scores[:] = error_log_p + scores.mismatch
    ins_scores[:] = error_log_p + scores.insertion

    del_scores[1] = error_log_p[1] + scores.deletion
    del_scores[end] = error_log_p[end] + scores.deletion
    for i=1:(length(error_log_p) - 1)
        del_scores[i+1] = max(error_log_p[i], error_log_p[i+1]) + scores.deletion
    end

    codon_ins_scores = Vector{LogProb}()
    codon_del_scores = Vector{LogProb}()
    if scores.codon_insertion > -Inf
        codon_ins_scores = Vector{LogProb}(length(error_log_p) - 2)
        for i=2:(length(error_log_p)-1)
            codon_ins_scores[i-1] = max(error_log_p[i - 1],
                                        error_log_p[i],
                                        error_log_p[i + 1]) + scores.codon_insertion
        end
    end
    if scores.codon_deletion > -Inf
        codon_del_scores = Vector{LogProb}(length(error_log_p) + 1)
        codon_del_scores[1] = error_log_p[1] + scores.codon_deletion
        codon_del_scores[end] = error_log_p[end] + scores.codon_deletion
        for i=1:(length(error_log_p) - 1)
            codon_del_scores[i+1] = max(error_log_p[i], error_log_p[i+1]) + scores.codon_deletion
        end
    end

    return RifrafSequence(seq, match_scores, mismatch_scores,
                          ins_scores, del_scores,
                          codon_ins_scores, codon_del_scores,
                          bandwidth)
end

function RifrafSequence(seq::DNASeq, phreds::Vector{Int8}, bandwidth::Int, scores::Scores)
    error_log_p = phred_to_log_p(phreds)
    return RifrafSequence(seq, error_log_p, bandwidth, scores)
end

do_codon_ins(s::RifrafSequence) = length(s.codon_ins_scores) > 0
do_codon_del(s::RifrafSequence) = length(s.codon_del_scores) > 0
do_codon_moves(s::RifrafSequence) = do_codon_ins(s) || do_codon_del(s)

function Base.reverse(s::RifrafSequence)
    return RifrafSequence(reverse(s.seq),
                          reverse(s.match_scores),
                          reverse(s.mismatch_scores),
                          reverse(s.ins_scores),
                          reverse(s.del_scores),
                          reverse(s.codon_ins_scores),
                          reverse(s.codon_del_scores),
                          s.bandwidth)
end

function Base.length(s::RifrafSequence)
    return length(s.seq)
end
