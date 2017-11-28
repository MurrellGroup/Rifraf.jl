"""Convenience type for passing around a sequence and its per-base
error probabilities.

"""
mutable struct RifrafSequence
    seq::DNASeq
    est_n_errors::Float64
    error_log_p::Vector{LogProb}
    match_scores::Vector{Score}
    mismatch_scores::Vector{Score}
    ins_scores::Vector{Score}
    del_scores::Vector{Score}
    codon_ins_scores::Vector{Score}
    codon_del_scores::Vector{Score}
    bandwidth::Int
    bandwidth_fixed::Bool
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
        return RifrafSequence()
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
    match_scores = Vector{Score}(length(error_log_p))
    mismatch_scores = Vector{Score}(length(error_log_p))
    ins_scores = Vector{Score}(length(error_log_p))
    del_scores = Vector{Score}(length(error_log_p) + 1)

    match_scores[:] = log10.(1.0 - exp10.(error_log_p))
    mismatch_scores[:] = error_log_p + scores.mismatch
    ins_scores[:] = error_log_p + scores.insertion

    del_scores[1] = error_log_p[1] + scores.deletion
    del_scores[end] = error_log_p[end] + scores.deletion
    for i = 1:(length(error_log_p) - 1)
        del_scores[i + 1] = max(error_log_p[i], error_log_p[i + 1]) + scores.deletion
    end

    codon_ins_scores = Vector{Score}()
    codon_del_scores = Vector{Score}()
    if scores.codon_insertion > -Inf
        codon_ins_scores = Vector{Score}(length(error_log_p) - 2)
        for i = 2:(length(error_log_p) - 1)
            codon_ins_scores[i - 1] = max(error_log_p[i - 1],
                                          error_log_p[i],
                                          error_log_p[i + 1]) + scores.codon_insertion
        end
    end
    if scores.codon_deletion > -Inf
        codon_del_scores = Vector{Score}(length(error_log_p) + 1)
        codon_del_scores[1] = error_log_p[1] + scores.codon_deletion
        codon_del_scores[end] = error_log_p[end] + scores.codon_deletion
        for i = 1:(length(error_log_p) - 1)
            codon_del_scores[i + 1] = max(error_log_p[i], error_log_p[i + 1]) + scores.codon_deletion
        end
    end

    n_errors = sum(exp10.(error_log_p))

    return RifrafSequence(seq, n_errors, error_log_p,
                          match_scores, mismatch_scores,
                          ins_scores, del_scores,
                          codon_ins_scores, codon_del_scores,
                          bandwidth, false)
end

"""Given PHRED scores instead of log error rates."""
function RifrafSequence(seq::DNASeq, phreds::Vector{Int8}, bandwidth::Int, scores::Scores)
    error_log_p = phred_to_log_p(phreds)
    return RifrafSequence(seq, error_log_p, bandwidth, scores)
end

"""Update scores for the given sequence."""
function RifrafSequence(seq::RifrafSequence, scores::Scores)
    result = RifrafSequence(seq.seq, seq.error_log_p, seq.bandwidth, scores)
    result.bandwidth_fixed = seq.bandwidth_fixed
    return result
end

"""Empty sequence"""
function RifrafSequence()
    RifrafSequence(DNASeq(), 0.0, LogProb[], Score[], Score[], Score[],
                   Score[], Score[], Score[], 0, false)
end

do_codon_ins(s::RifrafSequence) = length(s.codon_ins_scores) > 0
do_codon_del(s::RifrafSequence) = length(s.codon_del_scores) > 0
do_codon_moves(s::RifrafSequence) = do_codon_ins(s) || do_codon_del(s)

function Base.length(s::RifrafSequence)
    return length(s.seq)
end
