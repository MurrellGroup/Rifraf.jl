"""Convenience type for passing around a sequence and its per-base
error probabilities.

"""
type RifrafSequence
    seq::DNASeq
    match_log_p::Vector{LogProb}
    mismatch_log_p::Vector{LogProb}
    ins_log_p::Vector{LogProb}
    del_log_p::Vector{LogProb}
    bandwidth::Int

    function RifrafSequence(seq::DNASeq, error_log_p::Vector{LogProb},
                            bandwidth::Int, scores::Scores)
        if bandwidth < 1
            error("bandwidth must be positive")
        end

        if length(seq) != length(error_log_p)
            error("length mismatch")
        end
        if length(seq) == 0
            return new(seq, LogProb[], LogProb[])
        end
        if minimum(error_log_p) == -Inf
            error("a log error probability is negative infinity")
        end
        if maximum(error_log_p) > 0.0
            bad_value = maximum(error_log_p)
            error("a log error probability is > 0: $bad_value")
        end
        match_log_p = log10(1.0 - exp10(error_log_p))
        mismatch_log_p = error_log_p + scores.mismatch
        ins_log_p = error_log_p + scores.insertion
        del_log_p = Vector{LogProb}(length(error_log_p) + 1)
        del_log_p[1] = error_log_p[1]
        for i=1:length(error_log_p) - 1)
            del_log_p[i + 1] = max(error_log_p[i], error_log_p[i+1])
        end
        del_log_[end] = error_log_p[end]
        return new(seq, match_log_p, mismatch_log_p, ins_log_p, del_log_p, bandwidth)
    end
end

function RifrafSequence(seq::DNASeq, phreds::Vector{Int8}, bandwidth::Int, scores::Scores)
    error_log_p = phred_to_log_p(phreds)
    return RifrafSequence(seq, error_log_p, bandwidth, scores)
end

function length(s::RifrafSequence)
    return length(s.seq)
end
