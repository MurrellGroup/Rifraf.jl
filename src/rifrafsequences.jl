"""Convenience type for passing around a sequence and its per-base
error probabilities.

"""
type RifrafSequence
    seq::DNASeq
    error_log_p::Vector{LogProb}
    match_log_p::Vector{LogProb}
    bandwidth::Int

    function RifrafSequence(seq::DNASeq, error_log_p::Vector{LogProb},
                            bandwidth::Int)
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
        return new(seq, error_log_p, match_log_p, bandwidth)
    end
end

function RifrafSequence(seq::DNASeq, phreds::Vector{Int8}, bandwidth::Int)
    error_log_p = phred_to_log_p(phreds)
    return RifrafSequence(seq, error_log_p, bandwidth)
end

function length(s::RifrafSequence)
    return length(s.seq)
end

function reverse(s::RifrafSequence)
    return RifrafSequence(reverse(s.seq), reverse(s.error_log_p), s.bandwidth)
end
