const MIN_PHRED = Phred(1)
const MAX_PHRED = Phred(Int('~') - 33)

"""Convert error probability to PHRED score"""
function p_to_phred(p::Prob)
    return Phred(min(round(-10.0 * log10(p)), MAX_PHRED))
end

function p_to_phred(x::Vector{LogProb})
    return Phred[p_to_phred(p) for p in x]
end

"""Convert PHRED score to log error probability"""
@generated function phred_to_log_p(x)
    return quote
        return x / (-10.0)
    end
end

"""Convert PHRED score to error probability"""
function phred_to_p(q::Phred)
    return exp10(phred_to_log_p(q))
end

function phred_to_p(x::Vector{Phred})
    return exp10.(phred_to_log_p(x))
end


@generated function normalize(parts)
    return quote
        return parts / sum(parts)
    end
end

function cap_phreds(phreds::Vector{Phred}, max_phred::Phred)
    if max_phred < 1
        error("max phred value must be positive")
    end
    return Phred[min(p, max_phred) for p in phreds]
end
