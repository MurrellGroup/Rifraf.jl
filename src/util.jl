const MIN_PHRED = 1
const MAX_PHRED = Int('~') - 33

const BASES = dna"ACGT"
const BASEINTS = Dict(DNA_A => 1, DNA_C => 2, DNA_G => 3, DNA_T => 4)


"""Convert error probability to PHRED score"""
function p_to_phred(p::Float64)
    return Int8(min(round(-10.0 * log10(p)), MAX_PHRED))
end

function p_to_phred(x::Vector{Float64})
    return Int8[p_to_phred(p) for p in x]
end

"""Convert PHRED score to log error probability"""
@generated function phred_to_log_p(x)
    return quote
        return x / (-10.0)
    end
end

"""Convert PHRED score to error probability"""
function phred_to_p(q::Int8)
    return exp10(phred_to_log_p(q))
end

function phred_to_p(x::Vector{Int8})
    return exp10(phred_to_log_p(x))
end


function normalize(parts::Vector{Float64})
    return parts / sum(parts)
end

function cap_phreds(phreds::Vector{Int8}, max_phred::Int8)
    if max_phred < 1
        error("max phred value must be positive")
    end
    return Int8[min(p, max_phred) for p in phreds]
end

"""LogSumExp in base 10.

Borrowed from StatsFuns.jl.
"""
function logsumexp10{T<:Real}(x::AbstractArray{T})
    S = typeof(exp10(zero(T)))    # because of 0.4.0
    isempty(x) && return -S(Inf)
    u = maximum(x)
    abs(u) == Inf && return any(isnan, x) ? S(NaN) : u
    s = zero(S)
    for i = 1:length(x)
        @inbounds s += exp10(x[i] - u)
    end
    log10(s) + u
end
