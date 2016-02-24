module Quiver2

include("BandedArray.jl")

using BandedArrayModule

function update(A::BandedArray{Float64}, i::Int, j::Int, s_base::Char, t_base::Char,
                log_p::Float64, log_ins::Float64, log_del::Float64)
    score = typemin(Float64)
    if inband(A, i - 1, j)
        # insertion
        score = max(score, A[i - 1, j] + log_ins)
    end
    if inband(A, i, j - 1)
        # deletion
        score = max(score, A[i, j - 1] + log_del)
    end
    if inband(A, i - 1, j - 1)
        score = max(score, A[i - 1, j - 1] +( s_base == t_base ? 0.0 : log_p))
    end
    return score
end

function forward(s::AbstractString, log_p::Array{Float64, 1}, t::AbstractString,
                 log_ins::Float64, log_del::Float64, bandwidth::Int)
    result = BandedArray(Float64, (length(s) + 1, length(t) + 1), bandwidth)
    for i = 2:min(size(result)[1], result.v_offset + bandwidth + 1)
        result[i, 1] = log_ins * (i - 1)
    end
    for j = 2:min(size(result)[2], result.h_offset + bandwidth + 1)
        result[1, j] = log_del * (j - 1)
    end
    for j = 2:size(result)[2]
        start, stop = row_range(result, j)
        start = max(start, 2)
        for i = start:stop
            result[i, j] = update(result, i, j, s[i-1], t[j-1], log_p[i-1], log_ins, log_del)
        end
    end
    return result
end

function backward(s::AbstractString, log_p::Array{Float64, 1}, t::AbstractString,
                 log_ins::Float64, log_del::Float64, bandwidth::Int)
    s = reverse(s)
    log_p = flipdim(log_p, 1)
    t = reverse(t)
    result = forward(s, log_p, t, log_ins, log_del, bandwidth)
    return flip(result)
end

end
