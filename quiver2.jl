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

immutable Mutation
    m::AbstractString  # TODO: make custom type?
    pos::Int
    base::Char
end

immutable MCand
    mutation::Mutation
    score::Float64
end

function updated_col{T}(prev_col::Int, new_col::Int, base::Char,
                        template::AbstractString, seq::AbstractString,
                        log_p::Array{Float64, 1}, A::BandedArray{T},
                        log_ins::Float64, log_del::Float64)
    prev = sparsecol(A, prev_col)
    prev_start, prev_stop = row_range(A, prev_col)
    row_start, row_stop = row_range(A, new_col)
    offset = 1 - (row_start - prev_start)
    result = Array{T}((row_stop - row_start + 1))
    for real_i in row_start:row_stop
        i = real_i - row_start + 1
        score = typemin(Float64)
        if i > 1
            # insertion
            score = max(score, result[i - 1] + log_ins)
        end
        if real_i > prev_start && real_i <= prev_stop + 1
            # match / substitution
            score = max(score, prev[i - offset] + (base == seq[real_i-1] ? 0 : log_p[real_i-1]))
        end
        if real_i >= prev_start && real_i <= prev_stop
            # deletion
            score = max(score, prev[i - offset + 1] + log_del)
        end
        result[i] = score
    end
    return result
end

function score_mutation{T}(mutation::Mutation, template::AbstractString,
                           seq::AbstractString, log_p::Array{Float64, 1},
                           A::BandedArray{T}, B::BandedArray{T},
                           log_ins::Float64, log_del::Float64, bandwidth::Int)
    bj = mutation.pos + (mutation.m == "insertion" ? 0 : 1)
    Bcol = sparsecol(B, bj)
    if mutation.m == "deletion"
        aj = mutation.pos
        Acol = sparsecol(A, aj)
        # need to chop off beginning of Acol and end of Bcol and align
        a_start, a_stop = row_range(A, aj)
        b_start, b_stop = row_range(B, bj)
        amin = max(b_start - a_start + 1, 1)
        amax = length(Acol) - max(a_stop - b_stop + 1, 0)
        bmin = max(a_start - b_start + 1, 1)
        bmax = length(Bcol) - max(b_stop - a_stop + 1, 0)
        return maximum(Acol[amin:amax] + Bcol[bmin:bmax])
    end
    aj = mutation.pos + (mutation.m == "substitution" ? 1 : 0)
    Acol = updated_col(mutation.pos, aj, mutation.base, template, seq, log_p, A, log_ins, log_del)
    return maximum(Acol + Bcol)
end

# function mutations(template::AbstractString)
#     for j in 1:length(template)
#         for base in "ACGT"
#             if template[j] == base
#                 continue
#             end
#             produce(Mutation("substitution", j, base))
#         end
#         produce(Mutation("deletion", j, None))
#         for base in "ACGT"
#             produce(Mutation("insertion", j, base))
#         end
#     end
#     # insertion after last
#     for base in "ACGT"
#         produce(Mutation("insertion", length(template) + 1, base))
#     end
# end

function update_template(template::AbstractString, mutation::Mutation)
    if mutation.m == "substitution"
        return string(template[1:(mutation.pos - 1)], mutation.base, template[(mutation.pos + 1):end])
    elseif mutation.m == "insertion"
        return string(template[1:(mutation.pos - 1)], mutation.base, template[(mutation.pos):end])
    elseif mutation.m == "deletion"
        return string(template[1:(mutation.pos - 1)], mutation.base, template[(mutation.pos + 1):end])
    else
        # FIXME: solve with type system
        error("unknown mutation")
    end
end

# function apply_mutations!(template::AbstractString, mutations::Array(Mutation, 1))
#     # FIXME: make pure
#     for i in 1:length(mutations):
#         mutation = mutations[i]
#         template = update_template(template, mutation)
#         # TODO: solve with type system
#         offsets = Dict("insertion" => 1, "deletion" => -1, "substitution" => 0)
#         o = offsets[mutation.m]
#         for ii = (i + 1):length(mutations)
#             if mutations[ii].pos > mutation.pos
#                 mutations[ii] = Mutation(mutation.m, mutation.base, mutation..pos + o)
#             end
#         end
#     end
#     return template
# end

# function choose_candidates(candidates::Vector(MCand, 1), min_dist::Int)
#     final_cands = []
#     posns = set()
#     for c in reversed(sorted(candidates, key=lambda c: c[0]))
#         _, (_, posn, _) = c
#         if any(abs(posn - p) < min_dist for p in posns)
#             continue
#         end
#         posns.add(posn)
#         final_cands.append(c)
#     end
#     return final_cands
# end

end
