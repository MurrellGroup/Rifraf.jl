module Quiver2

include("BandedArray.jl")

using BandedArrayModule

export quiver2

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

function updated_col(mutation::Mutation,
                     template::AbstractString, seq::AbstractString,
                     log_p::Array{Float64, 1}, A::BandedArray{Float64},
                     log_ins::Float64, log_del::Float64)
    aj = mutation.pos + (mutation.m == "substitution" ? 1 : 0)
    prev::SubArray{Float64,1,Array{Float64,2},Tuple{UnitRange{Int64},Int64},2} = sparsecol(A, mutation.pos)
    prev_start, prev_stop = row_range(A, mutation.pos)
    row_start, row_stop = row_range(A, aj)
    offset = 1 - (row_start - prev_start)
    result = Array{Float64}(row_stop - row_start + 1)
    # do first position
    i = 1
    real_i = row_start
    if row_start > prev_start
        # deletion or substitution
        result[1] = max(prev[i - offset] + (mutation.base == seq[real_i - 1] ? 0.0 : log_p[real_i - 1]),
                        prev[i - offset + 1] + log_del)
    else
        if row_start != prev_start
            error("This should not happen")
        end
        # deletion only
        result[1] = prev[i - offset + 1] + log_del
    end
    # do middle positions
    stop = min(prev_stop, row_stop)
    for real_i in (row_start+1):stop
        seq_i = real_i - 1
        i = real_i - row_start + 1
        ii = i - offset + 1
        result[i] = max(result[i - 1] + log_ins,
                        prev[ii - 1] + (mutation.base == seq[seq_i] ? 0.0 : log_p[seq_i]),
                        prev[ii] + log_del)
    end
    # do last position(s)
    for real_i in (stop+1):row_stop
        i = real_i - row_start + 1
        ii = i - offset
        if real_i == prev_stop + 1
            # substitution or insertion
            result[i] = max(result[i - 1] + log_ins,
                            prev[ii] + (mutation.base == seq[real_i - 1] ? 0.0 : log_p[real_i - 1]))
        else
            # insertion only
            if real_i <= prev_stop
                error("This should not happen")
            end
            result[i] = result[i - 1] + log_ins
        end
    end
    return result
end

function equal_ranges(a_range::Tuple{Int64,Int64}, b_range::Tuple{Int64,Int64})
    a_start, a_stop = a_range
    b_start, b_stop = b_range
    alen = a_stop - a_start + 1
    blen = b_stop - b_start + 1
    amin = max(b_start - a_start + 1, 1)
    amax = alen - max(a_stop - b_stop, 0)
    bmin = max(a_start - b_start + 1, 1)
    bmax = blen - max(b_stop - a_stop, 0)
    return (amin, amax), (bmin, bmax)
end

function score_deletion(mutation::Mutation, A::BandedArray{Float64}, B::BandedArray{Float64})
    aj = mutation.pos
    bj = mutation.pos + 1
    Acol = sparsecol(A, aj)
    Bcol = sparsecol(B, bj)
    (amin, amax), (bmin, bmax) = equal_ranges(row_range(A, aj), row_range(B, bj))

    result = Acol[amin] + Bcol[bmin]
    for (i, j) in zip((amin+1):amax, (bmin+1):bmax)
        result = max(result, Acol[i] + Bcol[j])
    end
    return result
end

function score_mutation(mutation::Mutation, template::AbstractString,
                        seq::AbstractString, log_p::Array{Float64, 1},
                        A::BandedArray{Float64}, B::BandedArray{Float64},
                        log_ins::Float64, log_del::Float64, bandwidth::Int)
    if mutation.m == "deletion"
        return score_deletion(mutation, A, B)
    end
    Acol::Array{Float64} = updated_col(mutation, template, seq, log_p, A, log_ins, log_del)
    bj = mutation.pos + (mutation.m == "insertion" ? 0 : 1)
    Bcol::SubArray{Float64,1,Array{Float64,2},Tuple{UnitRange{Int64},Int64},2} = sparsecol(B, bj)

    result = Acol[1] + Bcol[1]
    for i = 2:length(Acol)
        result = max(result, Acol[i] + Bcol[i])
    end
    return result
end

function mutations(template::AbstractString)
    for j in 1:length(template)
        for base in "ACGT"
            if template[j] == base
                continue
            end
            produce(Mutation("substitution", j, base))
        end
        produce(Mutation("deletion", j, 'X'))  # FIXME: how to represent no base here?
        for base in "ACGT"
            produce(Mutation("insertion", j, base))
        end
    end
    # insertion after last
    for base in "ACGT"
        produce(Mutation("insertion", length(template) + 1, base))
    end
end

function update_template(template::AbstractString, mutation::Mutation)
    if mutation.m == "substitution"
        return string(template[1:(mutation.pos - 1)], mutation.base, template[(mutation.pos + 1):end])
    elseif mutation.m == "insertion"
        return string(template[1:(mutation.pos - 1)], mutation.base, template[(mutation.pos):end])
    elseif mutation.m == "deletion"
        return string(template[1:(mutation.pos - 1)], template[(mutation.pos + 1):end])
    else
        # FIXME: solve with type system
        error("unknown mutation: $(mutation.m)")
    end
end

function apply_mutations(template::AbstractString, mutations::Array{Mutation, 1})
    # check that mutations all have different positions.
    # this is a bit too strict, since there are some combinations of
    # mutations affecting the same spot that are unambiguous.
    if length(Set([m.pos for m in mutations])) != length(mutations)
        error("Cannot have multiple mutations affecting the same position")
    end
    remaining = [m for m in mutations]

    while length(remaining) > 0
        m = pop!(remaining)
        template = update_template(template, m)
        o = Dict("insertion" => 1, "deletion" => -1, "substitution" => 0)[m.m]
        for i in 1:length(remaining)
            m2 = remaining[i]
            if m2.pos >= m.pos
                remaining[i] = Mutation(m2.m, m2.pos + o, m2.base)
            end
        end
    end

    return template
end

function choose_candidates(candidates::Vector{MCand}, min_dist::Int)
    final_cands = MCand[]
    posns = Set()
    for c in sort(candidates, by=(c) -> c.score, rev=true)
        if any(Bool[(abs(c.mutation.pos - p) < min_dist) for p in posns])
            continue
        end
        push!(posns, c.mutation.pos)
        push!(final_cands, c)
    end
    return final_cands
end

function quiver2(template::AbstractString, sequences::Vector{ASCIIString},
                 phreds::Vector{Vector{Float64}},
                 log_ins::Float64, log_del::Float64;
                 bandwidth=10::Int, min_dist=9::Int, max_iters=100::Int,
                 verbose=false::Boolean)
    log_ps = [(-phred / 10) for phred in phreds]

    if bandwidth < 0
        error("bandwidth cannot be negative: $bandwidth")
    end

    if verbose
        print("computing initial alignments\n")
    end

    As = [forward(s, p, template, log_ins, log_del, bandwidth) for (s, p) in zip(sequences, log_ps)]
    Bs = [backward(s, p, template, log_ins, log_del, bandwidth) for (s, p) in zip(sequences, log_ps)]
    best_score = sum([A[end, end] for A in As])
    best_template = template
    current_score = best_score
    current_template = best_template
    if verbose
        print("initial score: $best_score\n")
    end
    converged = false
    total_mutations = 0
    iterations = 0
    for i in 1:max_iters
        iterations = i
        old_template = current_template
        old_score = current_score
        if verbose
            print("iteration $i\n")
        end
        candidates = MCand[]
        for mutation in Task(() -> mutations(current_template))
            score = 0
            for si in 1:length(sequences)
                score += score_mutation(mutation, current_template, sequences[si], log_ps[si], As[si], Bs[si], log_ins, log_del, bandwidth)
            end
            if score > current_score
                push!(candidates, MCand(mutation, score))
            end
        end
        if length(candidates) == 0
            converged = true
            break
        end
        if verbose
            print("  found $(length(candidates)) candidate mutations.\n")
        end
        chosen_cands = choose_candidates(candidates, min_dist)
        if verbose
            print("  filtered to $(length(chosen_cands)) candidate mutations\n")
        end
        current_template = apply_mutations(current_template, Mutation[c.mutation for c in chosen_cands])
        As = [forward(s, p, current_template, log_ins, log_del, bandwidth) for (s, p) in zip(sequences, log_ps)]
        Bs = [backward(s, p, current_template, log_ins, log_del, bandwidth) for (s, p) in zip(sequences, log_ps)]
        current_score = sum([A[end, end] for A in As])

        # detect if a single mutation is better
        # FIXME: code duplication
        # FIXME: this may not actually work, because score_mutation() is not exact
        if current_score < chosen_cands[1].score
            if verbose
                print("  rejecting multiple candidates in favor of best\n")
            end
            chosen_cands = MCand[chosen_cands[1]]
            current_template = apply_mutations(old_template, Mutation[c.mutation for c in chosen_cands])
            As = [forward(s, p, current_template, log_ins, log_del, bandwidth) for (s, p) in zip(sequences, log_ps)]
            Bs = [backward(s, p, current_template, log_ins, log_del, bandwidth) for (s, p) in zip(sequences, log_ps)]
            current_score = sum([A[end, end] for A in As])
        end
        total_mutations += length(chosen_cands)
        if verbose
            print("  score: $current_score\n")
        end
        if current_score > best_score
            best_template = current_template
            best_score = current_score
        end
        if old_score > current_score
            if verbose
                # this can happen because score_mutation() is not exact for insertions and deletions
                # FIXME: detect if this keeps happening and return best overall template
                print("Warning: score decreased.\n")
            end
        end
    end
    info = Dict("converged" => converged,
                "iterations" => iterations,
                "mutations" => total_mutations,
                )
    return best_template, info
end

end
