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

abstract Mutation

immutable Substitution <: Mutation
    pos::Int
    base::Char
end

immutable Insertion <: Mutation
    pos::Int
    base::Char
end

immutable Deletion <: Mutation
    pos::Int
end

immutable CandMutation
    mutation::Mutation
    score::Float64
end

function updated_col(mutation::Mutation,
                     template::AbstractString, seq::AbstractString,
                     log_p::Array{Float64, 1}, A::BandedArray{Float64},
                     log_ins::Float64, log_del::Float64)
    aj = mutation.pos + (typeof(mutation) == Substitution ? 1 : 0)
    prev::SubArray{Float64,1,Array{Float64,2},Tuple{UnitRange{Int64},Int64},2} = sparsecol(A, mutation.pos)
    prev_start, prev_stop = row_range(A, mutation.pos)
    row_start, row_stop = row_range(A, aj)
    offset = 1 - (row_start - prev_start)
    result = Array{Float64}(row_stop - row_start + 1)
    for real_i in row_start:row_stop
        seq_i = real_i - 1
        i = real_i - row_start + 1
        ii = i - offset + 1
        score = typemin(Float64)
        if i > 1
            score = max(score, result[i - 1] + log_ins)
        end
        if prev_start < real_i <= (prev_stop + 1)
            score = max(score, prev[ii - 1] + (mutation.base == seq[seq_i] ? 0.0 : log_p[seq_i]))
        end
        if prev_start <= real_i <= prev_stop
            score = max(score, prev[ii] + log_del)
        end
        result[i] = score
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

@generated function summax(a, b)
    return quote
        result::Float64 = a[1] + b[1]
        for i = 2:min(length(a), length(b))
            result = max(result, a[i] + b[i])
        end
        return result
    end
end

function score_mutation(mutation::Deletion, template::AbstractString,
                        seq::AbstractString, log_p::Array{Float64, 1},
                        A::BandedArray{Float64}, B::BandedArray{Float64},
                        log_ins::Float64, log_del::Float64, bandwidth::Int)
    aj = mutation.pos
    bj = mutation.pos + 1
    Acol = sparsecol(A, aj)
    Bcol = sparsecol(B, bj)
    (amin, amax), (bmin, bmax) = equal_ranges(row_range(A, aj), row_range(B, bj))

    asub = sub(Acol, amin:amax)
    bsub = sub(Bcol, bmin:bmax)

    return summax(asub, bsub)
end

function score_mutation(mutation::Union{Insertion,Substitution}, template::AbstractString,
                        seq::AbstractString, log_p::Array{Float64, 1},
                        A::BandedArray{Float64}, B::BandedArray{Float64},
                        log_ins::Float64, log_del::Float64, bandwidth::Int)
    Acol::Array{Float64} = updated_col(mutation, template, seq, log_p, A, log_ins, log_del)
    bj = mutation.pos + (typeof(mutation) == Insertion ? 0 : 1)
    Bcol::SubArray{Float64,1,Array{Float64,2},Tuple{UnitRange{Int64},Int64},2} = sparsecol(B, bj)

    return summax(Acol, Bcol)
end

function update_template(template::AbstractString, mutation::Substitution)
    return string(template[1:(mutation.pos - 1)], mutation.base, template[(mutation.pos + 1):end])
end

function update_template(template::AbstractString, mutation::Insertion)
    return string(template[1:(mutation.pos - 1)], mutation.base, template[(mutation.pos):end])
end

function update_template(template::AbstractString, mutation::Deletion)
    return string(template[1:(mutation.pos - 1)], template[(mutation.pos + 1):end])
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
        o = Dict(Insertion => 1, Deletion => -1, Substitution => 0)[typeof(m)]
        for i in 1:length(remaining)
            m2 = remaining[i]
            if m2.pos >= m.pos
                T = typeof(m2)
                if T == Deletion
                    remaining[i] = T(m2.pos + o)
                else
                    remaining[i] = T(m2.pos + o, m2.base)
                end
            end
        end
    end

    return template
end

function choose_candidates(candidates::Vector{CandMutation}, min_dist::Int)
    final_cands = CandMutation[]
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

# This macro is used in getcands() for performance reasons, to avoid
# using a variable of type `Mutation`.
#
# TODO: convert this to a generated function.
macro process_mutation(m)
    return quote
        score = 0.0
        for si in 1:length(sequences)
            score += score_mutation($m, template, sequences[si], log_ps[si],
                                    As[si], Bs[si], log_ins, log_del, bandwidth)
        end
        if score > current_score
            push!(candidates, CandMutation($m, score))
        end
    end
end

function getcands(template::AbstractString, current_score::Float64,
                  sequences::Vector{ASCIIString},
                  log_ps::Vector{Vector{Float64}},
                  As::Vector{BandedArray{Float64}}, Bs::Vector{BandedArray{Float64}},
                  log_ins::Float64, log_del::Float64, bandwidth::Int)
    candidates = CandMutation[]
    for j in 1:length(template)
        for base in "ACGT"
            mi = Insertion(j, base)
            @process_mutation mi
            if template[j] == base
                continue
            end
            ms = Substitution(j, base)
            @process_mutation ms
        end
        md = Deletion(j)
        @process_mutation md
    end
    # insertion after last
    for base in "ACGT"
        mi = Insertion(length(template) + 1, base)
        @process_mutation mi
    end
    return candidates
end

function quiver2(template::AbstractString, sequences::Vector{ASCIIString},
                 phreds::Vector{Vector{Float64}},
                 log_ins::Float64, log_del::Float64;
                 bandwidth::Int=10, min_dist::Int=9, batch::Int=10,
                 max_iters::Int=100, verbose::Bool=false)
    log_ps = [(-phred / 10.0) for phred in phreds]

    if bandwidth < 0
        error("bandwidth cannot be negative: $bandwidth")
    end

    if verbose
        print("computing initial alignments\n")
    end

    if batch < 0
        batch = length(sequences)
    end
    batch = min(batch, length(sequences))
    if batch < length(sequences)
        indices = rand(1:length(sequences), batch)
        seqs = sequences[indices]
        lps = log_ps[indices]
    else
        seqs = sequences
        lps = log_ps
    end

    As = [forward(s, p, template, log_ins, log_del, bandwidth) for (s, p) in zip(seqs, lps)]
    Bs = [backward(s, p, template, log_ins, log_del, bandwidth) for (s, p) in zip(seqs, lps)]
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

        candidates = getcands(current_template, current_score, seqs, lps, As, Bs, log_ins, log_del, bandwidth)
        if length(candidates) == 0
            if batch < length(sequences)
                if verbose
                    print("no candidates found. switching off batch mode.\n")
                end
                # start full runs
                batch = length(sequences)
                # TODO: this is necessary because unbatched scores will be worse.
                # what's the alternative?
                best_score = typemin(Float64)
                # TODO: instead of turning off batch mode, try increasing batch size
                # TODO: is there some fast way to detect convergence without a full run?
                # TODO: try multiple iterations before changing/disabling batch
            else
                converged = true
                break
            end
        end
        if length(candidates) > 0
            if verbose
                print("  found $(length(candidates)) candidate mutations.\n")
            end
            chosen_cands = choose_candidates(candidates, min_dist)
            if verbose
                print("  filtered to $(length(chosen_cands)) candidate mutations\n")
            end
            current_template = apply_mutations(current_template, Mutation[c.mutation for c in chosen_cands])
            current_score = 0.0
            for i = 1:length(seqs)
                current_score += forward(seqs[i], lps[i], current_template, log_ins, log_del, bandwidth)[end, end]
            end

            # detect if a single mutation is better
            # FIXME: code duplication
            # FIXME: this may not actually work, because score_mutation() is not exact
            if current_score < chosen_cands[1].score
                if verbose
                    print("  rejecting multiple candidates in favor of best\n")
                end
                chosen_cands = CandMutation[chosen_cands[1]]
                current_template = apply_mutations(old_template, Mutation[c.mutation for c in chosen_cands])
                current_score = 0.0
                for i = 1:length(seqs)
                    current_score += forward(seqs[i], lps[i], current_template, log_ins, log_del, bandwidth)[end, end]
                end
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
        if batch < length(sequences)
            indices = rand(1:length(sequences), batch)
            seqs = sequences[indices]
            lps = log_ps[indices]
        else
            seqs = sequences
            lps = log_ps
        end
        As = [forward(s, p, current_template, log_ins, log_del, bandwidth) for (s, p) in zip(seqs, lps)]
        Bs = [backward(s, p, current_template, log_ins, log_del, bandwidth) for (s, p) in zip(seqs, lps)]
        current_score = 0.0
        for i = 1:length(seqs)
            current_score += As[i][end, end]
        end

    end
    info = Dict("converged" => converged,
                "iterations" => iterations,
                "mutations" => total_mutations,
                )
    return best_template, info
end

end
