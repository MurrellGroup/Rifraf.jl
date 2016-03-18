# TODO: documentation

module Model

using Bio.Seq

using Quiver2.BandedArrays
using Quiver2.Mutations

export quiver2

function update(A::BandedArray{Float64}, i::Int, j::Int,
                s_base::Char, t_base::Char,
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

function forward(t::AbstractString, s::AbstractString, log_p::Vector{Float64},
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
            result[i, j] = update(result, i, j, s[i-1], t[j-1],
                                  log_p[i-1], log_ins, log_del)
        end
    end
    return result
end

function backward(t::AbstractString, s::AbstractString, log_p::Vector{Float64},
                  log_ins::Float64, log_del::Float64, bandwidth::Int)
    s = reverse(s)
    log_p = flipdim(log_p, 1)
    t = reverse(t)
    result = forward(t, s, log_p, log_ins, log_del, bandwidth)
    return flip(result)
end

function equal_ranges(a_range::Tuple{Int64,Int64},
                      b_range::Tuple{Int64,Int64})
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
                        seq::AbstractString, log_p::Vector{Float64},
                        A::BandedArray{Float64}, B::BandedArray{Float64},
                        log_ins::Float64, log_del::Float64, bandwidth::Int)
    aj = mutation.pos
    bj = mutation.pos + 1
    Acol = sparsecol(A, aj)
    Bcol = sparsecol(B, bj)
    (amin, amax), (bmin, bmax) = equal_ranges(row_range(A, aj),
                                              row_range(B, bj))

    asub = sub(Acol, amin:amax)
    bsub = sub(Bcol, bmin:bmax)

    return summax(asub, bsub)
end

function score_mutation(mutation::Union{Insertion,Substitution},
                        template::AbstractString,
                        seq::AbstractString, log_p::Vector{Float64},
                        A::BandedArray{Float64}, B::BandedArray{Float64},
                        log_ins::Float64, log_del::Float64, bandwidth::Int)
    bj = mutation.pos + (typeof(mutation) == Insertion ? 0 : 1)
    Bcol::SubArray{Float64,1,Array{Float64,2},Tuple{UnitRange{Int64},Int64},2} = sparsecol(B, bj)
    aj = mutation.pos + (typeof(mutation) == Substitution ? 1 : 0)
    prev::SubArray{Float64,1,Array{Float64,2},Tuple{UnitRange{Int64},Int64},2} = sparsecol(A, mutation.pos)
    prev_start, prev_stop = row_range(A, mutation.pos)
    row_start, row_stop = row_range(A, aj)
    offset = 1 - (row_start - prev_start)
    result::Float64 = typemin(Float64)
    # for efficiency, do not allocate and compute complete column of
    # A. Just keep running score.  `prev_score` is the value of
    # Acol[i - 1], for the insertion move
    prev_score::Float64 = 0.0
    for real_i in row_start:row_stop
        seq_i = real_i - 1
        i = real_i - row_start + 1
        ii = i - offset + 1
        score = typemin(Float64)
        if i > 1
            score = max(score, prev_score + log_ins)
        end
        if prev_start < real_i <= (prev_stop + 1)
            score = max(score, prev[ii - 1] + (mutation.base == seq[seq_i] ? 0.0 : log_p[seq_i]))
        end
        if prev_start <= real_i <= prev_stop
            score = max(score, prev[ii] + log_del)
        end
        result = max(result, score + Bcol[i])
        prev_score = score
    end
    return result
end

function update_template(template::AbstractString,
                         mutation::Substitution)
    return string(template[1:(mutation.pos - 1)],
                  mutation.base,
                  template[(mutation.pos + 1):end])
end

function update_template(template::AbstractString,
                         mutation::Insertion)
    return string(template[1:(mutation.pos - 1)],
                  mutation.base,
                  template[(mutation.pos):end])
end

function update_template(template::AbstractString,
                         mutation::Deletion)
    return string(template[1:(mutation.pos - 1)],
                  template[(mutation.pos + 1):end])
end

function apply_mutations(template::AbstractString,
                         mutations::Vector{Mutation})
    # check that mutations all have different positions. this is too
    # strict, since there are some combinations of mutations affecting
    # the same spot that are unambiguous, but combined with `min_dist`
    # it is fine.
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
                  As::Vector{BandedArray{Float64}},
                  Bs::Vector{BandedArray{Float64}},
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
                 log_ps::Vector{Vector{Float64}},
                 log_ins::Float64, log_del::Float64;
                 bandwidth::Int=10, min_dist::Int=9, batch::Int=10,
                 max_iters::Int=100, verbose::Bool=false)
    if bandwidth < 0
        error("bandwidth cannot be negative: $bandwidth")
    end

    if verbose
        print(STDERR, "computing initial alignments\n")
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

    As = [forward(template, s, p, log_ins, log_del, bandwidth)
          for (s, p) in zip(seqs, lps)]
    Bs = [backward(template, s, p, log_ins, log_del, bandwidth)
          for (s, p) in zip(seqs, lps)]
    current_score = sum([A[end, end] for A in As])
    current_template = template
    if verbose
        print(STDERR, "initial score: $current_score\n")
    end
    converged = false
    mutations = Int[]
    batch_iterations = 0
    full_iterations = 0
    for i in 1:max_iters
        if batch < length(sequences)
            batch_iterations +=1
        else
            full_iterations += 1
        end
        old_template = current_template
        old_score = current_score
        if verbose
            print(STDERR, "iteration $i\n")
        end

        # FIXME: this conversion should not be necessary
        As = convert(Vector{BandedArray{Float64}}, As)
        Bs = convert(Vector{BandedArray{Float64}}, Bs)
        candidates = getcands(current_template, current_score, seqs, lps, As, Bs,
                              log_ins, log_del, bandwidth)
        recompute_As = false
        if length(candidates) == 0
            push!(mutations, 0)
            if batch < length(sequences)
                if verbose
                    print(STDERR, "no candidates found. switching off batch mode.\n")
                end
                # start full runs
                batch = length(sequences)
                # TODO: instead of turning off batch mode, try increasing batch size
                # TODO: is there some fast way to detect convergence w/o full run?
                # TODO: try multiple iterations before changing/disabling batch
            else
                converged = true
                break
            end
        else
            if verbose
                print(STDERR, "  found $(length(candidates)) candidate mutations.\n")
            end
            chosen_cands = choose_candidates(candidates, min_dist)
            if verbose
                print(STDERR, "  filtered to $(length(chosen_cands)) candidate mutations\n")
            end
            current_template = apply_mutations(old_template,
                                            Mutation[c.mutation
                                                     for c in chosen_cands])
            As = [forward(current_template, s, p, log_ins, log_del, bandwidth)[end, end]
                  for (s, p) in zip(seqs, lps)]
            temp_score = sum([A[end, end] for A in As])

            # detect if a single mutation is better
            # note: this may not always be correct, because score_mutation() is not exact
            if temp_score < chosen_cands[1].score
                if verbose
                    print(STDERR, "  rejecting multiple candidates in favor of best\n")
                end
                chosen_cands = CandMutation[chosen_cands[1]]
                current_template = apply_mutations(old_template,
                                                   Mutation[c.mutation
                                                            for c in chosen_cands])
                recompute_As = true
            end
            push!(mutations, length(chosen_cands))
        end
        if batch < length(sequences)
            indices = rand(1:length(sequences), batch)
            seqs = sequences[indices]
            lps = log_ps[indices]
            recompute_As = true
        else
            seqs = sequences
            lps = log_ps
        end
        if recompute_As
            As = [forward(current_template, s, p, log_ins, log_del, bandwidth)
                  for (s, p) in zip(seqs, lps)]
        end
        Bs = [backward(current_template, s, p, log_ins, log_del, bandwidth)
              for (s, p) in zip(seqs, lps)]
        current_score = sum([A[end, end] for A in As])
        if verbose
            print(STDERR, "  score: $current_score\n")
        end
    end
    info = Dict("converged" => converged,
                "batch_iterations" => batch_iterations,
                "full_iterations" => full_iterations,
                "mutations" => mutations,
                )
    return current_template, info
end

"""
Alternate quiver2() using BioJulia types.

"""
function quiver2{T<:NucleotideSequence}(template::DNASequence,
                                        sequences::Vector{T},
                                        log_ps::Vector{Vector{Float64}},
                                        log_ins::Float64, log_del::Float64;
                                        bandwidth::Int=10, min_dist::Int=9, batch::Int=10,
                                        max_iters::Int=100, verbose::Bool=false)
    new_template = convert(AbstractString, template)
    new_sequences = ASCIIString[convert(AbstractString, s) for s in sequences]
    result, info = quiver2(new_template, new_sequences, log_ps,
                           log_ins, log_del,
                           bandwidth=bandwidth, min_dist=min_dist, batch=batch,
                           max_iters=max_iters, verbose=verbose)
    return DNASequence(result), info
end

end
