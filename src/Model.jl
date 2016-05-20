# TODO: documentation

# TODO: annealing schedule for log_ins and log_del penalties for
# aligning template to reference

# TODO: backtrace or sample reference alignments. If they contain
# length-1 or length-2 indels, increase penalty and keep running.

# FIXME: swap log_ins and log_del penalties for aligning template to
# reference

module Model

using Bio.Seq

using Quiver2.Sample
using Quiver2.BandedArrays
using Quiver2.Mutations

export quiver2

function update(A::BandedArray{Float64}, i::Int, j::Int,
                s_base::Char, t_base::Char,
                log_p::Float64, log_ins::Float64, log_del::Float64,
                allow_codon_indels::Bool,
                log_codon_ins::Float64, log_codon_del::Float64)
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
    if allow_codon_indels
        if i > 3  && inband(A, i - 3, j)
            # codon insertion
            score = max(score, A[i - 3, j] + log_codon_ins)
        end
        if j > 3 && inband(A, i, j - 3)
            # codon deletion
            score = max(score, A[i, j - 3] + log_codon_del)
        end
    end
    return score
end

"""
F[i, j] is the log probability of aligning s[1:i-1] to t[1:j-1].

"""
function forward_codon(t::AbstractString, s::AbstractString, log_p::Vector{Float64},
                       log_ins::Float64, log_del::Float64, bandwidth::Int,
                       allow_codon_indels::Bool,
                       log_codon_ins::Float64, log_codon_del::Float64)
    if length(s) != length(log_p)
        error("sequence length does not match quality score length")
    end
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
                                  log_p[i-1], log_ins, log_del,
                                  allow_codon_indels,
                                  log_codon_ins, log_codon_del)
        end
    end
    return result::BandedArray{Float64}
end

function forward(t::AbstractString, s::AbstractString, log_p::Vector{Float64},
                 log_ins::Float64, log_del::Float64, bandwidth::Int)
    return forward_codon(t, s, log_p, log_ins, log_del, bandwidth, false, 0.0, 0.0)
end

"""
B[i, j] is the log probability of aligning s[i:end] to t[j:end].

"""
function backward_codon(t::AbstractString, s::AbstractString, log_p::Vector{Float64},
                        log_ins::Float64, log_del::Float64, bandwidth::Int,
                        allow_codon_indels::Bool,
                        log_codon_ins::Float64, log_codon_del::Float64)
    s = reverse(s)
    log_p = flipdim(log_p, 1)
    t = reverse(t)
    result = forward_codon(t, s, log_p, log_ins, log_del, bandwidth,
                           allow_codon_indels, log_codon_ins, log_codon_del)
    return flip(result)
end

function backward(t::AbstractString, s::AbstractString, log_p::Vector{Float64},
                  log_ins::Float64, log_del::Float64, bandwidth::Int)
    return backward_codon(t, s, log_p, log_ins, log_del, bandwidth, false, 0.0, 0.0)
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

function score_mutation(mutation::Union{Deletion,CodonDeletion},
                        template::AbstractString,
                        seq::AbstractString, log_p::Vector{Float64},
                        A::BandedArray{Float64}, B::BandedArray{Float64},
                        log_ins::Float64, log_del::Float64, bandwidth::Int,
                        allow_codon_indels::Bool,
                        log_codon_ins::Float64, log_codon_del::Float64)
    offset = typeof(mutation) == Deletion ? 1 : 3
    aj = mutation.pos           # column before base to delete (new last column before change)
    bj = mutation.pos + offset  # column after bases to delete
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
                        log_ins::Float64, log_del::Float64, bandwidth::Int,
                        allow_codon_indels::Bool,
                        log_codon_ins::Float64, log_codon_del::Float64)
    bj = mutation.pos + 1 # column after base to mutate or insert
    Bcol::SubArray{Float64,1,Array{Float64,2},Tuple{UnitRange{Int64},Int64},2} = sparsecol(B, bj)
    b_start, b_stop = row_range(B, bj)  # row range of B column

    aj = mutation.pos + 1  # index of new column to compute

    ajprev = mutation.pos + (typeof(mutation) == Substitution ? 0 : 1)  # index of last unchanged column in A
    prev_start, prev_stop = row_range(A, ajprev)  # row range of column before base to mutate or insert

    ajprev3 = ajprev - 2
    prev_start::Int
    prev_stop::Int
    if ajprev3 > 0
        prev3_start, prev3_stop = row_range(A, ajprev3)
    end

    row_start, row_stop = row_range(A, aj)  # row range of column to compute
    offset = 1 - (row_start - prev_start)
    result::Float64 = typemin(Float64)

    if (row_start != b_start || row_stop != b_stop)
        error("Acol does not align to Bcol")
    end

    # for efficiency, do not allocate and compute complete column of
    # A. Just keep last four positions.
    # `scores[1]` is the score of the current entry
    # `scores[4]` is the score of the entry one codon above the current entry
    scores = [0.0, 0.0, 0.0, 0.0]
    for real_i in row_start:row_stop
        seq_i = real_i - 1  # position in the sequence
        i = real_i - row_start + 1  # position in this iteration
        prev_i = i - offset + 1 # position in `prev` matching i
        scores[1] = typemin(Float64)
        if i > 1
            # insertion
            scores[1] = max(scores[1], scores[2] + log_ins)
        end
        if prev_start < real_i <= (prev_stop + 1)
            # (mis)match
            scores[1] = max(scores[1], A[real_i - 1, ajprev] + (mutation.base == seq[seq_i] ? 0.0 : log_p[seq_i]))
        end
        if prev_start <= real_i <= prev_stop
            # deletion
            scores[1] = max(scores[1], A[real_i, ajprev] + log_del)
        end
        if allow_codon_indels
            if i > 3
                # insertion
                scores[1] = max(scores[1], scores[4] + log_codon_ins)
            end
            if ajprev3 > 0 && (prev3_start <= real_i <= prev3_stop) && inband(A, real_i, ajprev3)
                # deletion
                scores[1] = max(scores[1], A[real_i, ajprev-2] + log_codon_del)
            end
        end
        b_i = i  #  this is okay because we've chosen aj so it has the same range as bj
        result = max(result, scores[1] + Bcol[b_i])
        scores[2:4] = scores[1:3]
    end
    return result
end


function compute_subcol(row_start, row_stop,
                        prev, prev_start, prev_stop,
                        prev3, prev3_start, prev3_stop,
                        has_prev3,
                        base, seq, log_p,
                        log_ins, log_del,
                        allow_codon_indels,
                        log_codon_ins, log_codon_del)
    col = typemin(Float64) * ones(row_stop - row_start + 1)
    offset = 1 - (row_start - prev_start)
    for i in 1:length(col)
        real_i = i + row_start - 1
        seq_i = real_i - 1  # position in the sequence
        prev_i = i - offset + 1 # position in `prev` matching i
        if i > 1
            # insertion
            col[i] = max(col[i], col[i-1] + log_ins)
        end
        if prev_start < real_i <= (prev_stop + 1)
            # (mis)match
            col[i] = max(col[i], prev[prev_i-1] + (base == seq[seq_i] ? 0.0 : log_p[seq_i]))
        end
        if prev_start <= real_i <= prev_stop
            # deletion
            col[i] = max(col[i], prev[prev_i] + log_del)
        end
        if allow_codon_indels
            if i > 3
                # insertion
                col[i] = max(col[i], col[i-3] + log_codon_ins)
            end
            if has_prev3 && prev3_start <= real_i <= prev3_stop
                # deletion
                prev3_i = real_i - prev3_start + 1
                col[i] = max(col[i], prev3[prev3_i] + log_codon_del)
            end
        end
    end
    return col
end

function score_mutation(mutation::CodonInsertion,
                        template::AbstractString,
                        seq::AbstractString, log_p::Vector{Float64},
                        A::BandedArray{Float64}, B::BandedArray{Float64},
                        log_ins::Float64, log_del::Float64, bandwidth::Int,
                        allow_codon_indels::Bool,
                        log_codon_ins::Float64, log_codon_del::Float64)
    bj = mutation.pos + 1 # column after base to mutate or insert
    Bcol::SubArray{Float64,1,Array{Float64,2},Tuple{UnitRange{Int64},Int64},2} = sparsecol(B, bj)
    b_start, b_stop = row_range(B, bj)  # row range of B column

    aj = mutation.pos + 1  # index of new column to compute
    row_start, row_stop = row_range(A, aj)  # row range of column to compute
    if (row_start != b_start || row_stop != b_stop)
        error("Acol does not align to Bcol")
    end
    nrows = row_stop - row_start + 1
    Acols = Array(Float64, (nrows, 3))

    ajprev = aj
    prev::SubArray{Float64,1,Array{Float64,2},Tuple{UnitRange{Int64},Int64},2} = sparsecol(A, ajprev)
    prev_start, prev_stop = row_range(A, ajprev)  # row range of column before base to mutate or insert

    has_prev3 = false
    ajprev3 = ajprev - 2
    prev3::SubArray{Float64,1,Array{Float64,2},Tuple{UnitRange{Int64},Int64},2} = prev
    prev3_start::Int = -1
    prev3_stop::Int = -1
    if ajprev3 > 0
        has_prev3 = true
        prev3 = sparsecol(A, ajprev3)
        prev3_start, prev3_stop = row_range(A, ajprev3)
    end

    Acols[:, 1] = compute_subcol(row_start, row_stop,
                                 prev, prev_start, prev_stop,
                                 prev3, prev3_start, prev3_stop,
                                 has_prev3,
                                 mutation.bases[1], seq, log_p,
                                 log_ins,
                                 log_del,
                                 allow_codon_indels,
                                 log_codon_ins,
                                 log_codon_del)

    ajprev3 = ajprev - 1
    if ajprev3 > 0
        has_prev3 = true
        prev3 = sparsecol(A, ajprev3)
        prev3_start, prev3_stop = row_range(A, ajprev3)
    end
    Acols[:, 2] = compute_subcol(row_start, row_stop,
                                 Acols[:, 1], row_start, row_stop,
                                 prev3, prev3_start, prev3_stop,
                                 has_prev3,
                                 mutation.bases[2], seq, log_p,
                                 log_ins,
                                 log_del,
                                 allow_codon_indels,
                                 log_codon_ins,
                                 log_codon_del)

    ajprev3 = ajprev
    if ajprev3 > 0
        has_prev3 = true
        prev3 = sparsecol(A, ajprev3)
        prev3_start, prev3_stop = row_range(A, ajprev3)
    end
    Acols[:, 3] = compute_subcol(row_start, row_stop,
                                 Acols[:, 2], row_start, row_stop,
                                 prev3, prev3_start, prev3_stop,
                                 has_prev3,
                                 mutation.bases[3], seq, log_p,
                                 log_ins,
                                 log_del,
                                 allow_codon_indels,
                                 log_codon_ins,
                                 log_codon_del)

    result = maximum(Acols[:, end] + Bcol)
    return result
end


function choose_candidates(candidates::Vector{CandMutation}, min_dist::Int)
    final_cands = CandMutation[]
    posns = Set()
    for c in sort(candidates, by=(c) -> c.score, rev=true)
        if any(Bool[(abs(c.mutation.pos - p) < min_dist) for p in posns])
            continue
        end
        union!(posns, affected_positions(c.mutation))
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
                                    As[si], Bs[si], log_ins, log_del, bandwidth,
                                    false, log_codon_ins, log_codon_del)
        end
        if use_ref
            score += score_mutation($m, template, reference, reference_log_p,
                                    A_t, B_t, log_ref_ins, log_ref_del, bandwidth,
                                    true, log_codon_ins, log_codon_del)
        end
        if score > current_score
            push!(candidates, CandMutation($m, score))
        end
    end
end

function getcands(template::AbstractString, current_score::Float64,
                  use_ref::Bool,
                  reference::AbstractString,
                  reference_log_p::Vector{Float64},
                  A_t::BandedArray{Float64},
                  B_t::BandedArray{Float64},
                  sequences::Vector{ASCIIString},
                  log_ps::Vector{Vector{Float64}},
                  As::Vector{BandedArray{Float64}},
                  Bs::Vector{BandedArray{Float64}},
                  log_ins::Float64, log_del::Float64, bandwidth::Int,
                  log_ref_ins::Float64, log_ref_del::Float64,
                  log_codon_ins::Float64, log_codon_del::Float64)
    len = length(template)
    candidates = CandMutation[]
    # insertion at beginning
    for base in "ACGT"
        mi = Insertion(0, base)
        @process_mutation mi
    end
    for j in 1:len
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
    if use_ref
        mci = CodonInsertion(0, random_codon())
        @process_mutation mci
        for j in 1:len
            mci = CodonInsertion(j, random_codon())
            @process_mutation mci
            if j < len - 1
                mcd = CodonDeletion(j)
                @process_mutation mcd
            end
        end
    end
    return candidates
end


function quiver2(reference::AbstractString,
                 template::AbstractString,
                 sequences::Vector{ASCIIString},
                 log_ps::Vector{Vector{Float64}},
                 log_ins::Float64, log_del::Float64;
                 use_ref::Bool=true,
                 bandwidth::Int=10, min_dist::Int=9,
                 batch::Int=10, do_full::Bool=false,
                 max_iters::Int=100, verbose::Bool=false)
    if bandwidth < 0
        error("bandwidth cannot be negative: $bandwidth")
    end

    if verbose
        print(STDERR, "computing initial alignments\n")
    end

    # TODO: do not hardcode
    log_codon_ins = 0.0
    log_codon_del = 0.0
    log_mismatch = 0.0
    log_ref_ins = 0.0
    log_ref_del = 0.0
    if use_ref
        log_codon_ins = -9.0
        log_codon_del = -9.0
        log_mismatch = -3.0
        log_ref_ins = -1.0
        log_ref_del = -1.0
    end
    min_log_ref_ins = typemin(Float64)
    min_log_ref_del = typemin(Float64)

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

    reference_log_p = log_mismatch * ones(length(reference))

    As = [forward(template, s, p, log_ins, log_del, bandwidth)
          for (s, p) in zip(seqs, lps)]
    Bs = [backward(template, s, p, log_ins, log_del, bandwidth)
          for (s, p) in zip(seqs, lps)]
    current_score = sum([A[end, end] for A in As])
    A_t = BandedArray(Float64, (1, 1), 1)
    B_t = BandedArray(Float64, (1, 1), 1)
    if use_ref
        A_t = forward_codon(template, reference, reference_log_p, log_ref_ins, log_ref_del, bandwidth,
                            true, log_codon_ins, log_codon_del)
        B_t = backward_codon(template, reference, reference_log_p, log_ref_ins, log_ref_del, bandwidth,
                             true, log_codon_ins, log_codon_del)
        current_score += A_t[end, end]
    end
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

        candidates = getcands(current_template, current_score,
                              use_ref,
                              reference, reference_log_p, A_t, B_t,
                              seqs, lps, As, Bs,
                              log_ins, log_del, bandwidth,
                              log_ref_ins, log_ref_del,
                              log_codon_ins, log_codon_del)
        recompute_As = false
        if length(candidates) == 0
            push!(mutations, 0)
            if length(current_template) % 3 == 0
                if batch < length(sequences) && do_full
                    if verbose
                        print(STDERR, "no candidates found. switching off batch mode.\n")
                    end
                    # start full runs
                    batch = length(sequences)
                    recompute_As = true
                    # TODO: instead of turning off batch mode, try increasing batch size
                    # TODO: is there some fast way to detect convergence w/o full run?
                    # TODO: try multiple iterations before changing/disabling batch
                else
                    converged = true
                    break
                end
            elseif log_ref_ins > min_log_ref_ins || log_ref_del > min_log_ref_del
                # increase indel penalty
                log_ref_ins = min(log_ref_ins * 2, min_log_ref_ins)
                log_ref_del = min(log_ref_del * 2, min_log_ref_del)
            else
                # give up
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
            As = [forward(current_template, s, p, log_ins, log_del, bandwidth)
                  for (s, p) in zip(seqs, lps)]
            temp_score = sum([A[end, end] for A in As])
            if use_ref
                A_t = forward_codon(current_template, reference, reference_log_p, log_ref_ins, log_ref_del, bandwidth,
                                    true, log_codon_ins, log_codon_del)
                temp_score += A_t[end, end]
            end
            # detect if a single mutation is better
            # note: this may not always be correct, because score_mutation() is not exact
            if temp_score < chosen_cands[1].score
                if verbose
                    print(STDERR, "  rejecting multiple candidates in favor of best\n")
                    print(STDERR, "  $(chosen_cands[1])\n")
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
            if use_ref
                A_t = forward_codon(current_template, reference, reference_log_p, log_ref_ins, log_ref_del, bandwidth,
                                    true, log_codon_ins, log_codon_del)
            end
        end
        Bs = [backward(current_template, s, p, log_ins, log_del, bandwidth)
              for (s, p) in zip(seqs, lps)]
        if use_ref
            B_t = backward_codon(current_template, reference, reference_log_p, log_ref_ins, log_ref_del, bandwidth,
                                 true, log_codon_ins, log_codon_del)
        end
        current_score = sum([A[end, end] for A in As])
        if use_ref
            current_score += A_t[end, end]
        end
        if verbose
            print(STDERR, "  score: $current_score\n")
            print(STDERR, "  consensus: $(length(current_template))\n")
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
function quiver2{T<:NucleotideSequence}(reference::DNASequence,
                                        template::DNASequence,
                                        sequences::Vector{T},
                                        log_ps::Vector{Vector{Float64}},
                                        log_ins::Float64, log_del::Float64;
                                        use_ref::Bool=true,
                                        bandwidth::Int=10, min_dist::Int=9, batch::Int=10,
                                        max_iters::Int=100, verbose::Bool=false)
    new_reference = convert(AbstractString, reference)
    new_template = convert(AbstractString, template)
    new_sequences = ASCIIString[convert(AbstractString, s) for s in sequences]
    result, info = quiver2(new_reference, new_template, new_sequences, log_ps,
                           log_ins, log_del,
                           use_ref=use_ref,
                           bandwidth=bandwidth, min_dist=min_dist, batch=batch,
                           max_iters=max_iters, verbose=verbose)
    return DNASequence(result), info
end

end
