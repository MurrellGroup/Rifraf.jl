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


@enum DPMove match=1 ins=2 del=3 codon_ins=4 codon_del=5

const offsets = ([-1, -1],
                 [-1, 0],
                 [0, -1],
                 [-3, 0],
                 [0, -3])


function update_helper(A::BandedArray{Float64}, i::Int, j::Int, move::DPMove,
                       penalty::Float64, final_score::Float64, final_move::DPMove)
    offset = offsets[Int(move)]
    i = i + offset[1]
    j = j + offset[2]
    if inband(A, i, j)
        score = A[i, j] + penalty
        if score > final_score
            return score, move
        end
    end
    return final_score, final_move
end


function update(A::BandedArray{Float64}, i::Int, j::Int,
                s_base::Char, t_base::Char,
                log_p::Float64, log_ins::Float64, log_del::Float64,
                allow_codon_indels::Bool,
                log_codon_ins::Float64, log_codon_del::Float64)
    result = (typemin(Float64), match)
    log_match = (s_base == t_base ? 0.0 : log_p)
    result = update_helper(A, i, j, match, log_match, result...)
    result = update_helper(A, i, j, ins, log_ins, result...)
    result = update_helper(A, i, j, del, log_del, result...)
    if allow_codon_indels
        if i > 3
            result = update_helper(A, i, j, codon_ins, log_codon_ins, result...)
        end
        if j > 3
            result = update_helper(A, i, j, codon_del, log_codon_del, result...)
        end
    end
    return result
end


function backtrace(moves::BandedArray{Int}, t::AbstractString, s::AbstractString)
    aligned_t = Char[]
    aligned_s = Char[]
    i, j = moves.shape
    while i > 1 || j > 1
        m = moves[i, j]
        move = DPMove(m)
        si = i - 1
        tj = j - 1
        if move == match
            push!(aligned_t, t[tj])
            push!(aligned_s, s[si])
        elseif move == ins
            push!(aligned_t, '-')
            push!(aligned_s, s[si])
        elseif move == del
            push!(aligned_t, t[tj])
            push!(aligned_s, '-')
        elseif move == codon_ins
            append!(aligned_t, ['-', '-', '-'])
            append!(aligned_s, [s[si], s[si-1], s[si-2]])
        elseif move == codon_del
            append!(aligned_t, [t[tj], t[tj-1], t[tj-2]])
            append!(aligned_s, ['-', '-', '-'])
        end
        offset = offsets[m]
        i += offset[1]
        j += offset[2]
    end
    return join(reverse(aligned_t)), join(reverse(aligned_s))
end


function prepare_array(t::AbstractString, s::AbstractString,
                       log_p::Vector{Float64},
                       log_ins::Float64, log_del::Float64,
                       bandwidth::Int)
    # FIXME: consider codon moves in initialization
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
    return result
end


""" Does some work as forward_codon, but also keeps track of moves
does backtracing to find best alignment.

"""
function align(t::AbstractString, s::AbstractString,
               log_p::Vector{Float64},
               log_ins::Float64, log_del::Float64, bandwidth::Int,
               allow_codon_indels::Bool,
               log_codon_ins::Float64, log_codon_del::Float64)
    # FIXME: code duplication with forward_codon(). This is done in a
    # seperate function to keep return type stable and avoid
    # allocating the `moves` array unnecessarily.
    result = prepare_array(t, s, log_p, log_ins, log_del, bandwidth)
    moves = BandedArray(Int, result.shape, bandwidth)
    for i = 2:min(size(moves)[1], moves.v_offset + bandwidth + 1)
        moves[i, 1] = Int(ins)
    end
    for j = 2:min(size(moves)[2], moves.h_offset + bandwidth + 1)
        moves[1, j] = Int(del)
    end

    for j = 2:size(result)[2]
        start, stop = row_range(result, j)
        start = max(start, 2)
        for i = start:stop
            x = update(result, i, j, s[i-1], t[j-1],
                       log_p[i-1], log_ins, log_del,
                       allow_codon_indels,
                       log_codon_ins, log_codon_del)
            result[i, j] = x[1]
            moves[i, j] = x[2]
        end
    end
    return backtrace(moves, t, s)
end


"""
F[i, j] is the log probability of aligning s[1:i-1] to t[1:j-1].

"""
function forward_codon(t::AbstractString, s::AbstractString,
                       log_p::Vector{Float64},
                       log_ins::Float64, log_del::Float64, bandwidth::Int,
                       allow_codon_indels::Bool,
                       log_codon_ins::Float64, log_codon_del::Float64)
    result = prepare_array(t, s, log_p, log_ins, log_del, bandwidth)
    for j = 2:size(result)[2]
        start, stop = row_range(result, j)
        start = max(start, 2)
        for i = start:stop
            x = update(result, i, j, s[i-1], t[j-1],
                       log_p[i-1], log_ins, log_del,
                       allow_codon_indels,
                       log_codon_ins, log_codon_del)
            result[i, j] = x[1]
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
                  use_ref, codon_moves::Bool,
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
    if use_ref && codon_moves
        # TODO: make these optional, and test whether they are actually necessary
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


function only_codon_gaps(s::AbstractString)
    cur_gap_len = 0
    for i in 1:length(s)
        if s[i] == '-'
            cur_gap_len += 1
        else
            if cur_gap_len % 3 != 0
                return false
            end
            cur_gap_len = 0
        end
    end
    return cur_gap_len % 3 == 0
end


function is_inframe(allow_stutters::Bool,
                    template::AbstractString,
                    reference::AbstractString,
                    ref_log_p::Vector{Float64},
                    log_ref_ins::Float64, log_ref_del::Float64,
                    bandwidth::Int,
                    allow_codon_indels::Bool,
                    log_codon_ins::Float64, log_codon_del::Float64)
    if allow_stutters
        return length(template) % 3 == 0
    end
    alignment = align(template, reference, ref_log_p,
                      log_ref_ins, log_ref_del,
                      bandwidth, allow_codon_indels,
                      log_codon_ins, log_codon_del)
    t_aln = alignment[1]
    r_aln = alignment[2]
    return only_codon_gaps(t_aln) && only_codon_gaps(r_aln)
end


function quiver2(template::AbstractString,
                 sequences::Vector{ASCIIString},
                 log_ps::Vector{Vector{Float64}},
                 log_ins::Float64, log_del::Float64;
                 reference::AbstractString="",
                 codon_moves::Bool=true,
                 log_codon_ins::Float64=-9.0,
                 log_codon_del::Float64=-9.0,
                 log_mismatch::Float64=-3.0,
                 log_ref_ins::Float64=-9.0,
                 log_ref_del::Float64=-9.0,
                 allow_stutters::Bool=false,
                 cooling_rate::Float64=100.0,
                 bandwidth::Int=10, min_dist::Int=9,
                 batch::Int=10, do_full::Bool=false,
                 max_iters::Int=100, verbose::Int=0)
    if bandwidth < 0
        error("bandwidth cannot be negative: $bandwidth")
    end
    if log_ins >= 0 ||
        log_del >= 0 ||
        log_codon_ins >= 0 ||
        log_codon_del >= 0 ||
        log_mismatch >= 0 ||
        log_ref_ins >= 0 ||
        log_ref_del >= 0
        error("penalties must be < 0")
    end

    if cooling_rate <= 1
        error("cooling rate must be > 1")
    end
    if max_iters < 1
        error("invalid max iters: $max_iters")
    end

    if verbose > 1
        println(STDERR, "computing initial alignments")
    end

    min_log_ref_ins = typemin(Float64)
    min_log_ref_del = typemin(Float64)

    if batch < 0 || batch > length(sequences)
        batch = length(sequences)
    end
    seqs = sequences
    lps = log_ps
    batch = min(batch, length(sequences))
    if batch < length(sequences)
        indices = rand(1:length(sequences), batch)
        seqs = sequences[indices]
        lps = log_ps[indices]
    end

    As = [forward(template, s, p, log_ins, log_del, bandwidth)
          for (s, p) in zip(seqs, lps)]
    Bs = [backward(template, s, p, log_ins, log_del, bandwidth)
          for (s, p) in zip(seqs, lps)]
    current_score = sum([A[end, end] for A in As])

    # will need these later, but do not compute them on the first iteration
    A_t = BandedArray(Float64, (1, 1), 1)
    B_t = BandedArray(Float64, (1, 1), 1)

    use_ref = length(reference) > 0
    enable_ref = false  # ignore reference until first convergence
    reference_log_p = log_mismatch * ones(length(reference))

    current_template = template
    if verbose > 1
        println(STDERR, "initial score: $current_score")
    end
    converged = false
    n_single_mutations = Int[]
    n_codon_mutations = Int[]
    consensus_lengths = Int[]
    n_iterations = 0
    n_full_iterations = 0
    n_ref_iterations = 0
    for i in 1:max_iters
        n_iterations += 1
        if batch == length(sequences)
            n_full_iterations += 1
        end
        if enable_ref
            n_ref_iterations += 1
        end
        old_template = current_template
        old_score = current_score
        if verbose > 1
            println(STDERR, "iteration $i")
        end

        candidates = getcands(current_template, current_score,
                              enable_ref, codon_moves,
                              reference, reference_log_p, A_t, B_t,
                              seqs, lps, As, Bs,
                              log_ins, log_del, bandwidth,
                              log_ref_ins, log_ref_del,
                              log_codon_ins, log_codon_del)
        recompute_As = true
        if length(candidates) == 0
            push!(n_single_mutations, 0)
            push!(n_codon_mutations, 0)
            if !use_ref || is_inframe(allow_stutters,
                                      current_template, reference, reference_log_p,
                                      log_ref_ins, log_ref_del, bandwidth,
                                      true, log_codon_ins, log_codon_del)
                if batch < length(sequences) && do_full
                    if verbose > 0
                        println(STDERR, "converged. switching off batch mode and continuing.")
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
            elseif !enable_ref
                if verbose > 0
                    println(STDERR, "no candidates found. enabling reference.")
                end
                enable_ref = true
            elseif log_ref_ins > min_log_ref_ins || log_ref_del > min_log_ref_del
                if verbose > 1
                    println(STDERR, "no candidates found. increasing indel penalties.")
                end
                log_ref_ins = max(log_ref_ins * cooling_rate, min_log_ref_ins)
                log_ref_del = max(log_ref_del * cooling_rate, min_log_ref_del)
            else
                if verbose > 0
                    println(STDERR, "no candidates found, and penalties cannot change. giving up.\n")
                end
                break
            end
        else
            if verbose > 1
                print(STDERR, "  found $(length(candidates)) candidate mutations.\n")
            end
            chosen_cands = choose_candidates(candidates, min_dist)
            if verbose > 1
                print(STDERR, "  filtered to $(length(chosen_cands)) candidate mutations\n")
            end
            current_template = apply_mutations(old_template,
                                               Mutation[c.mutation
                                                        for c in chosen_cands])
            As = [forward(current_template, s, p, log_ins, log_del, bandwidth)
                  for (s, p) in zip(seqs, lps)]
            temp_score = sum([A[end, end] for A in As])
            if enable_ref
                A_t = forward_codon(current_template, reference, reference_log_p,
                                    log_ref_ins, log_ref_del, bandwidth,
                                    true, log_codon_ins, log_codon_del)
                temp_score += A_t[end, end]
            end
            # detect if a single mutation is better
            # note: this may not always be correct, because score_mutation() is not exact
            if temp_score < chosen_cands[1].score
                if verbose > 1
                    print(STDERR, "  rejecting multiple candidates in favor of best\n")
                end
                chosen_cands = CandMutation[chosen_cands[1]]
                current_template = apply_mutations(old_template,
                                                   Mutation[c.mutation
                                                            for c in chosen_cands])
            else
                # no need to recompute unless batch changes
                recompute_As = false
            end
            n_codon = length(filter(c -> (typeof(c.mutation) == CodonInsertion || typeof(c.mutation) == CodonDeletion),
                                    chosen_cands))
            push!(n_codon_mutations, n_codon)
            push!(n_single_mutations, length(chosen_cands) - n_codon)
        end
        push!(consensus_lengths, length(current_template))
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
            if enable_ref
                A_t = forward_codon(current_template, reference, reference_log_p, log_ref_ins, log_ref_del, bandwidth,
                                    true, log_codon_ins, log_codon_del)
            end
        end
        Bs = [backward(current_template, s, p, log_ins, log_del, bandwidth)
              for (s, p) in zip(seqs, lps)]
        if enable_ref
            B_t = backward_codon(current_template, reference, reference_log_p, log_ref_ins, log_ref_del, bandwidth,
                                 true, log_codon_ins, log_codon_del)
        end
        current_score = sum([A[end, end] for A in As])
        if enable_ref
            current_score += A_t[end, end]
        end
        if verbose > 1
            print(STDERR, "  score: $current_score\n")
        end
    end
    if verbose > 0
        print(STDERR, "done. converged: $converged\n")
    end
    exceeded = n_iterations >= max_iters

    info = Dict("converged" => converged,
                "n_iterations" => n_iterations,
                "n_full_iterations" => n_full_iterations,
                "n_reference_iterations" => n_ref_iterations,
                "exceeded_max_iterations" => exceeded,
                "n_single_mutations" => n_single_mutations,
                "n_codon_mutations" => n_codon_mutations,
                "consensus_lengths" => consensus_lengths,
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
                                        reference::DNASequence=DNASequence(""),
                                        kwargs...)
    new_reference = convert(AbstractString, reference)
    new_template = convert(AbstractString, template)
    new_sequences = ASCIIString[convert(AbstractString, s) for s in sequences]
    result, info = quiver2(new_template, new_sequences, log_ps,
                           log_ins, log_del;
                           reference=new_reference,
                           kwargs...)
    return DNASequence(result), info
end

end
