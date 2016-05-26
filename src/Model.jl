# TODO: documentation

# TODO: annealing schedule for log_ins and log_del penalties for
# aligning template to reference

# FIXME: swap log_ins and log_del penalties for aligning template to
# reference

module Model

using Bio.Seq
using DataStructures

using Quiver2.Sample
using Quiver2.BandedArrays
using Quiver2.Mutations

export quiver2, Penalties

# initial_stage:
#   - do not use reference.
#   - propose mismatches and single indels.
# frame_correction_stage:
#   - use reference.
#   - propose mismatches and single indels.
#   - increase reference single indel penalties.
# refinement_stage:
#   - use reference.
#   - propose mismatches and codon indels.
#   - maximize reference single indel penalties.

@enum Stage initial_stage=1 frame_correction_stage=2 refinement_stage=3

type State
    score::Float64
    template::AbstractString
    A_t::BandedArray{Float64}
    B_t::BandedArray{Float64}
    As::Vector{BandedArray{Float64}}
    Bs::Vector{BandedArray{Float64}}
    stage::Stage
    converged::Bool
end


immutable Penalties
    ins::Float64
    del::Float64
    codon_ins::Float64
    codon_del::Float64
end


function Penalties(ins, del, codon_ins, codon_del)
    if ins >= 0 || del >= 0 || codon_ins >= 0 || codon_del >= 0
        error("penalties must be < 0")
    end
    return Penalties()
end


function Penalties(ins, del)
    return Penalties(ins, del, typemin(Float64), typemin(Float64))
end



const default_penalties = Penalties(-2.0, -2.0)
const default_ref_penalties = Penalties(-9.0, -9.0, -9.0, -9.0)


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
                log_p::Float64, penalties::Penalties,
                allow_codon_indels::Bool)
    result = (typemin(Float64), match)
    log_match = 0.0
    if t_base != 'N'
        log_match = (s_base == t_base ? 0.0 : log_p)
    end
    result = update_helper(A, i, j, match, log_match, result...)
    result = update_helper(A, i, j, ins, penalties.ins, result...)
    result = update_helper(A, i, j, del, penalties.del, result...)
    if allow_codon_indels
        if i > 3
            result = update_helper(A, i, j, codon_ins, penalties.codon_ins, result...)
        end
        if j > 3
            result = update_helper(A, i, j, codon_del, penalties.codon_del, result...)
        end
    end
    return result
end


function backtrace(t::AbstractString, s::AbstractString, moves::BandedArray{Int})
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
                       penalties::Penalties,
                       bandwidth::Int)
    # FIXME: consider codon moves in initialization
    if length(s) != length(log_p)
        error("sequence length does not match quality score length")
    end
    result = BandedArray(Float64, (length(s) + 1, length(t) + 1), bandwidth)
    for i = 2:min(size(result)[1], result.v_offset + bandwidth + 1)
        result[i, 1] = penalties.ins * (i - 1)
    end
    for j = 2:min(size(result)[2], result.h_offset + bandwidth + 1)
        result[1, j] = penalties.del * (j - 1)
    end
    return result
end


""" Does some work as forward_codon, but also keeps track of moves
does backtracing to find best alignment.

"""
function forward_moves(t::AbstractString, s::AbstractString,
                       log_p::Vector{Float64},
                       penalties::Penalties,
                       bandwidth::Int,
                       allow_codon_indels::Bool)
    # FIXME: code duplication with forward_codon(). This is done in a
    # seperate function to keep return type stable and avoid
    # allocating the `moves` array unnecessarily.
    result = prepare_array(t, s, log_p, penalties, bandwidth)
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
                       log_p[i-1], penalties,
                       allow_codon_indels)
            result[i, j] = x[1]
            moves[i, j] = x[2]
        end
    end
    return moves
end


"""
F[i, j] is the log probability of aligning s[1:i-1] to t[1:j-1].

"""
function forward(t::AbstractString, s::AbstractString,
                 log_p::Vector{Float64},
                 penalties::Penalties,
                 bandwidth::Int;
                 allow_codon_indels::Bool=false)
    result = prepare_array(t, s, log_p, penalties, bandwidth)
    for j = 2:size(result)[2]
        start, stop = row_range(result, j)
        start = max(start, 2)
        for i = start:stop
            x = update(result, i, j, s[i-1], t[j-1],
                       log_p[i-1], penalties,
                       allow_codon_indels)
            result[i, j] = x[1]
        end
    end
    return result::BandedArray{Float64}
end


"""
B[i, j] is the log probability of aligning s[i:end] to t[j:end].

"""
function backward(t::AbstractString, s::AbstractString, log_p::Vector{Float64},
                  penalties::Penalties,
                  bandwidth::Int;
                  allow_codon_indels::Bool=false)
    s = reverse(s)
    log_p = flipdim(log_p, 1)
    t = reverse(t)
    result = forward(t, s, log_p, penalties, bandwidth; allow_codon_indels=allow_codon_indels)
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

function score_mutation(mutation::Union{Deletion,CodonDeletion},
                        A::BandedArray{Float64}, B::BandedArray{Float64},
                        args...)
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
                        A::BandedArray{Float64}, B::BandedArray{Float64},
                        template::AbstractString,
                        seq::AbstractString, log_p::Vector{Float64},
                        penalties::Penalties,
                        allow_codon_indels::Bool)
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
            scores[1] = max(scores[1], scores[2] + penalties.ins)
        end
        if prev_start < real_i <= (prev_stop + 1)
            # (mis)match
            scores[1] = max(scores[1], A[real_i - 1, ajprev] + (mutation.base == seq[seq_i] ? 0.0 : log_p[seq_i]))
        end
        if prev_start <= real_i <= prev_stop
            # deletion
            scores[1] = max(scores[1], A[real_i, ajprev] + penalties.del)
        end
        if allow_codon_indels
            if i > 3
                # insertion
                scores[1] = max(scores[1], scores[4] + penalties.codon_ins)
            end
            if ajprev3 > 0 && (prev3_start <= real_i <= prev3_stop) && inband(A, real_i, ajprev3)
                # deletion
                scores[1] = max(scores[1], A[real_i, ajprev-2] + penalties.codon_del)
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
                        penalties,
                        allow_codon_indels)
    col = typemin(Float64) * ones(row_stop - row_start + 1)
    offset = 1 - (row_start - prev_start)
    for i in 1:length(col)
        real_i = i + row_start - 1
        seq_i = real_i - 1  # position in the sequence
        prev_i = i - offset + 1 # position in `prev` matching i
        if i > 1
            # insertion
            col[i] = max(col[i], col[i-1] + penalties.ins)
        end
        if prev_start < real_i <= (prev_stop + 1)
            # (mis)match
            col[i] = max(col[i], prev[prev_i-1] + (base == seq[seq_i] ? 0.0 : log_p[seq_i]))
        end
        if prev_start <= real_i <= prev_stop
            # deletion
            col[i] = max(col[i], prev[prev_i] + penalties.del)
        end
        if allow_codon_indels
            if i > 3
                # insertion
                col[i] = max(col[i], col[i-3] + penalties.codon_ins)
            end
            if has_prev3 && prev3_start <= real_i <= prev3_stop
                # deletion
                prev3_i = real_i - prev3_start + 1
                col[i] = max(col[i], prev3[prev3_i] + penalties.codon_del)
            end
        end
    end
    return col
end

function score_mutation(mutation::CodonInsertion,
                        A::BandedArray{Float64}, B::BandedArray{Float64},
                        template::AbstractString,
                        seq::AbstractString, log_p::Vector{Float64},
                        penalties::Penalties,
                        allow_codon_indels::Bool)
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

    # FIXME: code duplication
    Acols[:, 1] = compute_subcol(row_start, row_stop,
                                 prev, prev_start, prev_stop,
                                 prev3, prev3_start, prev3_stop,
                                 has_prev3,
                                 mutation.bases[1], seq, log_p,
                                 penalties,
                                 allow_codon_indels)

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
                                 penalties,
                                 allow_codon_indels)

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
                                 penalties,
                                 allow_codon_indels)

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
            score += score_mutation($m, state.As[si], state.Bs[si], state.template,
                                    sequences[si], log_ps[si],
                                    penalties, false)
        end
        if state.stage > initial_stage
            score += score_mutation($m, state.A_t, state.B_t, state.template,
                                    reference, reference_log_p,
                                    ref_penalties, true)

        end
        if score > state.score
            push!(candidates, CandMutation($m, score))
        end
    end
end

function getcands(state::State,
                  sequences::Vector{ASCIIString},
                  log_ps::Vector{Vector{Float64}},
                  penalties::Penalties,
                  reference::AbstractString,
                  reference_log_p::Vector{Float64},
                  ref_penalties::Penalties,
                  bandwidth::Int)
    len = length(state.template)
    candidates = CandMutation[]
    # insertion at beginning
    if state.stage != refinement_stage
        for base in "ACGT"
            mi = Insertion(0, base)
            @process_mutation mi
        end
    end
    for j in 1:len
        for base in "ACGT"
            if state.stage != refinement_stage
                mi = Insertion(j, base)
                @process_mutation mi
            end
            if state.template[j] != base
                ms = Substitution(j, base)
                @process_mutation ms
            end
        end
        if state.stage != refinement_stage
            md = Deletion(j)
            @process_mutation md
        end
    end
    if state.stage == refinement_stage
        for j in 1:len
            mci = CodonInsertion(j, ('N', 'N', 'N'))
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


function is_inframe(check_alignment::Bool,
                    template::AbstractString,
                    reference::AbstractString,
                    ref_log_p::Vector{Float64},
                    penalties::Penalties,
                    bandwidth::Int)
    has_right_length = length(template) % 3 == 0
    if !check_alignment
        return has_right_length
    end
    moves = forward_moves(template, reference,
                          ref_log_p, penalties,
                          bandwidth, true)
    t_aln, r_aln = backtrace(template, reference, moves)
    result = only_codon_gaps(t_aln) && only_codon_gaps(r_aln)
    if result && !has_right_length
        error("template length is not a multiple of three")
    end
    return result
end


function initial_state(template, seqs, lps, penalties, bandwidth)
    As = [forward(template, s, p, penalties, bandwidth)
          for (s, p) in zip(seqs, lps)]
    Bs = [backward(template, s, p, penalties, bandwidth)
          for (s, p) in zip(seqs, lps)]
    score = sum([A[end, end] for A in As])

    # will need these later, but do not compute them on the first iteration
    A_t = BandedArray(Float64, (1, 1), 1)
    B_t = BandedArray(Float64, (1, 1), 1)

    return State(score, template, A_t, B_t, As, Bs, initial_stage, false)
end


function recompute!(state::State, seqs::Vector{ASCIIString},
                    lps::Vector{Vector{Float64}},
                    penalties::Penalties,
                    reference::AbstractString,
                    reference_log_p,
                    ref_penalties::Penalties,
                    bandwidth::Int, recompute_As::Bool, recompute_Bs::Bool)
    if recompute_As
        state.As = [forward(state.template, s, p, penalties, bandwidth)
                    for (s, p) in zip(seqs, lps)]
        if state.stage > initial_stage
            state.A_t = forward(state.template, reference, reference_log_p, ref_penalties,
                                bandwidth, allow_codon_indels=true)
        end
    end
    if recompute_Bs
        state.Bs = [backward(state.template, s, p, penalties, bandwidth)
                    for (s, p) in zip(seqs, lps)]
        if state.stage > initial_stage
            state.B_t = backward(state.template, reference, reference_log_p, ref_penalties, bandwidth,
                                 allow_codon_indels=true)
        end
    end
    state.score = sum([A[end, end] for A in state.As])
    if state.stage > initial_stage
        state.score += state.A_t[end, end]
    end
end


function best_base(counter::DataStructures.Accumulator{Char, Int})
    if '-' in keys(counter)
        pop!('-', counter)
    end
    c = rbase()
    if length(counter) > 0
        c, n = maximum(counter)
    end
    return c
end


function interleave(xs, ys)
    if length(xs) != length(ys) + 2
        error("`xs` must have two more elements than `ys`")
    end
    result = []
    for i in 1:length(ys)
        push!(result, xs[i])
        push!(result, ys[i])
    end
    push!(result, xs[end])
    return result
end


function replace_ns!(state::State, seqs::Vector{ASCIIString},
                    log_ps::Vector{Vector{Float64}}, penalties::Penalties,
                    bandwidth::Int)
    positions = []
    for i in 1:length(state.template)
        if state.template[i] == 'N'
            push!(positions, i)
        end
    end
    counters = [counter(Char) for i in positions]
    for (s, p) in zip(seqs, lps)
        moves = forward_moves(state.state.template, s, p, penalties, bandwidth, false)
        t_aln, s_aln = backtrace(state.template, s, moves)
        posn = 0
        for i in 1:length(t_aln)
            if t_aln[i] == 'N'
                posn += 1
                push!(counters[posn], s_aln[i])
            end
        end
    end
    bases = [best_base(c) for c in counters]
    parts = split(state.template, 'N')
    state.template = join(interleave(parts, bases))
end


function quiver2(template::AbstractString,
                 sequences::Vector{ASCIIString},
                 log_ps::Vector{Vector{Float64}};
                 reference::AbstractString="",
                 penalties::Penalties=default_penalties,
                 ref_penalties::Penalties=default_ref_penalties,
                 ref_mismatch::Float64=-3.0,
                 check_alignment::Bool=true,
                 cooling_rate::Float64=100.0,
                 bandwidth::Int=10, min_dist::Int=9,
                 batch::Int=10, do_full::Bool=false,
                 max_iters::Int=100, verbose::Int=0)
    if bandwidth < 0
        error("bandwidth cannot be negative: $bandwidth")
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

    if ref_mismatch >= 0
        error("reference mismatch penalty must be less than 0")
    end

    if batch < 0 || batch > length(sequences)
        batch = length(sequences)
    end
    seqs = sequences
    lps = log_ps
    if batch < length(sequences)
        indices = rand(1:length(sequences), batch)
        seqs = sequences[indices]
        lps = log_ps[indices]
    end

    state = initial_state(template, seqs, lps, penalties, bandwidth)
    empty_ref = length(reference) == 0
    reference_log_p = ref_mismatch * ones(length(reference))

    min_ref_ins = typemin(Float64)
    min_ref_del = typemin(Float64)

    if verbose > 1
        println(STDERR, "initial score: $(state.score)")
    end
    n_mutations = Vector{Int}[]
    consensus_lengths = Int[]

    stage_iterations = zeros(Int, Int(typemax(Stage)))
    for i in 1:max_iters
        stage_iterations[Int(state.stage)] += 1
        old_template = state.template
        old_score = state.score
        if verbose > 1
            println(STDERR, "iteration $i")
        end

        candidates = getcands(state, seqs, lps, penalties,
                              reference, reference_log_p, ref_penalties,
                              bandwidth)

        recompute_As = true
        if length(candidates) == 0
            push!(n_mutations, zeros(Int, Int(typemax(DPMove))))
            inframe = false
            if !empty_ref
                inframe = is_inframe(check_alignment, state.template,
                                     reference, reference_log_p,
                                     ref_penalties, bandwidth)
            end
            if state.stage == initial_stage
                if empty_ref
                    state.converged = true
                    break
                end
                if inframe
                    state.stage = refinement_stage
                    ref_penalties = Penalties(min_ref_ins,
                                              min_ref_del,
                                              ref_penalties.codon_ins,
                                              ref_penalties.codon_del)
                else
                    state.stage = frame_correction_stage
                end
            elseif state.stage == frame_correction_stage
                if inframe
                    state.stage = refinement_stage
                    ref_penalties = Penalties(min_ref_ins,
                                              min_ref_del,
                                              ref_penalties.codon_ins,
                                              ref_penalties.codon_del)
                elseif ref_penalties.ins > min_ref_ins || ref_penalties.del > min_ref_del
                    ref_penalties = Penalties(max(ref_penalties.ins * cooling_rate, min_ref_ins),
                                              max(ref_penalties.del * cooling_rate, min_ref_del),
                                              ref_penalties.codon_ins,
                                              ref_penalties.codon_del)
                else
                    # cannot increase penalties any more. give up.
                    break
                end
            elseif state.stage == refinement_stage
                if inframe
                    state.converged = true
                end
                break
            else
                error("unknown stage: $(state.stage)")
            end
        else
            if verbose > 1
                print(STDERR, "  found $(length(candidates)) candidate mutations.\n")
            end
            chosen_cands = choose_candidates(candidates, min_dist)
            if verbose > 1
                print(STDERR, "  filtered to $(length(chosen_cands)) candidate mutations\n")
            end
            state.template = apply_mutations(old_template,
                                             Mutation[c.mutation
                                                      for c in chosen_cands])
            recompute!(state, seqs, lps,
                       penalties, reference, reference_log_p, ref_penalties,
                       bandwidth, true, false)
            # detect if a single mutation is better
            # note: this may not always be correct, because score_mutation() is not exact
            if state.score < chosen_cands[1].score
                if verbose > 1
                    print(STDERR, "  rejecting multiple candidates in favor of best\n")
                end
                chosen_cands = CandMutation[chosen_cands[1]]
                state.template = apply_mutations(old_template,
                                                 Mutation[c.mutation
                                                          for c in chosen_cands])
            else
                # no need to recompute unless batch changes
                recompute_As = false
            end
            mutation_counts = [length(filter(c -> (typeof(c.mutation) == t),
                                             chosen_cands))
                               for t in [Substitution, Insertion, Deletion, CodonInsertion, CodonDeletion]]
            push!(n_mutations, mutation_counts)
        end
        push!(consensus_lengths, length(state.template))
        if batch < length(sequences)
            indices = rand(1:length(sequences), batch)
            seqs = sequences[indices]
            lps = log_ps[indices]
            recompute_As = true
        end
        recompute!(state, seqs, lps, penalties,
                   reference, reference_log_p, ref_penalties,
                   bandwidth, recompute_As, true)
        if 'N' in state.template
            if state.stage != refinement_stage
                error("'N' only allowed in template during refinement stage")
            end
            # replace 'N' with best base from alignments and recompute everything
            replace_ns!(state, seqs, lps, penalties, bandwidth)
            recompute!(state, seqs, lps, penalties,
                       reference, reference_log_p, ref_penalties,
                       bandwidth, true, true)
        end
        if verbose > 1
            print(STDERR, "  score: $(state.score)\n")
        end
    end
    if verbose > 0
        print(STDERR, "done. converged: $state.converged\n")
    end
    push!(consensus_lengths, length(state.template))
    exceeded = sum(stage_iterations) >= max_iters

    info = Dict("converged" => state.converged,
                "stage_iterations" => stage_iterations,
                "exceeded_max_iterations" => exceeded,
                "n_mutations" => n_mutations,
                "consensus_lengths" => consensus_lengths,
                )
    return state.template, info
end

"""
Alternate quiver2() using BioJulia types.

"""
function quiver2{T<:NucleotideSequence}(template::DNASequence,
                                        sequences::Vector{T},
                                        log_ps::Vector{Vector{Float64}};
                                        reference::DNASequence=DNASequence(""),
                                        kwargs...)
    new_reference = convert(AbstractString, reference)
    new_template = convert(AbstractString, template)
    new_sequences = ASCIIString[convert(AbstractString, s) for s in sequences]
    result, info = quiver2(new_template, new_sequences, log_ps;
                           reference=new_reference,
                           kwargs...)
    return DNASequence(result), info
end

end
