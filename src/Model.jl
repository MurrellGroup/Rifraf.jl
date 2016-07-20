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
#   - propose all mutations
#   - increase reference single indel penalties.
# refinement stage:
#   - do not use reference
#   - propose subsitutions only

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
    ins_multiplier::Float64
    del_multiplier::Float64
    codon_ins::Float64
    codon_del::Float64
end


function Penalties(ins, del, codon_ins, codon_del)
    if ins < 1 || del < 1
        error("multipliers must be >= 1")
    end
    if codon_ins >= 0 || codon_del >= 0
        error("penalties must be < 0")
    end
    return Penalties()
end

const default_penalties = Penalties(1.0, 1.0, -9.0, -9.0)

@enum DPMove match=1 ins=2 del=3 codon_ins=4 codon_del=5

const offsets = ([-1, -1],  # sub
                 [-1, 0],   # insertion
                 [0, -1],   # deletion
                 [-3, 0],   # codon insertion
                 [0, -3])   # codon deletion


const baseints = Dict('A' => 1,
                      'C' => 2,
                      'G' => 3,
                      'T' => 4,
                      '-' => 5,
                      )


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
                log_p::Float64, next_log_p::Float64,
                allow_codon_indels::Bool, penalties::Penalties)
    result = (typemin(Float64), match)
    # TODO: precompute inv_log_p
    match_penalty = (s_base == t_base ? log10(1.0 - exp10(log_p)) : log_p)
    if t_base == 'N'
        match_penalty = 0.0
    end
    ins_penalty = log_p * (allow_codon_indels ? penalties.ins_multiplier : 1.0)
    del_penalty = mean([log_p, next_log_p]) * (allow_codon_indels ? penalties.del_multiplier : 1.0)
    result = update_helper(A, i, j, match, match_penalty, result...)
    result = update_helper(A, i, j, ins, ins_penalty, result...)
    result = update_helper(A, i, j, del, del_penalty, result...)
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
                       bandwidth::Int)
    # FIXME: consider codon moves in initialization
    if length(s) != length(log_p)
        error("sequence length does not match quality score length")
    end
    result = BandedArray(Float64, (length(s) + 1, length(t) + 1), bandwidth)
    for i = 2:min(size(result)[1], result.v_offset + bandwidth + 1)
        result[i, 1] = result[i-1, 1] + log_p[i-1]
    end
    for j = 2:min(size(result)[2], result.h_offset + bandwidth + 1)
        result[1, j] = log_p[1] * (j - 1)
    end
    return result
end


"""Does some work as forward_codon, but also keeps track of moves does
backtracing to find best alignment.

"""
function forward_moves(t::AbstractString, s::AbstractString,
                       log_p::Vector{Float64},
                       bandwidth::Int,
                       allow_codon_indels::Bool, penalties::Penalties)
    # FIXME: code duplication with forward_codon(). This is done in a
    # seperate function to keep return type stable and avoid
    # allocating the `moves` array unnecessarily.
    result = prepare_array(t, s, log_p, bandwidth)
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
                       log_p[i-1], log_p[min(i, length(log_p))],
                       allow_codon_indels, penalties)
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
                 bandwidth::Int;
                 allow_codon_indels::Bool=false,
                 penalties::Penalties=default_penalties)
    result = prepare_array(t, s, log_p, bandwidth)
    for j = 2:size(result)[2]
        start, stop = row_range(result, j)
        start = max(start, 2)
        for i = start:stop
            x = update(result, i, j, s[i-1], t[j-1],
                       log_p[i-1], log_p[min(i, length(log_p))],
                       allow_codon_indels, penalties)
            result[i, j] = x[1]
        end
    end
    return result::BandedArray{Float64}
end


"""
B[i, j] is the log probability of aligning s[i:end] to t[j:end].

"""
function backward(t::AbstractString, s::AbstractString, log_p::Vector{Float64},
                  bandwidth::Int;
                  allow_codon_indels::Bool=false,
                  penalties::Penalties=default_penalties)
    s = reverse(s)
    log_p = flipdim(log_p, 1)
    t = reverse(t)
    result = forward(t, s, log_p, bandwidth;
                     allow_codon_indels=allow_codon_indels, penalties=penalties)
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
                        allow_codon_indels::Bool, penalties::Penalties)
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
        penalty = log_p[max(seq_i, 1)]
        if i > 1
            # insertion
            ins_penalty = penalty * (allow_codon_indels ? penalties.ins_multiplier : 1.0)
            scores[1] = max(scores[1], scores[2] + ins_penalty)
        end
        if prev_start < real_i <= (prev_stop + 1)
            # (mis)match
            mm_penalty = (mutation.base == seq[seq_i] ? log10(1.0 - exp10(penalty)) : penalty)
            scores[1] = max(scores[1], A[real_i - 1, ajprev] + mm_penalty)
        end
        if prev_start <= real_i <= prev_stop
            # deletion
            next_log_p = log_p[min(seq_i + 1, length(log_p))]
            del_penalty = mean([penalty, next_log_p]) * (allow_codon_indels ? penalties.del_multiplier : 1.0)
            scores[1] = max(scores[1], A[real_i, ajprev] + del_penalty)
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
                        allow_codon_indels, penalties)
    col = fill(typemin(Float64), row_stop - row_start + 1)
    offset = 1 - (row_start - prev_start)
    for i in 1:length(col)
        real_i = i + row_start - 1
        seq_i = real_i - 1  # position in the sequence
        prev_i = i - offset + 1 # position in `prev` matching i
        penalty = log_p[max(seq_i, 1)]
        if i > 1
            # insertion
            ins_penalty = penalty * (allow_codon_indels ? penalties.ins_multiplier : 1.0)
            col[i] = max(col[i], col[i-1] + ins_penalty)
        end
        if prev_start < real_i <= (prev_stop + 1)
            # (mis)match
            mm_penalty = (base == seq[seq_i] ? log10(1.0 - exp10(penalty)) : penalty)
            col[i] = max(col[i], prev[prev_i-1] + mm_penalty)
        end
        if prev_start <= real_i <= prev_stop
            # deletion
            next_log_p = log_p[min(seq_i + 1, length(log_p))]
            del_penalty = mean([penalty, next_log_p]) * (allow_codon_indels ? penalties.del_multiplier : 1.0)
            col[i] = max(col[i], prev[prev_i] + del_penalty)
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
                        allow_codon_indels::Bool, penalties::Penalties)
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
                                 allow_codon_indels, penalties)

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
                                 allow_codon_indels, penalties)

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
                                 allow_codon_indels, penalties)

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


function score_mutation(m::Mutation,
                        state::State,
                        sequences::Vector{ASCIIString},
                        log_ps::Vector{Vector{Float64}},
                        use_ref::Bool,
                        reference::ASCIIString,
                        reference_log_p::Vector{Float64},
                        penalties::Penalties)
    score = 0.0
    for si in 1:length(sequences)
        score += score_mutation(m, state.As[si], state.Bs[si], state.template,
                                sequences[si], log_ps[si],
                                false, penalties)
    end
    if use_ref
        score += score_mutation(m, state.A_t, state.B_t, state.template,
                                reference, reference_log_p, true, penalties)
    end
    return score
end


function candstask(stage::Stage,
                   template::AbstractString,
                   sequences::Vector{ASCIIString},
                   log_ps::Vector{Vector{Float64}})
    len = length(template)
    function _it()
        # substitutions
        for j in 1:len
            for base in "ACGT"
                if template[j] != base
                    produce(Substitution(j, base))
                end
            end
        end
        if stage == initial_stage ||
            stage == frame_correction_stage
            # single indels
            for base in "ACGT"
                produce(Insertion(0, base))
            end
            for j in 1:len
                for base in "ACGT"
                    produce(Insertion(j, base))
                end
                produce(Deletion(j))
            end
        end
        if stage == frame_correction_stage
            # codon indels
            produce(CodonInsertion(0, ('N', 'N', 'N')))
            for j in 1:len
                produce(CodonInsertion(j, ('N', 'N', 'N')))
                if j < len - 1
                    produce(CodonDeletion(j))
                end
            end
        end
    end
    Task(_it)
end


function getcands(state::State,
                  sequences::Vector{ASCIIString},
                  log_ps::Vector{Vector{Float64}},
                  reference::ASCIIString,
                  reference_log_p::Vector{Float64},
                  penalties::Penalties)
    candidates = CandMutation[]
    for m in candstask(state.stage, state.template,
                       sequences, log_ps)
        score = score_mutation(m, state, sequences, log_ps,
                               state.stage == frame_correction_stage,
                               reference, reference_log_p, penalties)
        if score > state.score && !isapprox(score, state.score)
            push!(candidates, CandMutation(m, score))
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


function align(t, s, log_p, bandwidth;
               allow_codon_indels::Bool=false,
               penalties::Penalties=default_penalties)
    moves = forward_moves(t, s, log_p, bandwidth, allow_codon_indels, penalties)
    return backtrace(t, s, moves)
end


function no_single_indels(template::AbstractString,
                          reference::AbstractString,
                          ref_log_p::Vector{Float64},
                          penalties::Penalties,
                          bandwidth::Int)
    has_right_length = length(template) % 3 == 0
    t_aln, r_aln = align(template, reference, ref_log_p, bandwidth,
                         allow_codon_indels=true, penalties=penalties)
    result = only_codon_gaps(t_aln) && only_codon_gaps(r_aln)
    if result && !has_right_length
        error("template length is not a multiple of three")
    end
    return result
end


function initial_state(template, seqs, lps, bandwidth)
    As = [forward(template, s, p, bandwidth)
          for (s, p) in zip(seqs, lps)]
    Bs = [backward(template, s, p, bandwidth)
          for (s, p) in zip(seqs, lps)]
    score = sum([A[end, end] for A in As])

    A_t = BandedArray(Float64, (1, 1), 1)
    B_t = BandedArray(Float64, (1, 1), 1)

    return State(score, template, A_t, B_t, As, Bs, initial_stage, false)
end


function recompute!(state::State, seqs::Vector{ASCIIString},
                    lps::Vector{Vector{Float64}},
                    reference::AbstractString,
                    reference_log_p,
                    penalties::Penalties,
                    bandwidth::Int, recompute_As::Bool, recompute_Bs::Bool)
    if recompute_As
        state.As = [forward(state.template, s, p, bandwidth)
                    for (s, p) in zip(seqs, lps)]
        if state.stage == frame_correction_stage
            state.A_t = forward(state.template, reference, reference_log_p,
                                bandwidth, allow_codon_indels=true, penalties=penalties)
        end
    end
    if recompute_Bs
        state.Bs = [backward(state.template, s, p, bandwidth)
                    for (s, p) in zip(seqs, lps)]
        if state.stage == frame_correction_stage
            state.B_t = backward(state.template, reference, reference_log_p, bandwidth,
                                 allow_codon_indels=true, penalties=penalties)
        end
    end
    state.score = sum([A[end, end] for A in state.As])
    if state.stage == frame_correction_stage
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
                    log_ps::Vector{Vector{Float64}},
                    bandwidth::Int)
    positions = []
    for i in 1:length(state.template)
        if state.template[i] == 'N'
            push!(positions, i)
        end
    end
    counters = [counter(Char) for i in positions]
    for (s, p) in zip(seqs, lps)
        t_aln, s_aln = align(state.state.template, s, p, bandwidth)
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


"""convert per-mutation log differences to a per-base error rate"""
function normalize_log_differences(position_scores, insertion_scores)
    # per-base insertion score is mean of neighboring insertions
    position_exp = exp10(position_scores)
    ins_exp = exp10(insertion_scores)
    position_probs = broadcast(/, position_exp, sum(position_exp, 2))
    ins_probs = broadcast(/, ins_exp, 1.0 + sum(ins_exp, 2))
    return position_probs, ins_probs
end


function estimate_probs(state::State,
                        sequences::Vector{ASCIIString},
                        log_ps::Vector{Vector{Float64}})
    # `position_scores[i]` gives the following log probabilities
    # for `template[i]`: [A, C, G, T, -]
    position_scores = zeros(length(state.template), 5)
    # `insertion_scores[i]` gives the following log probabilities for an
    # insertion before `template[i]` of [A, C, G, T]
    insertion_scores = zeros(length(state.template) + 1, 4)
    for m in candstask(initial_stage, state.template, sequences, log_ps)
        score = score_mutation(m, state, sequences, log_ps,
                               false, "", Float64[], default_penalties)
        score -= state.score
        if typeof(m) == Substitution
            position_scores[m.pos, baseints[m.base]] = score
        elseif typeof(m) == Deletion
            position_scores[m.pos, baseints['-']] = score
        elseif typeof(m) == Insertion
            insertion_scores[m.pos + 1, baseints[m.base]] = score
        end
    end
    return normalize_log_differences(position_scores, insertion_scores)
end


function estimate_point_probs(position_probs, insertion_probs)
    no_point_error_prob = maximum(position_probs, 2)
    # multiple by 0.5 to avoid double counting.
    # TODO: is this the right way to do this?
    no_ins_error_prob = 1.0 - 0.5 * sum(insertion_probs, 2)
    return 1.0 - broadcast(*, no_point_error_prob,
                           no_ins_error_prob[1:end-1],
                           no_ins_error_prob[2:end])
end


function quiver2(template::AbstractString,
                 sequences::Vector{ASCIIString},
                 phreds::Vector{Vector{Int8}};
                 reference::AbstractString="",
                 penalties::Penalties=default_penalties,
                 ref_mismatch::Float64=-3.0,
                 cooling_rate::Float64=100.0,
                 bandwidth::Int=10, min_dist::Int=9,
                 batch::Int=10,
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

    if verbose > 0
        println(STDERR, "computing initial alignments")
    end

    if ref_mismatch >= 0
        error("reference mismatch penalty must be less than 0")
    end

    log_ps = phreds / (-10.0)

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

    state = initial_state(template, seqs, lps, bandwidth)
    empty_ref = length(reference) == 0
    reference_log_p = fill(ref_mismatch, length(reference))

    min_ref_ins = typemin(Float64)
    min_ref_del = typemin(Float64)

    if verbose > 0
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
            println(STDERR, "iteration $i : $(state.stage)")
        end

        candidates = getcands(state, seqs, lps,
                              reference, reference_log_p, penalties)

        recompute_As = true
        if length(candidates) == 0
            if verbose > 1
                println(STDERR, "  no candidates found")
            end
            push!(n_mutations, zeros(Int, Int(typemax(DPMove))))
            if state.stage == initial_stage
                if empty_ref
                    state.converged = true
                    break
                end
                state.stage = frame_correction_stage
            elseif state.stage == frame_correction_stage
                if no_single_indels(state.template,
                                    reference, reference_log_p,
                                    penalties, bandwidth)
                    state.stage = refinement_stage
                else
                    penalties = Penalties(penalties.ins_multiplier * cooling_rate,
                                          penalties.del_multiplier * cooling_rate,
                                          penalties.codon_ins,
                                          penalties.codon_del)
                    if verbose > 1
                        println(STDERR, "  alignment to reference had single indels. increasing penalty.")
                    end
                    # FIXME: detect maximum penalty
                end
            elseif state.stage == refinement_stage
                state.converged = true
                break
            else
                error("unknown stage: $(state.stage)")
            end
        else
            if verbose > 1
                println(STDERR, "  found $(length(candidates)) candidate mutations.")
            end
            chosen_cands = choose_candidates(candidates, min_dist)
            if verbose > 1
                println(STDERR, "  filtered to $(length(chosen_cands)) candidate mutations")
            end
            state.template = apply_mutations(old_template,
                                             Mutation[c.mutation
                                                      for c in chosen_cands])
            recompute!(state, seqs, lps,
                       reference, reference_log_p, penalties,
                       bandwidth, true, false)
            # detect if a single mutation is better
            # note: this may not always be correct, because score_mutation() is not exact
            if state.score < chosen_cands[1].score || isapprox(state.score, chosen_cands[1].score)
                if verbose > 1
                    println(STDERR, "  rejecting multiple candidates in favor of best")
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
        recompute!(state, seqs, lps,
                   reference, reference_log_p, penalties,
                   bandwidth, recompute_As, true)
        if 'N' in state.template
            if state.stage != frame_correction_stage
                error("'N' not allowed in template during this stage")
            end
            # replace 'N' with best base from alignments and recompute everything
            replace_ns!(state, seqs, lps, bandwidth)
            recompute!(state, seqs, lps,
                       reference, reference_log_p, penalties,
                       bandwidth, true, true)
        end
        if verbose > 1
            println(STDERR, "  score: $(state.score)")
        end
    end
    if verbose > 0
        println(STDERR, "done. converged: $(state.converged)")
    end
    push!(consensus_lengths, length(state.template))
    exceeded = sum(stage_iterations) >= max_iters

    info = Dict("converged" => state.converged,
                "stage_iterations" => stage_iterations,
                "exceeded_max_iterations" => exceeded,
                "n_mutations" => transpose(hcat(n_mutations...)),
                "consensus_lengths" => consensus_lengths,
                )

    # FIXME: recomputing for all sequences may be costly
    recompute!(state, sequences, log_ps,
               reference, reference_log_p, penalties,
               bandwidth, true, true)
    base_probs, insertion_probs = estimate_probs(state, sequences, log_ps)
    point_probs = estimate_point_probs(base_probs, insertion_probs)
    return state.template, base_probs, insertion_probs, point_probs, info
end

"""
Alternate quiver2() using BioJulia types.

"""
function quiver2(template::DNASequence,
                 sequences::Vector{DNASequence},
                 phreds::Vector{Vector{Int8}};
                 reference::DNASequence=DNASequence(""),
                 kwargs...)
    new_reference = convert(ASCIIString, reference)
    new_template = convert(ASCIIString, template)
    new_sequences = ASCIIString[convert(ASCIIString, s) for s in sequences]
    (result, base_probs, insertion_probs,
     point_probs, info) = quiver2(new_template, new_sequences, phreds;
                                  reference=new_reference,
                                  kwargs...)
    return (DNASequence(result), base_probs,
            insertion_probs, point_probs, info)
end

end
