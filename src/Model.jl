module Model

using Bio.Seq
using Iterators

using Quiver2.BandedArrays
using Quiver2.Mutations
using Quiver2.Util

export quiver2, ErrorModel

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
# scoring stage:
#   - do not change consensus
#   - use reference
#   - propose subsitutions and indels

@enum Stage initial_stage=1 frame_correction_stage=2 refinement_stage=3 scoring_stage=4

immutable ErrorModel
    mismatch::Float64
    insertion::Float64
    deletion::Float64
    codon_insertion::Float64
    codon_deletion::Float64

    function ErrorModel(mismatch::Float64,
                        insertion::Float64,
                        deletion::Float64,
                        codon_insertion::Float64,
                        codon_deletion::Float64)
        args = Float64[mismatch, insertion, deletion, codon_insertion, codon_deletion]
        m, i, d, ci, cd = log10(normalize(args))
        return new(m, i, d, ci, cd)
    end
end

const default_ref_errors = ErrorModel(10.0, 0.1, 0.1, 1.0, 1.0)

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

@enum DPMove dp_none=0 dp_match=1 dp_ins=2 dp_del=3 dp_codon_ins=4 dp_codon_del=5

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


function update_helper(A::BandedArray{Float64}, i::Int, j::Int,
                       move::DPMove, move_score::Float64,
                       final_score::Float64, final_move::DPMove)
    offset = offsets[Int(move)]
    i = i + offset[1]
    j = j + offset[2]
    if inband(A, i, j)
        score = A[i, j] + move_score
        if score > final_score
            return score, move
        end
    end
    return final_score, final_move
end


function move_scores(t_base::Char,
                     s_base::Char,
                     seq_i::Int,
                     log_p::Vector{Float64},
                     errors::ErrorModel)
    cur_log_p = log_p[max(seq_i, 1)]
    next_log_p = log_p[min(seq_i + 1, length(log_p))]
    match_score = s_base == t_base ? inv_log10(cur_log_p) : cur_log_p + errors.mismatch
    ins_score = cur_log_p + errors.insertion
    del_score = max(cur_log_p, next_log_p) + errors.deletion
    return match_score, ins_score, del_score
end


function codon_move_scores(t_base::Char,
                           s_base::Char,
                           seq_i::Int,
                           log_p::Vector{Float64},
                           errors::ErrorModel)
    log_p_1 = log_p[max(seq_i, 1)]
    log_p_2 = log_p[min(seq_i + 1, length(log_p))]
    log_p_3 = log_p[min(seq_i + 2, length(log_p))]
    max_p = maximum([log_p_1, log_p_2, log_p_3])
    codon_ins_score = max_p + errors.codon_insertion
    codon_del_score = max(log_p_1, log_p_2) + errors.codon_deletion
    return codon_ins_score, codon_del_score
end


function update(A::BandedArray{Float64}, i::Int, j::Int,
                s_base::Char, t_base::Char,
                log_p::Vector{Float64},
                errors::ErrorModel)
    result = (typemin(Float64), dp_none)
    match_score, ins_score, del_score = move_scores(t_base, s_base, i-1, log_p, errors)
    result = update_helper(A, i, j, dp_match, match_score, result...)
    result = update_helper(A, i, j, dp_ins, ins_score, result...)
    result = update_helper(A, i, j, dp_del, del_score, result...)
    codon_ins_score, codon_del_score = codon_move_scores(t_base, s_base, i-1,
                                                         log_p, errors)
    if errors.codon_insertion > -Inf && i > 3
        result = update_helper(A, i, j, dp_codon_ins,
                               codon_ins_score, result...)
    end
    if errors.codon_deletion > -Inf && j > 3
        result = update_helper(A, i, j, dp_codon_del,
                               codon_del_score, result...)
    end
    if result[1] == -Inf
        error("new score is invalid")
    end
    if result[2] == dp_none
        error("failed to find a move")
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
        if move == dp_match
            push!(aligned_t, t[tj])
            push!(aligned_s, s[si])
        elseif move == dp_ins
            push!(aligned_t, '-')
            push!(aligned_s, s[si])
        elseif move == dp_del
            push!(aligned_t, t[tj])
            push!(aligned_s, '-')
        elseif move == dp_codon_ins
            append!(aligned_t, ['-', '-', '-'])
            append!(aligned_s, [s[si], s[si-1], s[si-2]])
        elseif move == dp_codon_del
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
                       errors::ErrorModel,
                       bandwidth::Int)
    # FIXME: consider codon moves in initialization
    if length(s) != length(log_p)
        error("sequence length does not match quality score length")
    end
    result = BandedArray(Float64, (length(s) + 1, length(t) + 1), bandwidth)
    if errors.insertion == -Inf || errors.deletion == -Inf
        error("indel probabilities cannot be 0.0")
    end
    for i = 2:min(size(result)[1], result.v_offset + bandwidth + 1)
        ins_score = log_p[i-1] + errors.insertion
        result[i, 1] = result[i-1, 1] + ins_score
    end
    del_score = log_p[1] + + errors.deletion
    for j = 2:min(size(result)[2], result.h_offset + bandwidth + 1)
        result[1, j] = result[1, j-1] + del_score
    end
    return result
end


"""Does some work as forward_codon, but also keeps track of moves

Does backtracing to find best alignment.

"""
function forward_moves(t::AbstractString, s::AbstractString,
                       log_p::Vector{Float64},
                       errors::ErrorModel,
                       bandwidth::Int)
    # FIXME: code duplication with forward_codon(). This is done in a
    # seperate function to keep return type stable and avoid
    # allocating the `moves` array unnecessarily.
    result = prepare_array(t, s, log_p, errors, bandwidth)
    moves = BandedArray(Int, result.shape, bandwidth)
    moves[1, 1] = Int(dp_none)
    for i = 2:min(size(moves)[1], moves.v_offset + bandwidth + 1)
        moves[i, 1] = Int(dp_ins)
    end
    for j = 2:min(size(moves)[2], moves.h_offset + bandwidth + 1)
        moves[1, j] = Int(dp_del)
    end

    for j = 2:size(result)[2]
        start, stop = row_range(result, j)
        start = max(start, 2)
        for i = start:stop
            x = update(result, i, j, s[i-1], t[j-1], log_p, errors)
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
                 log_p::Vector{Float64}, errors::ErrorModel,
                 bandwidth::Int)
    result = prepare_array(t, s, log_p, errors, bandwidth)
    for j = 2:size(result)[2]
        start, stop = row_range(result, j)
        start = max(start, 2)
        for i = start:stop
            x = update(result, i, j, s[i-1], t[j-1],
                       log_p, errors)
            result[i, j] = x[1]
        end
    end
    return result::BandedArray{Float64}
end


"""
B[i, j] is the log probability of aligning s[i:end] to t[j:end].

"""
function backward(t::AbstractString, s::AbstractString,
                  log_p::Vector{Float64},
                  errors::ErrorModel,
                  bandwidth::Int)
    s = reverse(s)
    log_p = flipdim(log_p, 1)
    t = reverse(t)
    result = forward(t, s, log_p, errors, bandwidth)
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

function seq_score_mutation(mutation::Union{Deletion,CodonDeletion},
                            A::BandedArray{Float64}, B::BandedArray{Float64},
                            template::AbstractString,
                            seq::AbstractString, log_p::Vector{Float64},
                            errors::ErrorModel, codon_moves::Bool)
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

function seq_score_mutation(mutation::Union{Insertion,Substitution},
                            A::BandedArray{Float64}, B::BandedArray{Float64},
                            template::AbstractString,
                            seq::AbstractString, log_p::Vector{Float64},
                            errors::ErrorModel, codon_moves::Bool)
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
        sbase = seq_i > 0 ? seq[seq_i] : 'X'
        match_score, ins_score, del_score = move_scores(mutation.base, sbase, seq_i, log_p, errors)
        if i > 1
            # insertion
            scores[1] = max(scores[1], scores[2] + ins_score)
        end
        if prev_start < real_i <= (prev_stop + 1)
            # (mis)match
            scores[1] = max(scores[1], A[real_i - 1, ajprev] + match_score)
        end
        if prev_start <= real_i <= prev_stop
            # deletion
            scores[1] = max(scores[1], A[real_i, ajprev] + del_score)
        end
        if codon_moves
            codon_ins_score, codon_del_score = codon_move_scores(mutation.base, sbase, seq_i, log_p, errors)
            if i > 3
                # insertion
                scores[1] = max(scores[1], scores[4] + codon_ins_score)
            end
            if ajprev3 > 0 && (prev3_start <= real_i <= prev3_stop) && inband(A, real_i, ajprev3)
                # deletion
                scores[1] = max(scores[1], A[real_i, ajprev-2] + codon_del_score)
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
                        base, seq, log_p, errors, codon_moves)
    col = fill(typemin(Float64), row_stop - row_start + 1)
    offset = 1 - (row_start - prev_start)
    for i in 1:length(col)
        real_i = i + row_start - 1
        seq_i = real_i - 1  # position in the sequence
        prev_i = i - offset + 1 # position in `prev` matching i
        sbase = seq_i > 0 ? seq[seq_i] : 'X'
        match_score, ins_score, del_score = move_scores(base, sbase, seq_i, log_p, errors)
        if i > 1
            # insertion
            col[i] = max(col[i], col[i-1] + ins_score)
        end
        if prev_start < real_i <= (prev_stop + 1)
            # (mis)match
            col[i] = max(col[i], prev[prev_i-1] + match_score)
        end
        if prev_start <= real_i <= prev_stop
            # deletion
            col[i] = max(col[i], prev[prev_i] + del_score)
        end
        if codon_moves
            codon_ins_score, codon_del_score = codon_move_scores(base, sbase, seq_i, log_p, errors)
            if i > 3
                # insertion
                col[i] = max(col[i], col[i-3] + codon_ins_score)
            end
            if has_prev3 && prev3_start <= real_i <= prev3_stop
                # deletion
                prev3_i = real_i - prev3_start + 1
                col[i] = max(col[i], prev3[prev3_i] + codon_del_score)
            end
        end
    end
    return col
end

function seq_score_mutation(mutation::CodonInsertion,
                            A::BandedArray{Float64}, B::BandedArray{Float64},
                            template::AbstractString,
                            seq::AbstractString, log_p::Vector{Float64},
                            errors::ErrorModel, codon_moves::Bool)
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
                                 errors, codon_moves)

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
                                 errors, codon_moves)

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
                                 errors, codon_moves)

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
                        errors::ErrorModel,
                        use_ref::Bool,
                        reference::ASCIIString,
                        ref_log_p::Vector{Float64},
                        ref_errors::ErrorModel)
    score = 0.0
    for si in 1:length(sequences)
        score += seq_score_mutation(m, state.As[si], state.Bs[si], state.template,
                                    sequences[si], log_ps[si], errors, false)
    end
    if use_ref
        score += seq_score_mutation(m, state.A_t, state.B_t, state.template,
                                    reference, ref_log_p, ref_errors, true)
    end
    return score
end


function candstask(stage::Stage,
                   template::AbstractString)
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
            stage == frame_correction_stage ||
            stage == scoring_stage
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
            # codon deletions
            for j in 1:(len-2)
                produce(CodonDeletion(j))
            end
            # codon insertions
            for codon in product("ACGT", "ACGT", "ACGT")
                for j in 0:len
                    produce(CodonInsertion(j, codon))
                end
            end
        end
    end
    Task(_it)
end


function getcands(state::State,
                  sequences::Vector{ASCIIString},
                  log_ps::Vector{Vector{Float64}},
                  errors::ErrorModel,
                  reference::ASCIIString,
                  ref_log_p::Vector{Float64},
                  ref_errors::ErrorModel)
    candidates = CandMutation[]
    use_ref = (state.stage == frame_correction_stage)
    for m in candstask(state.stage, state.template)
        score = score_mutation(m, state,
                               sequences, log_ps, errors,
                               use_ref,
                               reference, ref_log_p, ref_errors)
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


function align(t::AbstractString, s::AbstractString,
               log_p::Vector{Float64},
               errors::ErrorModel,
               bandwidth::Int)
    moves = forward_moves(t, s, log_p, errors, bandwidth)
    return backtrace(t, s, moves)
end


function no_single_indels(template::AbstractString,
                          reference::AbstractString,
                          ref_log_p::Vector{Float64},
                          ref_errors::ErrorModel,
                          bandwidth::Int)
    has_right_length = length(template) % 3 == 0
    t_aln, r_aln = align(template, reference, ref_log_p, ref_errors, bandwidth)
    result = only_codon_gaps(t_aln) && only_codon_gaps(r_aln)
    if result && !has_right_length
        error("template length is not a multiple of three")
    end
    return result
end


function initial_state(template, seqs, lps, errors, bandwidth)
    As = [forward(template, s, p, errors, bandwidth)
          for (s, p) in zip(seqs, lps)]
    Bs = [backward(template, s, p, errors, bandwidth)
          for (s, p) in zip(seqs, lps)]
    score = sum([A[end, end] for A in As])

    A_t = BandedArray(Float64, (1, 1), 1)
    B_t = BandedArray(Float64, (1, 1), 1)

    return State(score, template, A_t, B_t, As, Bs, initial_stage, false)
end


function recompute!(state::State, seqs::Vector{ASCIIString},
                    lps::Vector{Vector{Float64}},
                    errors::ErrorModel,
                    reference::AbstractString,
                    ref_log_p::Vector{Float64},
                    ref_errors::ErrorModel,
                    bandwidth::Int, recompute_As::Bool, recompute_Bs::Bool)
    if recompute_As
        state.As = [forward(state.template, s, p, errors, bandwidth)
                    for (s, p) in zip(seqs, lps)]
        if ((state.stage == frame_correction_stage ||
             state.stage == scoring_stage) &&
            length(reference) > 0)
            state.A_t = forward(state.template, reference, ref_log_p,
                                errors, bandwidth)
        end
    end
    if recompute_Bs
        state.Bs = [backward(state.template, s, p, errors, bandwidth)
                    for (s, p) in zip(seqs, lps)]
        if ((state.stage == frame_correction_stage ||
             state.stage == scoring_stage) &&
            length(reference) > 0)
            state.B_t = backward(state.template, reference, ref_log_p,
                                 ref_errors, bandwidth)
        end
    end
    state.score = sum([A[end, end] for A in state.As])
    if ((state.stage == frame_correction_stage ||
         state.stage == scoring_stage) &&
        length(reference) > 0)
        state.score += state.A_t[end, end]
    end
end


"""convert per-mutation log differences to a per-base error rate"""
function normalize_log_differences(position_scores, insertion_scores, state_score)
    # per-base insertion score is mean of neighboring insertions
    position_exp = exp10(position_scores)
    position_probs = broadcast(/, position_exp, sum(position_exp, 2))
    ins_exp = exp10(insertion_scores)
    ins_probs = broadcast(/, ins_exp, exp10(state_score) + sum(ins_exp, 2))
    return position_probs, ins_probs
end


function estimate_probs(state::State,
                        sequences::Vector{ASCIIString},
                        log_ps::Vector{Vector{Float64}},
                        errors::ErrorModel,
                        reference::AbstractString,
                        ref_log_p::Vector{Float64},
                        ref_errors::ErrorModel)
    # `position_scores[i]` gives the following log probabilities
    # for `template[i]`: [A, C, G, T, -]
    position_scores = zeros(length(state.template), 5) + state.score
    # `insertion_scores[i]` gives the following log probabilities for an
    # insertion before `template[i]` of [A, C, G, T]
    insertion_scores = zeros(length(state.template) + 1, 4)

    # TODO: should we modify penalties before using reference?
    # - do not penalize mismatches
    # - use max indel penalty

    use_ref = (length(reference) > 0)
    for m in candstask(scoring_stage, state.template)
        score = score_mutation(m, state,
                               sequences, log_ps, errors,
                               use_ref,
                               reference, ref_log_p, ref_errors)
        if typeof(m) == Substitution
            position_scores[m.pos, baseints[m.base]] = score
        elseif typeof(m) == Deletion
            position_scores[m.pos, baseints['-']] = score
        elseif typeof(m) == Insertion
            insertion_scores[m.pos + 1, baseints[m.base]] = score
        end
    end
    max_score = max(maximum(position_scores), maximum(insertion_scores))
    position_scores -= max_score
    insertion_scores -= max_score
    if maximum(position_scores) > 0.0
        error("position scores cannot be positive")
    end
    if maximum(insertion_scores) > 0.0
        error("insertion scores cannot be positive")
    end
    return normalize_log_differences(position_scores, insertion_scores, state.score - max_score)
end


function estimate_point_probs(position_probs, insertion_probs)
    no_point_error_prob = maximum(position_probs, 2)
    # multiple by 0.5 to avoid double counting.
    # TODO: is this the right way to do this?
    no_ins_error_prob = 1.0 - 0.5 * sum(insertion_probs, 2)
    result = 1.0 - broadcast(*, no_point_error_prob,
                             no_ins_error_prob[1:end-1],
                             no_ins_error_prob[2:end])
    return reshape(result, length(result))
end


function estimate_indel_probs(position_probs, insertion_probs)
    no_del_prob = 1.0 - position_probs[:, end]
    no_ins_error_prob = 1.0 - 0.5 * sum(insertion_probs, 2)
    result = 1.0 - broadcast(*, no_del_prob,
                             no_ins_error_prob[1:end-1],
                             no_ins_error_prob[2:end])
    return reshape(result, length(result))
end


function quiver2(template::AbstractString,
                 sequences::Vector{ASCIIString},
                 phreds::Vector{Vector{Int8}},
                 errors::ErrorModel;
                 reference::AbstractString="",
                 ref_log_p::Float64=0.0,
                 ref_errors::ErrorModel=default_ref_errors,
                 cooling_rate::Float64=2.0,
                 bandwidth::Int=10, min_dist::Int=9,
                 batch::Int=10, batch_threshold::Float64=0.05,
                 max_iters::Int=100, verbose::Int=0)
    if errors.codon_insertion > -Inf || errors.codon_deletion > -Inf
        error("error model allows codon indels")
    end
    if errors.insertion == -Inf || errors.deletion == -Inf
        error("indel probabilities cannot be 0.0")
    end

    if bandwidth < 0
        error("bandwidth cannot be negative: $bandwidth")
    end

    if length(reference) > 0
        if (ref_log_p == -Inf || ref_log_p >= 0.0)
            error("ref_log_p=$ref_log_p but should be less than 0.0")
        end
        if ref_errors.insertion == -Inf || ref_errors.deletion == -Inf
            error("ref indel probabilities cannot be 0.0")
        end
    end

    ref_log_p_vec = fill(ref_log_p, length(reference))

    if cooling_rate <= 1
        error("cooling rate must be > 1")
    end
    if max_iters < 1
        error("invalid max iters: $max_iters")
    end

    if any([minimum(p) < 0 for p in phreds])
        error("phred score cannot be negative")
    end
    log_ps = phred_to_log_p(phreds)
    if any([minimum(p) == -Inf for p in log_ps])
        error("a log error probability is negative infinity")
    end

    if batch < 0 || batch > length(sequences)
        batch = length(sequences)
    end
    base_batch = batch
    seqs = sequences
    lps = log_ps
    if batch < length(sequences)
        indices = rand(1:length(sequences), batch)
        seqs = sequences[indices]
        lps = log_ps[indices]
    end

    if verbose > 1
        println(STDERR, "computing initial alignments")
    end
    state = initial_state(template, seqs, lps, errors, bandwidth)
    empty_ref = length(reference) == 0

    if verbose > 1
        println(STDERR, "initial score: $(state.score)")
    end
    n_mutations = Vector{Int}[]
    consensus_lengths = Int[length(template)]
    consensus_noref = ""
    consensus_ref = ""

    stage_iterations = zeros(Int, Int(typemax(Stage)))
    for i in 1:max_iters
        stage_iterations[Int(state.stage)] += 1
        old_template = state.template
        old_score = state.score
        if verbose > 1
            println(STDERR, "iteration $i : $(state.stage)")
        end

        candidates = getcands(state, seqs, lps, errors,
                              reference, ref_log_p_vec, ref_errors)

        recompute_As = true
        if length(candidates) == 0
            if verbose > 1
                println(STDERR, "  no candidates found")
            end
            push!(n_mutations, zeros(Int, Int(typemax(DPMove))))
            if state.stage == initial_stage
                consensus_noref = state.template
                if empty_ref
                    state.converged = true
                    break
                end
                state.stage = frame_correction_stage
            elseif state.stage == frame_correction_stage
                if no_single_indels(state.template, reference, ref_log_p_vec, ref_errors, bandwidth)
                    consensus_ref = state.template
                    state.stage = refinement_stage
                else
                    # TODO: decrease indel probability
                    if verbose > 1
                        println(STDERR, "  alignment had single indels but scores already minimized.")
                    end
                    state.stage = refinement_stage
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
            recompute!(state, seqs, lps, errors,
                       reference, ref_log_p_vec, ref_errors,
                       bandwidth, true, false)
            # detect if a single mutation is better
            # note: this may not always be correct, because score_mutation() is not exact
            if length(chosen_cands) > 1 && (state.score < chosen_cands[1].score
                                            || isapprox(state.score, chosen_cands[1].score))
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
        recompute!(state, seqs, lps, errors,
                   reference, ref_log_p_vec, ref_errors,
                   bandwidth, recompute_As, true)
        if verbose > 1
            println(STDERR, "  score: $(state.score)")
        end
        if ((state.score - old_score) / old_score > batch_threshold &&
            batch < length(sequences))
            batch = min(batch + base_batch, length(sequences))
            if verbose > 1
                println(STDERR, "  increased batch size to $batch")
            end
            indices = rand(1:length(sequences), batch)
            seqs = sequences[indices]
            lps = log_ps[indices]
            recompute!(state, seqs, lps, errors,
                       reference, ref_log_p_vec, ref_errors,
                       bandwidth, true, true)
            if verbose > 1
                println(STDERR, "  new score: $(state.score)")
            end

        end
    end
    state.stage = scoring_stage
    if verbose > 0
        println(STDERR, "done. converged: $(state.converged)")
    end
    push!(consensus_lengths, length(state.template))
    exceeded = sum(stage_iterations) >= max_iters

    info = Dict("converged" => state.converged,
                "stage_iterations" => stage_iterations,
                "exceeded_max_iterations" => exceeded,
                "ref_error_model" => ref_errors,
                "consensus_noref" => consensus_noref,
                "consensus_ref" => consensus_ref,
                "n_mutations" => transpose(hcat(n_mutations...)),
                "consensus_lengths" => consensus_lengths,
                )

    # FIXME: recomputing for all sequences may be costly
    recompute!(state, sequences, log_ps, errors,
               reference, ref_log_p_vec, ref_errors,
               bandwidth, true, true)
    base_probs, ins_probs = estimate_probs(state, sequences, log_ps, errors,
                                           reference, ref_log_p_vec, ref_errors)
    return state.template, base_probs, ins_probs, info
end

"""
Alternate quiver2() using BioJulia types.

"""
function quiver2(template::DNASequence,
                 sequences::Vector{DNASequence},
                 phreds::Vector{Vector{Int8}},
                 errors::ErrorModel;
                 reference::DNASequence=DNASequence(""),
                 kwargs...)
    new_reference = convert(ASCIIString, reference)
    new_template = convert(ASCIIString, template)
    new_sequences = ASCIIString[convert(ASCIIString, s) for s in sequences]
    (result, base_probs,
     insertion_probs, info) = quiver2(new_template, new_sequences, phreds, errors;
                                      reference=new_reference,
                                      kwargs...)
    info["consensus_noref"] = DNASequence(info["consensus_noref"])
    info["consensus_ref"] = DNASequence(info["consensus_ref"])
    return (DNASequence(result), base_probs, insertion_probs, info)
end

end
