module Model

using Bio.Seq
using Iterators

using Quiver2.BandedArrays
using Quiver2.Proposals
using Quiver2.Util

export quiver2, ErrorModel, Scores, normalize, modified_emissions

# initial_stage:
#   - do not use reference.
#   - propose all proposals
# frame_correction_stage:
#   - use reference. (allow codon moves in reference alignment)
#   - propose all proposals
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
end

function ErrorModel(mismatch::Float64,
                    insertion::Float64,
                    deletion::Float64)
    return ErrorModel(mismatch, insertion, deletion, 0.0, 0.0)
end

function normalize(errors::ErrorModel)
    args = Float64[errors.mismatch,
                   errors.insertion,
                   errors.deletion,
                   errors.codon_insertion,
                   errors.codon_deletion]
    m, i, d, ci, cd = Util.normalize(args)
    return ErrorModel(m, i, d, ci, cd)
end

immutable Scores
    mismatch::Float64
    insertion::Float64
    deletion::Float64
    codon_insertion::Float64
    codon_deletion::Float64
end

function Scores(errors::ErrorModel)
    args = Float64[errors.mismatch,
                   errors.insertion,
                   errors.deletion,
                   errors.codon_insertion,
                   errors.codon_deletion]
    m, i, d, ci, cd = log10(Util.normalize(args))
    return Scores(m, i, d, ci, cd)
end

function modified_emissions(scores::Scores;
                            mismatch_emission::Float64=(1.0 / 3.0),
                            insertion_emission::Float64=0.25)
    return Scores(scores.mismatch + log10(mismatch_emission),
                  scores.insertion + log10(insertion_emission),
                  scores.deletion,
                  scores.codon_insertion + log10(insertion_emission^3),
                  scores.codon_deletion)
end

# default to invalid scores, to force user to set them
const default_ref_scores = Scores(0.0, 0.0, 0.0, 0.0, 0.0)

# just to avoid magical constants in code
const codon_length = 3

type State
    score::Float64
    template::AbstractString
    A_t::BandedArray{Float64}
    B_t::BandedArray{Float64}
    As::Vector{BandedArray{Float64}}
    Amoves::Vector{BandedArray{Int}}
    Bs::Vector{BandedArray{Float64}}
    stage::Stage
    converged::Bool
    bandwidth::Int
end

@enum DPMove dp_none=0 dp_match=1 dp_ins=2 dp_del=3 dp_codon_ins=4 dp_codon_del=5

const offsets = ([-1, -1],  # sub
                 [-1, 0],   # insertion
                 [0, -1],   # deletion
                 [-3, 0],   # codon insertion
                 [0, -3])   # codon deletion

const bases = "ACGT-"
const baseints = Dict('A' => 1,
                      'C' => 2,
                      'G' => 3,
                      'T' => 4,
                      '-' => 5,
                      )

const empty_array = Array(Float64, (0, 0))

function move_scores(t_base::Char,
                     s_base::Char,
                     seq_i::Int,
                     log_p::Vector{Float64},
                     scores::Scores)
    cur_log_p = log_p[max(seq_i, 1)]
    next_log_p = log_p[min(seq_i + 1, length(log_p))]
    match_score = s_base == t_base ? inv_log10(cur_log_p) : cur_log_p + scores.mismatch
    ins_score = cur_log_p + scores.insertion
    del_score = max(cur_log_p, next_log_p) + scores.deletion
    return match_score, ins_score, del_score
end

function codon_move_scores(seq_i::Int,
                           log_p::Vector{Float64},
                           scores::Scores)
    # we're moving INTO seq_i. so need previous three
    start = max(1, seq_i-2)
    stop = min(seq_i, length(log_p))
    max_p = start <= stop ? maximum(log_p[start:stop]) : -Inf
    codon_ins_score =  max_p + scores.codon_insertion
    cur_log_p = log_p[max(seq_i, 1)]
    next_log_p = log_p[min(seq_i + 1, length(log_p))]
    codon_del_score = max(cur_log_p, next_log_p) + scores.codon_deletion
    return codon_ins_score, codon_del_score
end

function update_helper(newcols::Array{Float64, 2},
                       A::BandedArray{Float64},
                       i::Int, j::Int, acol::Int,
                       move::DPMove, move_score::Float64,
                       final_score::Float64, final_move::DPMove)
    offset = offsets[Int(move)]
    prev_i = i + offset[1]
    prev_j = j + offset[2]

    rangecol = min(prev_j, size(A)[2])
    if inband(A, prev_i, rangecol)
        score = -Inf
        if acol < 1 || prev_j <= acol
            score = A[prev_i, prev_j] + move_score
        else
            score = newcols[prev_i, prev_j - acol] + move_score
        end
        if score > final_score
            return score, move
        end
    end
    return final_score, final_move
end

function update(A::BandedArray{Float64},
                i::Int, j::Int,
                s_base::Char, t_base::Char,
                log_p::Vector{Float64},
                scores::Scores;
                newcols::Array{Float64, 2}=empty_array,
                acol=-1)
    result = (-Inf, dp_none)
    match_score, ins_score, del_score = move_scores(t_base, s_base, i-1, log_p, scores)
    result = update_helper(newcols, A, i, j, acol, dp_match, match_score, result...)
    result = update_helper(newcols, A, i, j, acol, dp_ins, ins_score, result...)
    result = update_helper(newcols, A, i, j, acol, dp_del, del_score, result...)
    codon_ins_score, codon_del_score = codon_move_scores(i-1, log_p, scores)
    if scores.codon_insertion > -Inf && i > codon_length
        result = update_helper(newcols, A, i, j, acol, dp_codon_ins,
                               codon_ins_score, result...)
    end
    if scores.codon_deletion > -Inf && j > codon_length
        result = update_helper(newcols, A, i, j, acol, dp_codon_del,
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


"""Does backtracing to find best alignment.

"""
function forward_moves(t::AbstractString, s::AbstractString,
                       log_p::Vector{Float64},
                       scores::Scores,
                       bandwidth::Int)
    # FIXME: code duplication with forward_codon(). This is done in a
    # seperate function to keep return type stable and avoid
    # allocating the `moves` array unnecessarily.
    result = BandedArray(Float64, (length(s) + 1, length(t) + 1), bandwidth)
    moves = BandedArray(Int, result.shape, bandwidth)
    moves[1, 1] = Int(dp_none)
    nrows, ncols = size(result)
    for j = 1:ncols
        start, stop = row_range(result, j)
        for i = start:stop
            if i == 1 && j == 1
                continue
            end
            sbase = i > 1 ? s[i-1] : 'X'
            tbase = j > 1 ? t[j-1] : 'X'
            x = update(result, i, j, sbase, tbase,
                       log_p, scores)
            result[i, j] = x[1]
            moves[i, j] = x[2]
        end
    end
    return result, moves
end


"""
F[i, j] is the log probability of aligning s[1:i-1] to t[1:j-1].

"""
function forward(t::AbstractString, s::AbstractString,
                 log_p::Vector{Float64}, scores::Scores,
                 bandwidth::Int)
    result = BandedArray(Float64, (length(s) + 1, length(t) + 1), bandwidth)
    nrows, ncols = size(result)
    for j = 1:ncols
        start, stop = row_range(result, j)
        for i = start:stop
            if i == 1 && j == 1
                continue
            end
            sbase = i > 1 ? s[i-1] : 'X'
            tbase = j > 1 ? t[j-1] : 'X'
            x = update(result, i, j, sbase, tbase,
                       log_p, scores)
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
                  scores::Scores,
                  bandwidth::Int)
    result = forward(reverse(t), reverse(s), reverse(log_p),
                     scores, bandwidth)
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


function get_sub_template(proposal::Proposal, seq::AbstractString,
                          next_posn::Int, n_after::Int)
    # next valid position in sequence after this proposal
    t = typeof(proposal)
    pos = proposal.pos
    prefix = ""
    stop = min(next_posn + n_after - 1, length(seq))
    suffix = seq[next_posn:stop]
    if t in (Substitution, Insertion)
        prefix = string(proposal.base)
    elseif t == CodonInsertion
        prefix = string(proposal.bases...)
    end
    return string(prefix, suffix)
end

function seq_score_deletion(A::BandedArray{Float64}, B::BandedArray{Float64},
                            acol::Int, bcol::Int)
    Acol = sparsecol(A, acol)
    Bcol = sparsecol(B, bcol)
    (amin, amax), (bmin, bmax) = equal_ranges(row_range(A, acol),
                                              row_range(B, bcol))
    asub = sub(Acol, amin:amax)
    bsub = sub(Bcol, bmin:bmax)
    return summax(asub, bsub)
end


const n_proposal_bases = Dict(Substitution => 1,
                              Insertion => 1,
                              Deletion => 0,
                              CodonInsertion => 3,
                              CodonDeletion => 0)

const boffsets = Dict(Substitution => 2,
                      Insertion => 1,
                      Deletion => 2,
                      CodonInsertion => 1,
                      CodonDeletion => 4)

function score_nocodon(proposal::Proposal,
                       A::BandedArray{Float64}, B::BandedArray{Float64},
                       template::AbstractString,
                       seq::AbstractString, log_p::Vector{Float64},
                       scores::Scores,
                       newcols::Array{Float64, 2})
    target_col = proposal.pos + 1
    t = typeof(proposal)
    if t in (Deletion, CodonDeletion)
        # nothing to recompute
        n_to_del = t == Deletion ? 1 : 3
        acol = proposal.pos
        bcol = proposal.pos + n_to_del
        return seq_score_deletion(A, B, acol, bcol)
    end
    # need to compute new columns
    # TODO: reuse an array for `newcols`
    n_new = (t == CodonInsertion? 3 : 1)
    nrows, ncols = size(A)

    # last valid A column
    acol = proposal.pos + (t == Substitution ? 0 : 1)
    # new bases
    sub_template = t == CodonInsertion ? string(proposal.bases...) : string(proposal.base)
    for j in 1:n_new
        amin, amax = row_range(A, min(acol + j, ncols))
        for i in amin:amax
            seq_base = i > 1 ? seq[i-1] : 'X'
            x = update(A, i, acol + j,
                       seq_base, sub_template[j],
                       log_p, scores;
                       newcols=newcols, acol=acol)
            newcols[i, j] = x[1]
        end
    end

    # add up results
    imin, imax = row_range(A, min(acol + n_new, ncols))
    Acol = newcols[amin:amax, n_new]

    bj = proposal.pos + 1
    Bcol = sparsecol(B, bj)
    (amin, amax), (bmin, bmax) = equal_ranges((imin, imax),
                                              row_range(B, bj))
    asub = sub(Acol, amin:amax)
    bsub = sub(Bcol, bmin:bmax)
    score = summax(asub, bsub)
    if score == -Inf
        error("failed to compute a valid score")
    end
    return score
end

function seq_score_proposal(proposal::Proposal,
                            A::BandedArray{Float64}, B::BandedArray{Float64},
                            template::AbstractString,
                            seq::AbstractString, log_p::Vector{Float64},
                            scores::Scores,
                            newcols::Array{Float64, 2})
    codon_moves = (scores.codon_insertion > -Inf ||
                   scores.codon_deletion > -Inf)
    if !codon_moves
        return score_nocodon(proposal, A, B,
                             template, seq, log_p,
                             scores, newcols)
    end
    t = typeof(proposal)
    # last valid column of A
    acol_offset = t in (Insertion, CodonInsertion) ? 0 : -1
    acol = proposal.pos + acol_offset + 1

    # first column of B to use
    first_bcol = acol + boffsets[t]
    # last column of B to use
    last_bcol = first_bcol + codon_length - 1

    if t in (Deletion, CodonDeletion)
        n_del = (t == Deletion ? 1 : codon_length)
        if acol == (size(A)[2] - n_del)
            # suffix deletions do not need recomputation
            return A[end, end - n_del]
        end
    end

    # number of bases changed/inserted
    n_bases = n_proposal_bases[t]
    # number of columns after recomputed columns to also recompute.
    n_after = codon_length

    # if we'd go to or past the last column of B, just recompute the
    # rest of A
    nrows, ncols = size(A)
    just_a = last_bcol >= ncols
    next_posn = proposal.pos + (t == CodonDeletion ? 3 : 1)
    if just_a
        # go to end of template
        n_after = length(template) - next_posn + 1
    end

    if n_bases == 0 && n_after == 0
        error("no new columns need to be recomputed.")
    end

    sub_template = get_sub_template(proposal, template,
                                    next_posn, n_after)

    # TODO: reuse an array for `newcols`
    n_new = n_bases + n_after

    # compute new columns
    for j in 1:n_new
        range_col = min(acol + j, ncols)
        amin, amax = row_range(A, range_col)
        for i in amin:amax
            seq_base = i > 1 ? seq[i-1] : 'X'
            x = update(A, i, acol + j,
                       seq_base, sub_template[j],
                       log_p, scores;
                       newcols=newcols, acol=acol)
            newcols[i, j] = x[1]
        end
    end

    if just_a
        return newcols[nrows, n_new]
    end

    # add up results
    best_score = -Inf
    for j in 1:codon_length
        new_j = n_new - codon_length + j
        imin, imax = row_range(A, min(acol + new_j, ncols))
        Acol = newcols[imin:imax, new_j]
        bj = first_bcol + j - 1
        if bj > size(B)[2]
            error("wrong column")
        else
            Bcol = sparsecol(B, bj)
            (amin, amax), (bmin, bmax) = equal_ranges((imin, imax),
                                                      row_range(B, bj))
            asub = sub(Acol, amin:amax)
            bsub = sub(Bcol, bmin:bmax)
            score = summax(asub, bsub)
        end
        if score > best_score
            best_score = score
        end
    end
    if best_score == -Inf
        error("failed to compute a valid score")
    end
    return best_score
end

function choose_candidates(candidates::Vector{CandProposal}, min_dist::Int)
    final_cands = CandProposal[]
    posns = Set()
    for c in sort(candidates, by=(c) -> c.score, rev=true)
        if any(Bool[(abs(c.proposal.pos - p) < min_dist) for p in posns])
            continue
        end
        union!(posns, affected_positions(c.proposal))
        push!(final_cands, c)
    end
    return final_cands
end


function score_proposal(m::Proposal,
                        state::State,
                        sequences::Vector{ASCIIString},
                        log_ps::Vector{Vector{Float64}},
                        scores::Scores,
                        use_ref::Bool,
                        reference::ASCIIString,
                        ref_log_p::Vector{Float64},
                        ref_scores::Scores,
                        newcols::Array{Float64, 2})
    score = 0.0
    for si in 1:length(sequences)
        score += seq_score_proposal(m, state.As[si], state.Bs[si], state.template,
                                    sequences[si], log_ps[si], scores, newcols)
    end
    if use_ref
        score += seq_score_proposal(m, state.A_t, state.B_t, state.template,
                                    reference, ref_log_p, ref_scores, newcols)
    end
    return score
end


function candstask(stage::Stage,
                   template::AbstractString,
                   sequences::Vector{ASCIIString},
                   log_ps::Vector{Vector{Float64}},
                   Amoves::Vector{BandedArray{Int}},
                   scores::Scores,
                   bandwidth::Int)
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
            for proposal in best_codons(template, sequences, log_ps,
                                        scores, bandwidth,
                                        Amoves=Amoves)
                produce(proposal)
            end
        end
    end
    Task(_it)
end


function getcands(state::State,
                  sequences::Vector{ASCIIString},
                  log_ps::Vector{Vector{Float64}},
                  scores::Scores,
                  reference::ASCIIString,
                  ref_log_p::Vector{Float64},
                  ref_scores::Scores)
    candidates = CandProposal[]
    use_ref = (state.stage == frame_correction_stage)

    maxlen = maximum([length(s) for s in sequences])
    nrows = max(maxlen, length(reference)) + 1
    newcols = zeros(Float64, (nrows, 6))

    for m in candstask(state.stage, state.template,
                       sequences, log_ps, state.Amoves,
                       scores, state.bandwidth)
        score = score_proposal(m, state,
                               sequences, log_ps, scores,
                               use_ref,
                               reference, ref_log_p, ref_scores,
                               newcols)
        if score > state.score && !isapprox(score, state.score)
            push!(candidates, CandProposal(m, score))
        end
    end
    return candidates
end


function backtrace(moves::BandedArray{Int})
    taken_moves = DPMove[]
    i, j = moves.shape
    while i > 1 || j > 1
        m = moves[i, j]
        push!(taken_moves, DPMove(m))
        ii, jj = offsets[m]
        i += ii
        j += jj
    end
    return reverse(taken_moves)
end

function band_tolerance(Amoves::BandedArray{Int})
    nrows, ncols = size(Amoves)
    dist = nrows
    i, j = nrows, ncols
    while i > 1 || j > 1
        start, stop = row_range(Amoves, j)
        if start > 1
            dist = min(dist, abs(i - start))
        end
        if stop < nrows
            dist = min(dist, abs(i - stop))
        end
        m = Amoves[i, j]
        ii, jj = offsets[m]
        i += ii
        j += jj
    end
    start, stop = row_range(Amoves, j)
    if start > 1
        dist = min(dist, abs(i - start))
    end
    if stop < nrows
        dist = min(dist, abs(i - stop))
    end
    return dist
end

function moves_to_alignment_strings(moves::Vector{DPMove},
                                    t::AbstractString, s::AbstractString)
    aligned_t = Char[]
    aligned_s = Char[]
    i, j = (0, 0)
    for move in moves
        (ii, jj) = offsets[Int(move)]
        i -= ii
        j -= jj
        if move == dp_match
            push!(aligned_t, t[j])
            push!(aligned_s, s[i])
        elseif move == dp_ins
            push!(aligned_t, '-')
            push!(aligned_s, s[i])
        elseif move == dp_del
            push!(aligned_t, t[j])
            push!(aligned_s, '-')
        elseif move == dp_codon_ins
            append!(aligned_t, ['-', '-', '-'])
            append!(aligned_s, [s[i], s[i-1], s[i-2]])
        elseif move == dp_codon_del
            append!(aligned_t, [t[j], t[j-1], t[j-2]])
        end
    end
    return string(aligned_t...), string(aligned_s...)
end


function moves_to_indices(moves::Vector{DPMove},
                          t::AbstractString, s::AbstractString)
    tlen = length(t)
    slen = length(s)
    result = zeros(Int, tlen + 1)
    i, j = (1, 1)
    last_j = 0
    for move in moves
        if j > last_j
            result[j] = i
            last_j = j
        end
        (ii, jj) = offsets[Int(move)]
        i -= ii
        j -= jj
    end
    if j > last_j
        result[j] = i
        last_j = j
    end
    return result
end

function align_moves(t::AbstractString, s::AbstractString,
                     log_p::Vector{Float64},
                     scores::Scores,
                     bandwidth::Int)
    A, Amoves = forward_moves(t, s, log_p, scores, bandwidth)
    return backtrace(Amoves)
end

function align(t::AbstractString, s::AbstractString,
               log_p::Vector{Float64},
               scores::Scores,
               bandwidth::Int)
    moves = align_moves(t, s, log_p, scores, bandwidth)
    return moves_to_alignment_strings(moves, t, s)
end

function best_codons(template::AbstractString,
                     sequences::Vector{ASCIIString},
                     log_ps::Vector{Vector{Float64}},
                     scores::Scores,
                     bandwidth::Int;
                     Amoves::Vector{BandedArray{Int}}=BandedArray{Int}[])
    if length(Amoves) == length(sequences)
        moves = [backtrace(Am) for Am in Amoves]
    else
        moves = [align_moves(template, s, p, scores, bandwidth)
                 for (s, p) in zip(sequences, log_ps)]
    end
    indices = [moves_to_indices(m, template, s)
               for (m, s) in zip(moves, sequences)]

    results = []
    for j in 1:(length(template) + 1)
        cands = [Dict{Char, Float64}('A' => 0.0,
                                     'C' => 0.0,
                                     'G' => 0.0,
                                     'T' => 0.0)
                 for i in 1:3]
        for (seq, log_p, index) in zip(sequences, log_ps, indices)
            i = index[j]
            seq_i = i - 1
            if seq_i < length(seq) - 2
                for ii in 1:3
                    cands[ii][seq[seq_i + ii]] += log_p[seq_i + ii]
                end
            end
        end
        if maximum([minimum(values(c)) for c in cands]) == 0.0
            break
        end
        codon = [base_consensus(c) for c in cands]
        push!(results, CodonInsertion(j - 1, tuple(codon...)))
    end
    return results
end

function base_consensus(d::Dict{Char, Float64})
    return minimum([(v, k) for (k, v) in d])[2]
end


function has_single_indels(template::AbstractString,
                           reference::AbstractString,
                           ref_log_p::Vector{Float64},
                           ref_scores::Scores,
                           bandwidth::Int)
    has_right_length = length(template) % codon_length == 0
    moves = align_moves(template, reference, ref_log_p,
                        ref_scores, bandwidth)
    result = dp_ins in moves || dp_del in moves
    if !result && !has_right_length
        error("template length is not a multiple of three")
    end
    return result
end


function initial_state(template, seqs, lps, scores,
                       bandwidth)
    As = [forward(template, s, p, scores, bandwidth)
          for (s, p) in zip(seqs, lps)]
    Bs = [backward(template, s, p, scores, bandwidth)
          for (s, p) in zip(seqs, lps)]
    score = sum([A[end, end] for A in As])

    A_t = BandedArray(Float64, (1, 1), 1)
    B_t = BandedArray(Float64, (1, 1), 1)

    Amoves = BandedArray{Int}[]

    return State(score, template, A_t, B_t, As, Amoves,
                 Bs, initial_stage, false, bandwidth)
end


function recompute!(state::State, seqs::Vector{ASCIIString},
                    lps::Vector{Vector{Float64}},
                    scores::Scores,
                    reference::AbstractString,
                    ref_log_p::Vector{Float64},
                    ref_scores::Scores,
                    bandwidth_delta::Int,
                    recompute_As::Bool, recompute_Bs::Bool,
                    verbose::Int)
    tolerance = 0
    if recompute_As
        while tolerance < codon_length
            state.As = BandedArray{Float64}[]
            state.Amoves = BandedArray{Int}[]
            for (s, p) in zip(seqs, lps)
                As, Amoves = forward_moves(state.template, s, p, scores, state.bandwidth)
                push!(state.As, As)
                push!(state.Amoves, Amoves)
            end
            tolerance = minimum([band_tolerance(Am) for Am in state.Amoves])
            if ((state.stage == frame_correction_stage ||
                 state.stage == scoring_stage) &&
                length(reference) > 0)
                state.A_t, Amoves_t = forward_moves(state.template, reference, ref_log_p,
                                                    ref_scores, state.bandwidth)
                tolerance = min(tolerance, band_tolerance(Amoves_t))
            end
            if tolerance < codon_length
                state.bandwidth += bandwidth_delta
                if verbose > 1
                    println("NEW BANDWIDTH: $(state.bandwidth)")
                end
            end
        end
    end
    if recompute_Bs
        state.Bs = [backward(state.template, s, p, scores, state.bandwidth)
                    for (s, p) in zip(seqs, lps)]
        if ((state.stage == frame_correction_stage ||
             state.stage == scoring_stage) &&
            length(reference) > 0)
            state.B_t = backward(state.template, reference, ref_log_p,
                                 ref_scores, state.bandwidth)
        end
    end
    state.score = sum([A[end, end] for A in state.As])
    if ((state.stage == frame_correction_stage ||
         state.stage == scoring_stage) &&
        length(reference) > 0)
        state.score += state.A_t[end, end]
    end
end


"""convert per-proposal log differences to a per-base error rate"""
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
                        scores::Scores,
                        reference::AbstractString,
                        ref_log_p::Vector{Float64},
                        ref_scores::Scores)
    # `position_scores[i]` gives the following log probabilities
    # for `template[i]`: [A, C, G, T, -]
    position_scores = zeros(length(state.template), 5) + state.score
    # `insertion_scores[i]` gives the following log probabilities for an
    # insertion before `template[i]` of [A, C, G, T]
    insertion_scores = zeros(length(state.template) + 1, 4)

    # TODO: should we modify penalties before using reference?
    # - do not penalize mismatches
    # - use max indel penalty

    maxlen = maximum([length(s) for s in sequences])
    nrows = max(maxlen, length(reference)) + 1
    newcols = zeros(Float64, (nrows, 6))

    use_ref = (length(reference) > 0)
    for m in candstask(scoring_stage, state.template,
                       sequences, log_ps, state.Amoves, scores, state.bandwidth)
        score = score_proposal(m, state,
                               sequences, log_ps, scores,
                               use_ref,
                               reference, ref_log_p, ref_scores,
                               newcols)
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
                 scores::Scores;
                 reference::AbstractString="",
                 ref_log_p::Float64=0.0,
                 ref_scores::Scores=default_ref_scores,
                 ref_indel_penalty::Float64=-3.0,
                 min_ref_indel_score::Float64=-15.0,
                 bandwidth::Int=10, bandwidth_delta::Int=5,
                 min_dist::Int=15,
                 batch::Int=10, batch_threshold::Float64=0.05,
                 max_iters::Int=100, verbose::Int=0)
    if (scores.mismatch >= 0.0 || scores.mismatch == -Inf ||
        scores.insertion >= 0.0 || scores.insertion == -Inf ||
        scores.deletion >= 0.0 || scores.deletion == -Inf)
        error("scores must be between -Inf and 0.0")
    end
    if scores.codon_insertion > -Inf || scores.codon_deletion > -Inf
        error("error model cannot allow codon indels")
    end

    if bandwidth < 0
        error("bandwidth cannot be negative: $bandwidth")
    end

    if length(reference) > 0
        if ref_indel_penalty >= 0.0
            error("ref_indel_penalty must be < 0.0")
        end
        if min_ref_indel_score >= 0.0
            error("min_ref_indel_score must be < 0.0")
        end
        if (ref_log_p == -Inf || ref_log_p >= 0.0)
            error("ref_log_p=$ref_log_p but should be less than 0.0")
        end
        if (ref_scores.mismatch >= 0.0 ||
            ref_scores.insertion >= 0.0 ||
            ref_scores.deletion >= 0.0 ||
            ref_scores.codon_insertion >= 0.0 ||
            ref_scores.codon_deletion >= 0.0)
            error("ref scores cannot be >= 0")
        end
        if (ref_scores.mismatch == -Inf ||
            ref_scores.insertion == -Inf ||
            ref_scores.deletion == -Inf ||
            ref_scores.codon_insertion == -Inf ||
            ref_scores.codon_deletion == -Inf)
            error("ref scores cannot be -Inf")
        end
        if (ref_scores.insertion < min_ref_indel_score ||
            ref_scores.deletion < min_ref_indel_score)
            error("ref indel scores are less than specified minimum")
        end

    end

    ref_log_p_vec = fill(ref_log_p, length(reference))

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
    state = initial_state(template, seqs, lps, scores, bandwidth)
    empty_ref = length(reference) == 0

    if verbose > 1
        println(STDERR, "initial score: $(state.score)")
    end
    n_proposals = Vector{Int}[]
    consensus_lengths = Int[length(template)]
    consensus_stages = ["", "", ""]

    stage_iterations = zeros(Int, Int(typemax(Stage)))
    for i in 1:max_iters
        stage_iterations[Int(state.stage)] += 1
        old_template = state.template
        old_score = state.score
        if verbose > 1
            println(STDERR, "iteration $i : $(state.stage)")
        end

        penalties_increased = false

        candidates = getcands(state, seqs, lps, scores,
                              reference, ref_log_p_vec, ref_scores)

        recompute_As = true
        if length(candidates) == 0
            if verbose > 1
                println(STDERR, "  no candidates found")
            end
            push!(n_proposals, zeros(Int, Int(typemax(DPMove))))
            consensus_stages[Int(state.stage)] = state.template
            if state.stage == initial_stage
                if empty_ref
                    state.converged = true
                    break
                end
                state.stage = frame_correction_stage
            elseif state.stage == frame_correction_stage
                if !has_single_indels(state.template, reference, ref_log_p_vec, ref_scores, bandwidth)
                    consensus_ref = state.template
                    state.stage = refinement_stage
                elseif (ref_scores.insertion == min_ref_indel_score ||
                        ref_scores.deletion == min_ref_indel_score)
                    if verbose > 1
                        println(STDERR, "  alignment had single indels but indel scores already minimized.")
                    end
                    state.stage = refinement_stage
                else
                    penalties_increased = true
                    # TODO: this is not probabilistically correct
                    ref_scores = Scores(ref_scores.mismatch,
                                        max(ref_scores.insertion + ref_indel_penalty,
                                            min_ref_indel_score),
                                        max(ref_scores.deletion + ref_indel_penalty,
                                            min_ref_indel_score),
                                        ref_scores.codon_insertion,
                                        ref_scores.codon_deletion)
                    if verbose > 1
                        println(STDERR, "  alignment to reference had single indels. increasing penalty.")
                    end
                end
            elseif state.stage == refinement_stage
                state.converged = true
                break
            else
                error("unknown stage: $(state.stage)")
            end
        else
            if verbose > 1
                println(STDERR, "  found $(length(candidates)) candidate")
            end
            chosen_cands = choose_candidates(candidates, min_dist)
            if verbose > 1
                println(STDERR, "  filtered to $(length(chosen_cands)) candidates")
            end
            state.template = apply_proposals(old_template,
                                             Proposal[c.proposal
                                                      for c in chosen_cands])
            recompute!(state, seqs, lps, scores,
                       reference, ref_log_p_vec, ref_scores,
                       bandwidth_delta, true, false, verbose)
            # detect if a single proposal is better
            # note: this may not always be correct, because score_proposal() is not exact
            if length(chosen_cands) > 1 &&
                (state.score < chosen_cands[1].score ||
                 isapprox(state.score, chosen_cands[1].score))
                if verbose > 1
                    println(STDERR, "  rejecting multiple candidates in favor of best")
                end
                chosen_cands = CandProposal[chosen_cands[1]]
                state.template = apply_proposals(old_template,
                                                 Proposal[c.proposal
                                                          for c in chosen_cands])
            else
                # no need to recompute unless batch changes
                recompute_As = false
            end
            proposal_counts = [length(filter(c -> (typeof(c.proposal) == t),
                                             chosen_cands))
                               for t in [Substitution, Insertion, Deletion, CodonInsertion, CodonDeletion]]
            push!(n_proposals, proposal_counts)
        end
        push!(consensus_lengths, length(state.template))
        if batch < length(sequences)
            indices = rand(1:length(sequences), batch)
            seqs = sequences[indices]
            lps = log_ps[indices]
            recompute_As = true
        end
        recompute!(state, seqs, lps, scores,
                   reference, ref_log_p_vec, ref_scores,
                   bandwidth_delta, recompute_As, true, verbose)
        if verbose > 1
            println(STDERR, "  score: $(state.score) ($(state.stage))")
        end
        if (state.score < old_score &&
            stage_iterations[Int(state.stage)] > 0 &&
            !penalties_increased &&
            batch == length(sequences))
            if verbose > 1
                # caused by bandwidth being too small.
                # this should not usually happen, because bandwidth
                # should be increased before this.
                println(STDERR, "  WARNING: not using batches, but score decreased.")
            end
        end
        if ((state.score - old_score) / old_score > batch_threshold &&
            !penalties_increased &&
            batch < length(sequences) &&
            stage_iterations[Int(state.stage)] > 0)
            batch = min(batch + base_batch, length(sequences))
            indices = rand(1:length(sequences), batch)
            seqs = sequences[indices]
            lps = log_ps[indices]
            recompute!(state, seqs, lps, scores,
                       reference, ref_log_p_vec, ref_scores,
                       bandwidth_delta, true, true, verbose)
            if verbose > 1
                println(STDERR, "  increased batch size to $batch. new score: $(state.score)")
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
                "ref_scores" => ref_scores,
                "consensus_stages" => consensus_stages,
                "n_proposals" => transpose(hcat(n_proposals...)),
                "consensus_lengths" => consensus_lengths,
                "bandwidth" => state.bandwidth,
                )

    # FIXME: recomputing for all sequences is costly, but using batch is less accurate
    recompute!(state, seqs, lps, scores,
               reference, ref_log_p_vec, ref_scores,
               bandwidth_delta, true, true, verbose)
    base_probs, ins_probs = estimate_probs(state, seqs, lps, scores,
                                           reference, ref_log_p_vec, ref_scores)
    return state.template, base_probs, ins_probs, info
end

"""
Alternate quiver2() using BioJulia types.

"""
function quiver2(template::DNASequence,
                 sequences::Vector{DNASequence},
                 phreds::Vector{Vector{Int8}},
                 scores::Scores;
                 reference::DNASequence=DNASequence(""),
                 kwargs...)
    new_reference = convert(ASCIIString, reference)
    new_template = convert(ASCIIString, template)
    new_sequences = ASCIIString[convert(ASCIIString, s) for s in sequences]
    (result, base_probs,
     insertion_probs, info) = quiver2(new_template, new_sequences, phreds, scores;
                                      reference=new_reference,
                                      kwargs...)
    info["consensus_stages"] = DNASequence[DNASequence(s) for s in info["consensus_stages"]]
    return (DNASequence(result), base_probs, insertion_probs, info)
end

end
