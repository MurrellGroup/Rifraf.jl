module Model

using Bio.Seq

using Levenshtein

using Quiver2.BandedArrays
using Quiver2.Proposals
using Quiver2.Util

import Base.length
import Base.reverse

export quiver2, ErrorModel, Scores, normalize, estimate_point_probs, estimate_indel_probs

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

"""Convenience type for passing around a sequence and its per-base
error probabilities.

"""
type PString
    seq::String
    error_log_p::Vector{Float64}
    match_log_p::Vector{Float64}
    bandwidth::Int

    function PString(seq::String, error_log_p::Vector{Float64}, bandwidth::Int)
        if bandwidth < 1
            error("bandwidth must be positive")
        end

        if length(seq) != length(error_log_p)
            error("length mismatch")
        end
        if length(seq) == 0
            return new(seq, Float64[], Float64[])
        end
        if minimum(error_log_p) == -Inf
            error("a log error probability is negative infinity")
        end
        if maximum(error_log_p) > 0.0
            bad_value = maximum(error_log_p)
            error("a log error probability is > 0: $bad_value")
        end
        match_log_p = log10(1.0 - exp10(error_log_p))
        return new(seq, error_log_p, match_log_p, bandwidth)
    end
end

function PString(seq::String, phreds::Vector{Int8}, bandwidth::Int)
    error_log_p = phred_to_log_p(phreds)
    return PString(seq, error_log_p, bandwidth)
end

function length(s::PString)
    return length(s.seq)
end

function reverse(s::PString)
    return PString(reverse(s.seq), reverse(s.error_log_p), s.bandwidth)
end


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

function Scores(errors::ErrorModel;
                mismatch::Float64=0.0,
                insertion::Float64=0.0,
                deletion::Float64=0.0)
    args = Float64[errors.mismatch,
                   errors.insertion,
                   errors.deletion,
                   errors.codon_insertion,
                   errors.codon_deletion]
    m, i, d, ci, cd = log10(Util.normalize(args))
    return Scores(m + mismatch,
                  i + insertion,
                  d + deletion,
                  ci + 3 * insertion,
                  cd + 3 * deletion)
end


immutable EstErrorProbs
    sub::Array{Float64, 2}
    del::Array{Float64, 1}
    ins::Array{Float64, 2}
end


# default to invalid scores, to force user to set them
const default_ref_scores = Scores(0.0, 0.0, 0.0, 0.0, 0.0)

# just to avoid magical constants in code
const codon_length = 3

type State
    score::Float64
    consensus::String
    A_t::BandedArray{Float64}
    B_t::BandedArray{Float64}
    As::Vector{BandedArray{Float64}}
    Amoves::Vector{BandedArray{Int}}
    Bs::Vector{BandedArray{Float64}}
    stage::Stage
    converged::Bool
end

@enum DPMove dp_none=0 dp_match=1 dp_ins=2 dp_del=3 dp_codon_ins=4 dp_codon_del=5


# all offsets are relative to the consensus.
# so an insertion is a base NOT in the consensus.
const offsets = ([1, 1],  # sub
                 [1, 0],  # insertion
                 [0, 1],  # deletion
                 [3, 0],  # codon insertion
                 [0, 3])  # codon deletion

const bases = "ACGT"
const bases_set = Set{Char}(bases)
const baseints = Dict('A' => 1,
                      'C' => 2,
                      'G' => 3,
                      'T' => 4,
                      )

const empty_array = Array(Float64, (0, 0))

function move_scores(t_base::Char,
                     s_base::Char,
                     seq_i::Int,
                     error_log_p::Vector{Float64},
                     match_log_p::Vector{Float64},
                     scores::Scores;
                     match_mult::Float64=0.0)
    cur_i = max(seq_i, 1)
    cur_log_p = error_log_p[cur_i]
    next_log_p = error_log_p[min(seq_i + 1, length(error_log_p))]

    match_score = 0.0
    if s_base == t_base
        match_score = match_log_p[cur_i]
    else
        # match_mult makes mismatches slightly better, proportional
        # to that base's qv score.
        # this is to break the symmetry of insertions and mismatches, and
        # encourage insertions of low-quality bases.
        # TODO: think of a less-hacky solution
        match_score = cur_log_p + scores.mismatch
        if match_mult > 0.0
            match_score -= cur_log_p * match_mult
        end
    end
    ins_score = cur_log_p + scores.insertion
    del_score = max(cur_log_p, next_log_p) + scores.deletion
    return match_score, ins_score, del_score
end

function codon_move_scores(seq_i::Int,
                           error_log_p::Vector{Float64},
                           scores::Scores)
    # we're moving INTO seq_i. so need previous three
    start = max(1, seq_i-2)
    stop = min(seq_i, length(error_log_p))
    max_p = start <= stop ? maximum(error_log_p[start:stop]) : -Inf
    codon_ins_score =  max_p + scores.codon_insertion
    cur_log_p = error_log_p[max(seq_i, 1)]
    next_log_p = error_log_p[min(seq_i + 1, length(error_log_p))]
    codon_del_score = max(cur_log_p, next_log_p) + scores.codon_deletion
    return codon_ins_score, codon_del_score
end

function update_helper(newcols::Array{Float64, 2},
                       A::BandedArray{Float64},
                       i::Int, j::Int, acol::Int,
                       move::DPMove, move_score::Float64,
                       final_score::Float64, final_move::DPMove)
    offset = offsets[Int(move)]
    prev_i = i - offset[1]
    prev_j = j - offset[2]

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
                pseq::PString,
                scores::Scores;
                newcols::Array{Float64, 2}=empty_array,
                acol=-1, trim=false,
                skew_matches=false)
    result = (-Inf, dp_none)
    # TODO: this cannot make mismatches preferable to codon indels
    match_mult = skew_matches ? 0.1 : 0.0
    match_score, ins_score, del_score = move_scores(t_base, s_base, i-1,
                                                    pseq.error_log_p,
                                                    pseq.match_log_p,
                                                    scores;
                                                    match_mult=match_mult)
    # allow terminal insertions for free
    if trim && (j == 1)
        ins_score = 0.0
    end
    if trim && (j == size(A)[2])
        ins_score = 0.0
    end
    result = update_helper(newcols, A, i, j, acol, dp_match, match_score, result...)
    result = update_helper(newcols, A, i, j, acol, dp_ins, ins_score, result...)
    result = update_helper(newcols, A, i, j, acol, dp_del, del_score, result...)
    codon_ins_score, codon_del_score = codon_move_scores(i-1, pseq.error_log_p, scores)

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


"""Does backtracing to find best alignment."""
function forward_moves(t::String, s::PString,
                       scores::Scores;
                       trim::Bool=false,
                       skew_matches::Bool=false)
    # FIXME: code duplication with forward_codon(). This is done in a
    # seperate function to keep return type stable and avoid
    # allocating the `moves` array unnecessarily.
    result = BandedArray(Float64, (length(s) + 1, length(t) + 1), s.bandwidth)
    moves = BandedArray(Int, result.shape, s.bandwidth)
    moves[1, 1] = Int(dp_none)
    nrows, ncols = size(result)
    for j = 1:ncols
        start, stop = row_range(result, j)
        for i = start:stop
            if i == 1 && j == 1
                continue
            end
            sbase = i > 1 ? s.seq[i-1] : 'X'
            tbase = j > 1 ? t[j-1] : 'X'
            x = update(result, i, j, sbase, tbase,
                       s, scores; trim=trim,
                       skew_matches=skew_matches)
            result[i, j] = x[1]
            moves[i, j] = x[2]
        end
    end
    return result, moves
end


"""
F[i, j] is the log probability of aligning s[1:i-1] to t[1:j-1].

"""
function forward(t::String, s::PString,
                 scores::Scores)
    result = BandedArray(Float64, (length(s) + 1, length(t) + 1), s.bandwidth)
    nrows, ncols = size(result)
    for j = 1:ncols
        start, stop = row_range(result, j)
        for i = start:stop
            if i == 1 && j == 1
                continue
            end
            sbase = i > 1 ? s.seq[i-1] : 'X'
            tbase = j > 1 ? t[j-1] : 'X'
            x = update(result, i, j, sbase, tbase,
                       s, scores)
            result[i, j] = x[1]
        end
    end
    return result::BandedArray{Float64}
end


"""
B[i, j] is the log probability of aligning s[i:end] to t[j:end].

"""
function backward(t::String, s::PString,
                  scores::Scores)
    result = forward(reverse(t), reverse(s), scores)
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


function get_sub_consensus(proposal::Proposal, seq::String,
                           next_posn::Int, n_after::Int)
    # next valid position in sequence after this proposal
    t = typeof(proposal)
    pos = proposal.pos
    prefix = ""
    stop = min(next_posn + n_after - 1, length(seq))
    suffix = seq[next_posn:stop]
    if t in (Substitution, Insertion)
        prefix = string(proposal.base)
    end
    return string(prefix, suffix)
end

function seq_score_deletion(A::BandedArray{Float64}, B::BandedArray{Float64},
                            acol::Int, bcol::Int)
    Acol = sparsecol(A, acol)
    Bcol = sparsecol(B, bcol)
    (amin, amax), (bmin, bmax) = equal_ranges(row_range(A, acol),
                                              row_range(B, bcol))
    asub = view(Acol, amin:amax)
    bsub = view(Bcol, bmin:bmax)
    return summax(asub, bsub)
end


const boffsets = Dict(Substitution => 2,
                      Insertion => 1,
                      Deletion => 2)

function score_nocodon(proposal::Proposal,
                       A::BandedArray{Float64}, B::BandedArray{Float64},
                       pseq::PString,
                       scores::Scores,
                       newcols::Array{Float64, 2})
    t = typeof(proposal)
    if t == Deletion
        # nothing to recompute
        acol = proposal.pos
        bcol = proposal.pos + 1
        return seq_score_deletion(A, B, acol, bcol)
    end
    # need to compute new columns
    nrows, ncols = size(A)

    # last valid A column
    acol = proposal.pos + (t == Substitution ? 0 : 1)
    new_acol = acol + 1
    # new base
    amin, amax = row_range(A, min(new_acol, ncols))
    for i in amin:amax
        seq_base = i > 1 ? pseq.seq[i-1] : 'X'
        x = update(A, i, new_acol,
                   seq_base, proposal.base,
                   pseq, scores;
                   newcols=newcols, acol=acol)
        newcols[i, 1] = x[1]
    end

    # add up results
    imin, imax = row_range(A, min(new_acol, ncols))
    Acol = view(newcols, imin:imax, 1)

    bj = proposal.pos + 1
    Bcol = sparsecol(B, bj)
    (amin, amax), (bmin, bmax) = equal_ranges((imin, imax),
                                              row_range(B, bj))
    asub = view(Acol, amin:amax)
    bsub = view(Bcol, bmin:bmax)
    score = summax(asub, bsub)
    if score == -Inf
        error("failed to compute a valid score")
    end
    return score
end

function seq_score_proposal(proposal::Proposal,
                            A::BandedArray{Float64}, B::BandedArray{Float64},
                            consensus::String,
                            pseq::PString,
                            scores::Scores,
                            newcols::Array{Float64, 2})
    codon_moves = (scores.codon_insertion > -Inf ||
                   scores.codon_deletion > -Inf)
    if !codon_moves
        return score_nocodon(proposal, A, B,
                             pseq, scores, newcols)
    end
    t = typeof(proposal)
    # last valid column of A
    acol_offset = t == Insertion ? 0 : -1
    acol = proposal.pos + acol_offset + 1

    # first column of B to use
    first_bcol = acol + boffsets[t]
    # last column of B to use
    last_bcol = first_bcol + codon_length - 1

    if t == Deletion
        n_del = (t == Deletion ? 1 : codon_length)
        if acol == (size(A)[2] - n_del)
            # suffix deletions do not need recomputation
            return A[end, end - n_del]
        end
    end

    # number of bases changed/inserted
    n_bases = (t == Deletion ? 0 : 1)
    # number of columns after recomputed columns to also recompute.
    n_after = codon_length

    # if we'd go to or past the last column of B, just recompute the
    # rest of A
    nrows, ncols = size(A)
    just_a = last_bcol >= ncols
    next_posn = proposal.pos + 1
    if just_a
        # go to end of consensus
        n_after = length(consensus) - next_posn + 1
    end

    if n_bases == 0 && n_after == 0
        error("no new columns need to be recomputed.")
    end

    sub_consensus = get_sub_consensus(proposal, consensus,
                                      next_posn, n_after)
    n_new = n_bases + n_after
    # compute new columns
    for j in 1:n_new
        range_col = min(acol + j, ncols)
        amin, amax = row_range(A, range_col)
        for i in amin:amax
            seq_base = i > 1 ? pseq.seq[i-1] : 'X'
            x = update(A, i, acol + j,
                       seq_base, sub_consensus[j],
                       pseq, scores;
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
            asub = view(Acol, amin:amax)
            bsub = view(Bcol, bmin:bmax)
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
    posns = Set{Int}()
    for c in sort(candidates, by=(c) -> c.score, rev=true)
        if any(Bool[(abs(c.proposal.pos - p) < min_dist) for p in posns])
            continue
        end
        push!(posns, c.proposal.pos)
        push!(final_cands, c)
    end
    return final_cands
end


function score_proposal(m::Proposal,
                        state::State,
                        sequences::Vector{PString},
                        scores::Scores,
                        use_ref::Bool,
                        reference::PString,
                        ref_scores::Scores,
                        newcols::Array{Float64, 2})
    score = 0.0
    for si in 1:length(sequences)
        score += seq_score_proposal(m, state.As[si], state.Bs[si], state.consensus,
                                    sequences[si], scores, newcols)
    end
    if use_ref
        score += seq_score_proposal(m, state.A_t, state.B_t, state.consensus,
                                    reference, ref_scores, newcols)
    end
    return score
end


function all_proposals(stage::Stage,
                       consensus::String,
                       sequences::Vector{PString},
                       Amoves::Vector{BandedArray{Int}},
                       scores::Scores,
                       indel_correction_only::Bool)
    len = length(consensus)
    function _it()
        if stage != frame_correction_stage || !indel_correction_only
            # substitutions
            for j in 1:len
                for base in "ACGT"
                    if consensus[j] != base
                        produce(Substitution(j, base))
                    end
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
    end
    Task(_it)
end


function moves_to_proposals(moves::Vector{DPMove},
                            consensus::String, seq::PString)
    proposals = Proposal[]
    i, j = (0, 0)
    for move in moves
        (ii, jj) = offsets[Int(move)]
        i += ii
        j += jj

        score = seq.match_log_p[max(i, 1)]
        next_score = seq.match_log_p[min(i + 1, length(seq))]
        del_score = min(score, next_score)

        if move == dp_match
            if seq.seq[i] != consensus[j]
                push!(proposals, Substitution(j, seq.seq[i]))
            end
        elseif move == dp_ins
            push!(proposals, Insertion(j, seq.seq[i]))
        elseif move == dp_del
            push!(proposals, Deletion(j))
        end
    end
    return proposals
end

"""Only get proposals that appear in at least one alignment"""
function alignment_proposals(state::State,
                             sequences::Vector{PString},
                             do_subs::Bool,
                             do_indels::Bool)
    result = Set{Proposal}()
    for (Amoves, seq) in zip(state.Amoves, sequences)
        moves = backtrace(Amoves)
        for proposal in moves_to_proposals(moves, state.consensus, seq)
            if (typeof(proposal) == Substitution && do_subs) || do_indels
                push!(result, proposal)
            end
        end
    end
    return collect(result)
end


"""Use model surgery heuristic to get proposals with positive score deltas."""
function surgery_proposals(state::State,
                           sequences::Vector{PString},
                           scores::Scores,
                           do_subs::Bool,
                           do_indels::Bool)
    # FIXME: modularize this function and test each part
    sub_deltas = zeros(Float64, (length(state.consensus), 4))
    del_deltas = zeros(Float64, (length(state.consensus)))
    ins_deltas = zeros(Float64, (length(state.consensus) + 1, 4))
    for (Amoves, seq) in zip(state.Amoves, sequences)
        moves = backtrace(Amoves)
        seq_idx = 0
        cons_idx = 0

        # map each base to its maximum delta
        # start with deletion scores for each
        del_score = seq.error_log_p[1] + scores.deletion

        # push insertions and deletions to the end, on the fly
        # insertions that match consensus should *keep* getting pushed
        # deletions in a poly-base run should also get pushed
        insertion_bases = Dict{Char, Float64}(
                                              'A' => del_score,
                                              'C' => del_score,
                                              'G' => del_score,
                                              'T' => del_score,
                                              )

        del_base = '-'
        del_idx = 0
        pushed_del_score = -Inf

        for (aln_idx, move) in enumerate(moves)
            cbase = '-'
            sbase = '-'
            if move in (dp_match, dp_ins)
                seq_idx += 1
                sbase = seq.seq[seq_idx]
            end
            if move in (dp_match, dp_del)
                cons_idx += 1
                cbase = state.consensus[cons_idx]
            end

            match_score = -Inf
            error_score = Inf
            if seq_idx > 0
                match_score = seq.match_log_p[seq_idx]
                error_score = seq.error_log_p[seq_idx]
            end
            mismatch_score = error_score + scores.mismatch
            insertion_score = error_score + scores.insertion
            deletion_score = max(seq.error_log_p[max(seq_idx, 1)],
                                 seq.error_log_p[min(seq_idx + 1, length(seq.error_log_p))])

            # handle pushed deletions
            if del_base != '-' && del_base != cbase
                del_deltas[del_idx] += pushed_del_score
                del_base = '-'
                del_idx = 0
                pushed_del_score = -Inf
            end
            if move != dp_ins
                del_score = max(seq.error_log_p[max(seq_idx, 1)],
                                seq.error_log_p[min(seq_idx + 1, length(seq.error_log_p))])

                # handle insertions before this position.
                for (base, delta) in insertion_bases
                    # if current base equals attempted insertion,
                    # keep pushing it
                    if cbase != base
                        ins_deltas[cons_idx, baseints[base]] += delta
                        insertion_bases[base] = del_score
                    end
                end
            end

            if move == dp_ins
                # consider insertion proposals.
                # push all insertion proposals to the end
                for (baseint, new_base) in enumerate(bases)
                    # try all insertions, converting dp_ins to (mis)match
                    sbase = seq.seq[seq_idx]
                    new_score = (sbase == new_base) ? match_score : mismatch_score
                    insertion_bases[new_base] = max(insertion_bases[new_base],
                                                    new_score - insertion_score)
                end
            elseif move == dp_match
                # consider all substitution proposals
                for (baseint, new_base) in enumerate(bases)
                    if new_base == cbase
                        continue
                    end
                    if sbase == new_base
                        # proposal would help
                        sub_deltas[cons_idx, baseint] += (match_score - mismatch_score)
                    elseif sbase == cbase
                        # proposal would hurt
                        sub_deltas[cons_idx, baseint] += (mismatch_score - match_score)
                    else
                        # proposal would still be a mismatch. do nothing.
                    end
                end
                # consider deleting the consensus base.
                # This would convert (mis)match to insertion
                del_base = cbase
                del_idx = cons_idx
                pushed_del_score = max(pushed_del_score,
                                       insertion_score - (sbase == cbase ? match_score : mismatch_score))
            elseif move == dp_del
                # no reason to consider substitutions. delta = 0.
                # consider deletion proposal. would convert deletion
                # to no-op.
                del_base = cbase
                del_idx = cons_idx
                pushed_del_score = max(pushed_del_score,
                                       0 - deletion_score)
            end
        end
        # handle deletion at end
        if del_base != '-'
            del_deltas[end] += pushed_del_score
        end

        # handle insertions at the end
        for (base, delta) in insertion_bases
            ins_deltas[end, baseints[base]] += delta
        end
    end

    # only return proposals with positive deltas
    result = Proposal[]
    deltas = Float64[]
    nrows, ncols = size(sub_deltas)
    if do_subs
        for i in 1:nrows, j in 1:ncols
            if sub_deltas[i, j] > 0
                push!(result, Substitution(i, bases[j]))
                push!(deltas, sub_deltas[i, j])
            end
        end
    end
    if do_indels
        for i in 1:nrows
            if del_deltas[i] > 0
                push!(result, Deletion(i))
                push!(deltas, del_deltas[i])
            end
        end
        for i in 1:(nrows+1), j in 1:ncols
            if ins_deltas[i, j] > 0
                push!(result, Insertion(i - 1, bases[j]))
                push!(deltas, ins_deltas[i, j])
            end
        end
    end
    return result, deltas
end

function get_candidate_proposals(state::State,
                                 sequences::Vector{PString},
                                 scores::Scores,
                                 reference::PString,
                                 ref_scores::Scores,
                                 do_alignment_proposals::Bool,
                                 do_surgery_proposals::Bool,
                                 trust_proposals::Bool,
                                 indel_correction_only::Bool)
    candidates = CandProposal[]
    use_ref = (state.stage == frame_correction_stage)

    maxlen = maximum(length(s) for s in sequences)
    nrows = max(maxlen, length(reference)) + 1
    newcols = zeros(Float64, (nrows, codon_length + 1))

    proposals = all_proposals(state.stage, state.consensus, sequences,
                              state.Amoves, scores,
                              indel_correction_only)
    if (state.stage == initial_stage ||
        state.stage == refinement_stage) && do_surgery_proposals
        do_indels = state.stage in (initial_stage, frame_correction_stage)
        do_subs = true
        if state.stage == frame_correction_stage && indel_correction_only
            # FIXME: this never happens
            do_subs = false
        end
        proposals, deltas = surgery_proposals(state, sequences, scores, do_subs, do_indels)
        if trust_proposals
            for (p, d) in zip(proposals, deltas)
                push!(candidates, CandProposal(p, state.score + d))
            end
            return candidates
        end
    elseif (state.stage == initial_stage ||
        state.stage == refinement_stage) && do_alignment_proposals
        do_indels = state.stage in (initial_stage, frame_correction_stage)
        do_subs = true
        if state.stage == frame_correction_stage && indel_correction_only
            do_subs = false
        end
        proposals = alignment_proposals(state, sequences, do_subs, do_indels)
    end

    for p in proposals
        score = score_proposal(p, state, sequences, scores, use_ref,
                               reference, ref_scores, newcols)
        if score > state.score
            push!(candidates, CandProposal(p, score))
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
        i -= ii
        j -= jj
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
        i -= ii
        j -= jj
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
                                    t::String, s::String)
    aligned_t = Char[]
    aligned_s = Char[]
    i, j = (0, 0)
    for move in moves
        (ii, jj) = offsets[Int(move)]
        i += ii
        j += jj
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
            append!(aligned_s, ['-', '-', '-'])
        end
    end
    return string(aligned_t...), string(aligned_s...)
end


""" Compute index vector mapping from position in `t` to position in
`s`.

"""
function moves_to_indices(moves::Vector{DPMove},
                          tlen::Int, slen::Int)
    result = zeros(Int, tlen + 1)
    i, j = (1, 1)
    last_j = 0
    for move in moves
        if j > last_j
            result[j] = i
            last_j = j
        end
        (ii, jj) = offsets[Int(move)]
        i += ii
        j += jj
    end
    if j > last_j
        result[j] = i
        last_j = j
    end
    return result
end

function align_moves(t::String, s::PString,
                     scores::Scores;
                     trim::Bool=false,
                     skew_matches::Bool=false)
    A, Amoves = forward_moves(t, s, scores;
                              trim=trim, skew_matches=skew_matches)
    return backtrace(Amoves)
end

function align(t::String, s::PString,
               scores::Scores;
               trim::Bool=false,
               skew_matches::Bool=false)
    moves = align_moves(t, s, scores, trim=trim,
                        skew_matches=skew_matches)
    return moves_to_alignment_strings(moves, t, s.seq)
end

function align(t::String, s::String, phreds::Vector{Int8},
               scores::Scores,
               bandwidth::Int;
               trim::Bool=false,
               skew_matches::Bool=false)
    moves = align_moves(t, PString(s, phreds, bandwidth), scores,
                        trim=trim, skew_matches=skew_matches)
    return moves_to_alignment_strings(moves, t, s)
end


function base_consensus(d::Dict{Char, Float64})
    return minimum((v, k) for (k, v) in d)[2]
end


function has_single_indels(consensus::String,
                           reference::PString,
                           ref_scores::Scores)
    has_right_length = length(consensus) % codon_length == 0
    moves = align_moves(consensus, reference, ref_scores)
    result = dp_ins in moves || dp_del in moves
    if !result && !has_right_length
        error("consensus length is not a multiple of three")
    end
    return result
end


function single_indel_proposals(reference::String,
                                consensus::PString,
                                ref_scores::Scores)
    moves = align_moves(reference, consensus, ref_scores;
                        skew_matches=true)
    results = Proposal[]
    cons_idx = 0
    ref_idx = 0
    for move in moves
        if move == dp_match
            cons_idx += 1
            ref_idx += 1
        elseif move == dp_ins
            cons_idx += 1
            push!(results, Deletion(cons_idx))
        elseif move == dp_del
            ref_idx += 1
            push!(results, Insertion(cons_idx, reference[ref_idx]))
        elseif move == dp_codon_ins
            cons_idx += 3
        elseif move == dp_codon_del
            ref_idx += 3
        end
    end
    return results
end


function initial_state(consensus::String, seqs::Vector{PString})
    if length(consensus) == 0
        # choose highest-quality sequence
        idx = indmin([Util.logsumexp10(s.error_log_p)
                      for s in seqs])
        consensus = seqs[idx].seq
    end

    A_t = BandedArray(Float64, (1, 1), 1)
    B_t = BandedArray(Float64, (1, 1), 1)

    As = BandedArray{Float64}[]
    Bs = BandedArray{Float64}[]
    Amoves = BandedArray{Int}[]

    return State(0.0, consensus, A_t, B_t, As, Amoves,
                 Bs, initial_stage, false)
end


function recompute!(state::State, seqs::Vector{PString},
                    scores::Scores, reference::PString,
                    ref_scores::Scores, bandwidth_mult::Int,
                    recompute_As::Bool, recompute_Bs::Bool,
                    verbose::Int, use_ref_for_qvs::Bool)
    if recompute_As
        state.As = BandedArray{Float64}[]
        state.Amoves = BandedArray{Int}[]
        for s in seqs
            As, Amoves = forward_moves(state.consensus, s, scores)
            while band_tolerance(Amoves) < codon_length
                s.bandwidth *= bandwidth_mult
                As, Amoves = forward_moves(state.consensus, s, scores)
            end
            push!(state.As, As)
            push!(state.Amoves, Amoves)
        end
        if (((state.stage == frame_correction_stage) ||
             ((state.stage == scoring_stage) &&
              (use_ref_for_qvs))) &&
              (length(reference) > 0))
            state.A_t, Amoves_t = forward_moves(state.consensus, reference, ref_scores)
            while band_tolerance(Amoves_t) < codon_length
                reference.bandwidth *= bandwidth_mult
                state.A_t, Amoves_t = forward_moves(state.consensus, reference, ref_scores)
            end
        end
    end
    if recompute_Bs
        state.Bs = [backward(state.consensus, s, scores) for s in seqs]
        if (((state.stage == frame_correction_stage) ||
             ((state.stage == scoring_stage) &&
              (use_ref_for_qvs))) &&
            (length(reference) > 0))
            state.B_t = backward(state.consensus, reference, ref_scores)
        end
    end
    state.score = sum(A[end, end] for A in state.As)
    if (((state.stage == frame_correction_stage) ||
         ((state.stage == scoring_stage) &&
          (use_ref_for_qvs))) &&
        (length(reference) > 0))
        state.score += state.A_t[end, end]
    end
end


"""convert per-proposal log differences to a per-base error rate"""
function normalize_log_differences(sub_scores,
                                   del_scores,
                                   ins_scores,
                                   state_score)
    # per-base insertion score is mean of neighboring insertions
    pos_scores = hcat(sub_scores, del_scores)
    pos_exp = exp10(pos_scores)
    pos_probs = broadcast(/, pos_exp, sum(pos_exp, 2))
    ins_exp = exp10(ins_scores)
    ins_probs = broadcast(/, ins_exp, exp10(state_score) + sum(ins_exp, 2))
    sub_probs = pos_probs[:, 1:4]
    del_probs = pos_probs[:, 5]
    return EstErrorProbs(sub_probs, del_probs, ins_probs)
end


function estimate_probs(state::State,
                        sequences::Vector{PString},
                        scores::Scores,
                        reference::PString,
                        ref_scores::Scores,
                        use_ref_for_qvs::Bool)
    # `sub_scores[i]` gives the following log probabilities
    # for `consensus[i]`: [A, C, G, T, -]
    sub_scores = zeros(length(state.consensus), 4) + state.score
    del_scores = zeros(length(state.consensus)) + state.score
    # `ins_scores[i]` gives the following log probabilities for an
    # insertion before `consensus[i]` of [A, C, G, T]
    ins_scores = zeros(length(state.consensus) + 1, 4)

    # TODO: should we modify penalties before using reference?
    # - do not penalize mismatches
    # - use max indel penalty

    maxlen = maximum(length(s) for s in sequences)
    nrows = max(maxlen, length(reference)) + 1
    newcols = zeros(Float64, (nrows, codon_length + 1))

    use_ref = (length(reference) > 0) && use_ref_for_qvs
    # do not want to penalize indels too harshly, otherwise consensus
    # appears too likely.
    # TODO: how to set scores correctly for estimating qvs?
    qv_ref_scores = Scores(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    for m in all_proposals(scoring_stage, state.consensus, sequences,
                           state.Amoves, scores, false)
        score = score_proposal(m, state,
                               sequences, scores,
                               use_ref,
                               reference, qv_ref_scores,
                               newcols)
        if typeof(m) == Substitution
            sub_scores[m.pos, baseints[m.base]] = score
        elseif typeof(m) == Deletion
            del_scores[m.pos] = score
        elseif typeof(m) == Insertion
            ins_scores[m.pos + 1, baseints[m.base]] = score
        end
    end
    max_score = maximum([maximum(sub_scores),
                         maximum(del_scores),
                         maximum(ins_scores)])
    sub_scores -= max_score
    del_scores -= max_score
    ins_scores -= max_score
    if maximum(sub_scores) > 0.0
        error("sub scores cannot be positive")
    end
    if maximum(del_scores) > 0.0
        error("deletion scores cannot be positive")
    end
    if maximum(ins_scores) > 0.0
        error("insertion scores cannot be positive")
    end
    return normalize_log_differences(sub_scores,
                                     del_scores,
                                     ins_scores,
                                     state.score - max_score)
end


function estimate_point_probs(probs::EstErrorProbs)
    pos_probs = hcat(probs.sub, probs.del)
    no_point_error_prob = maximum(pos_probs, 2)
    # multiple by 0.5 to avoid double counting.
    no_ins_error_prob = 1.0 - 0.5 * sum(probs.ins, 2)
    result = 1.0 - broadcast(*, no_point_error_prob,
                             no_ins_error_prob[1:end-1],
                             no_ins_error_prob[2:end])
    return reshape(result, length(result))
end


function base_distribution(base::Char, lp, ilp)
    result = fill(lp - log10(3), 4)
    result[baseints[base]] = ilp
    return result
end


function posterior_error_probs(tlen::Int,
                               seqs::Vector{PString},
                               Amoves::Vector{BandedArray{Int}})
    # FIXME: incorporate scores
    # FIXME: account for indels
    probs = zeros(tlen, 4)
    for (s, Am) in zip(seqs, Amoves)
        moves = backtrace(Am)
        result = zeros(Int, tlen + 1)
        i, j = (1, 1)
        for move in moves
            (ii, jj) = offsets[Int(move)]
            i += ii
            j += jj
            s_i = i - 1
            t_j = j - 1
            if move == dp_match
                probs[t_j, 1:4] += base_distribution(s.seq[s_i],
                                                     s.error_log_p[s_i],
                                                     s.match_log_p[s_i])
            end
        end
    end
    probs = exp10(probs)
    probs = 1.0 - maximum(probs ./ sum(probs, 2), 2)
    return vec(probs)
end


function quiver2(seqstrings::Vector{String},
                 error_log_ps::Vector{Vector{Float64}},
                 scores::Scores;
                 consensus::String="",
                 reference::String="",
                 ref_scores::Scores=default_ref_scores,
                 ref_indel_penalty::Float64=-3.0,
                 min_ref_indel_score::Float64=-15.0,
                 do_alignment_proposals::Bool=true,
                 do_surgery_proposals::Bool=true,
                 trust_proposals::Bool=true,
                 fix_indels_stat::Bool=true,
                 skip_frame_correction::Bool=true,
                 indel_correction_only::Bool=true,
                 use_ref_for_qvs::Bool=false,
                 bandwidth::Int=10, bandwidth_mult::Int=2,
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

    if length(reference) > 0
        if ref_indel_penalty >= 0.0
            error("ref_indel_penalty must be < 0.0")
        end
        if min_ref_indel_score >= 0.0
            error("min_ref_indel_score must be < 0.0")
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

    sequences = PString[PString(s, p, bandwidth)
                        for (s, p) in zip(seqstrings, error_log_ps)]

    # will need to update after initial stage
    ref_error_rate = 1.0
    ref_error_log_p = fill(log10(ref_error_rate), length(reference))
    ref_pstring = PString(reference, ref_error_log_p, bandwidth)

    if max_iters < 1
        error("invalid max iters: $max_iters")
    end

    if batch < 0 || batch > length(sequences)
        batch = length(sequences)
    end
    base_batch = batch
    seqs = sequences
    if batch < length(sequences)
        indices = rand(1:length(sequences), batch)
        seqs = sequences[indices]
    end

    if verbose > 1
        println(STDERR, "computing initial alignments")
    end
    state = initial_state(consensus, seqs)
    recompute_Bs = !(do_surgery_proposals && trust_proposals)
    recompute!(state, seqs, scores,
               ref_pstring, ref_scores,
               bandwidth_mult, true, recompute_Bs, verbose,
               use_ref_for_qvs)
    empty_ref = length(reference) == 0

    if verbose > 1
        println(STDERR, "initial score: $(state.score)")
    end
    if verbose > 2
        println(STDERR, "  consensus: $(state.consensus)")
    end

    n_proposals = Vector{Int}[]
    consensus_lengths = Int[length(consensus)]
    consensus_stages = ["", "", ""]

    stage_iterations = zeros(Int, Int(typemax(Stage)) - 1)
    stage_times = Float64[0, 0, 0]
    tic()
    for i in 1:max_iters
        stage_iterations[Int(state.stage)] += 1
        old_consensus = state.consensus
        old_score = state.score
        if verbose > 1
            println(STDERR, "iteration $i : $(state.stage). score: $(state.score)")
        end

        penalties_increased = false
        candidates = get_candidate_proposals(state, seqs, scores, ref_pstring,
                                             ref_scores,
                                             do_alignment_proposals,
                                             do_surgery_proposals,
                                             trust_proposals,
                                             indel_correction_only)
        recompute_As = true
        if length(candidates) == 0
            if verbose > 1
                println(STDERR, "  no candidates found")
            end
            push!(n_proposals, zeros(Int, 3))
            consensus_stages[Int(state.stage)] = state.consensus
            stage_times[Int(state.stage)] = toq()
            tic()
            if state.stage == initial_stage
                if empty_ref
                    state.converged = true
                    break
                end
                state.stage = frame_correction_stage

                # fix distant single indels right away
                if fix_indels_stat
                    cons_errors = posterior_error_probs(length(state.consensus),
                                                      seqs, state.Amoves)
                    # ensure none are 0.0
                    cons_errors = [max(p, 1e-10) for p in cons_errors]
                    cons_pstring = PString(state.consensus, log10(cons_errors), bandwidth)
                    indel_proposals = single_indel_proposals(reference, cons_pstring, ref_scores)
                    if verbose > 1
                        println(STDERR, "  fixing $(length(indel_proposals)) single indels")
                    end
                    if verbose > 2
                        println(STDERR, indel_proposals)
                    end
                    state.consensus = apply_proposals(state.consensus, indel_proposals)
                    if skip_frame_correction
                        if verbose > 1
                            println(STDERR, "  skipping straight to refinement stage")
                        end
                        state.stage = refinement_stage
                    end
                end
                if state.stage == frame_correction_stage
                    # estimate reference error rate
                    # TODO: use consensus estimated error rate here too
                    edit_dist = levenshtein(reference, state.consensus)
                    ref_error_rate = edit_dist / max(length(reference), length(state.consensus))
                    # needs to be < 0.5, otherwise matches aren't rewarded at all
                    ref_error_rate = min(max(ref_error_rate, 1e-10), 0.5)
                    ref_error_log_p = fill(log10(ref_error_rate), length(reference))
                    ref_pstring = PString(reference, ref_error_log_p, bandwidth)
                end
            elseif state.stage == frame_correction_stage
                if !has_single_indels(state.consensus, ref_pstring, ref_scores)
                    consensus_ref = state.consensus
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
                println(STDERR, "  found $(length(candidates)) candidates")
            end
            chosen_cands = choose_candidates(candidates, min_dist)
            if verbose > 1
                println(STDERR, "  filtered to $(length(chosen_cands)) candidates")
            end
            if verbose > 2
                println(STDERR, "  $(chosen_cands)")
            end
            state.consensus = apply_proposals(old_consensus,
                                              Proposal[c.proposal
                                                       for c in chosen_cands])
            recompute_Bs = (!(do_surgery_proposals && trust_proposals) ||
                            state.stage == frame_correction_stage)
            recompute!(state, seqs, scores,
                       ref_pstring, ref_scores,
                       bandwidth_mult, true, recompute_Bs, verbose,
                       use_ref_for_qvs)
            # detect if a single proposal is better
            # note: this may not always be correct, because score_proposal() is not exact
            if length(chosen_cands) > 1 &&
                (state.score < chosen_cands[1].score ||
                 isapprox(state.score, chosen_cands[1].score))
                if verbose > 1
                    println(STDERR, "  rejecting multiple candidates in favor of best")
                end
                chosen_cands = CandProposal[chosen_cands[1]]
                if verbose > 2
                    println(STDERR, "  $(chosen_cands)")
                end
                state.consensus = apply_proposals(old_consensus,
                                                  Proposal[c.proposal
                                                           for c in chosen_cands])
            else
                # no need to recompute unless batch changes
                recompute_As = false
            end
            proposal_counts = [length(filter(c -> (typeof(c.proposal) == t),
                                             chosen_cands))
                               for t in [Substitution, Insertion, Deletion]]
            push!(n_proposals, proposal_counts)
        end
        push!(consensus_lengths, length(state.consensus))
        if batch < length(sequences)
            indices = rand(1:length(sequences), batch)
            seqs = sequences[indices]
            recompute_As = true
        end
        recompute_Bs = (!(do_surgery_proposals && trust_proposals) ||
                        state.stage == frame_correction_stage)
        recompute!(state, seqs, scores,
                   ref_pstring, ref_scores,
                   bandwidth_mult, recompute_As, recompute_Bs, verbose,
                   use_ref_for_qvs)
        if verbose > 1
            println(STDERR, "  score: $(state.score) ($(state.stage))")
        end
        if verbose > 2
            println(STDERR, "  consensus: $(state.consensus)")
        end
        if (state.score < old_score &&
            stage_iterations[Int(state.stage)] > 0 &&
            !penalties_increased &&
            batch == length(sequences))
            if verbose > 1
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
            recompute_Bs = (!(do_surgery_proposals && trust_proposals) ||
                            state.stage == frame_correction_stage)
            recompute!(state, seqs, scores,
                       ref_pstring, ref_scores,
                       bandwidth_mult, true, recompute_Bs, verbose,
                       use_ref_for_qvs)
            if verbose > 1
                println(STDERR, "  increased batch size to $batch. new score: $(state.score)")
            end
        end
    end
    state.stage = scoring_stage
    if verbose > 0
        println(STDERR, "done. converged: $(state.converged)")
    end
    push!(consensus_lengths, length(state.consensus))
    exceeded = sum(stage_iterations) >= max_iters

    info = Dict("converged" => state.converged,
                "stage_iterations" => stage_iterations,
                "exceeded_max_iterations" => exceeded,
                "ref_scores" => ref_scores,
                "consensus_stages" => consensus_stages,
                "n_proposals" => transpose(hcat(n_proposals...)),
                "consensus_lengths" => consensus_lengths,
                "ref_error_rate" => ref_error_rate,
                "stage_times" => stage_times,
                )
    # FIXME: recomputing for all sequences is costly, but using batch
    # is less accurate
    recompute!(state, seqs, scores, ref_pstring, ref_scores,
               bandwidth_mult, true, true, verbose,
               use_ref_for_qvs)
    info["error_probs"] = estimate_probs(state, seqs, scores,
                                         ref_pstring, ref_scores,
                                         use_ref_for_qvs)
    return state.consensus, info
end

function quiver2(sequences::Vector{String},
                 phreds::Vector{Vector{Int8}},
                 scores::Scores;
                 kwargs...)
    if any(minimum(p) < 0 for p in phreds)
        error("phred score cannot be negative")
    end
    error_log_ps = Util.phred_to_log_p(phreds)
    return quiver2(sequences, error_log_ps, scores; kwargs...)
end


"""
Alternate quiver2() using BioJulia types.

"""
function quiver2(sequences::Vector{DNASequence},
                 phreds::Vector{Vector{Int8}},
                 scores::Scores;
                 consensus::DNASequence=DNASequence(""),
                 reference::DNASequence=DNASequence(""),
                 kwargs...)
    new_reference = convert(String, reference)
    new_consensus = convert(String, consensus)
    new_sequences = String[convert(String, s) for s in sequences]
    (result, info) = quiver2(new_sequences, phreds, scores;
                             consensus=new_consensus,
                             reference=new_reference,
                             kwargs...)
    info["consensus_stages"] = DNASequence[DNASequence(s) for s in info["consensus_stages"]]
    return DNASequence(result), info
end

end
