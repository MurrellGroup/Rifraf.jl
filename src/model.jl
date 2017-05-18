@enum(Stage,
      STAGE_INIT=1,    # no reference; all proposals
      STAGE_FRAME=2,   # reference; indel proposals
      STAGE_REFINE=3,  # no reference; substitutions
      STAGE_SCORE=4)

function next_stage(s::Stage)
    return Stage(Int(s) + 1)
end

immutable EstimatedProbs
    sub::Array{Prob, 2}
    del::Array{Prob, 1}
    ins::Array{Prob, 2}
end

@with_kw type RifrafParams
    ref_scores::Scores = Scores(0.0, 0.0, 0.0, 0.0, 0.0)
    ref_indel_penalty::Score = -3.0
    min_ref_indel_score::Score = -15.0
    enabled_stages::Vector{Stage} = [STAGE_INIT,
                                     STAGE_FRAME,
                                     STAGE_REFINE,
                                     STAGE_SCORE]
    do_alignment_proposals::Bool = true
    do_surgery_proposals::Bool = false
    seed_indels::Bool = true
    indel_correction_only::Bool = true
    use_ref_for_qvs::Bool = false
    bandwidth::Int = (3*CODON_LENGTH)
    min_dist::Int = (5 * CODON_LENGTH)
    batch_size::Int = 10
    batch_threshold::Float64 = 0.05
    max_iters::Int = 100
    verbose::Int = 0
end

@with_kw type RifrafState
    score::Score = -Inf
    consensus::DNASeq
    ref_scores::Scores
    ref_error_rate::Float64 = -Inf
    reference::RifrafSequence
    batch_size::Int
    base_batch_size::Int
    sequences::Vector{RifrafSequence}
    batch_seqs::Vector{RifrafSequence} = RifrafSequence[]
    maxlen::Int
    As::Vector{BandedArray{Score}}
    Bs::Vector{BandedArray{Score}}
    Amoves::Vector{BandedArray{Trace}}
    A_ref::BandedArray{Score}
    B_ref::BandedArray{Score}
    Amoves_ref::BandedArray{Trace}
    realign_As::Bool = true
    realign_Bs::Bool = true
    penalties_increased::Bool = false
    stage::Stage = STAGE_INIT
    stage_iterations::Vector{Int} = zeros(Int, Int(typemax(Stage)))
    converged::Bool = false
end

@with_kw type RifrafResult
    consensus::DNASeq
    params::RifrafParams
    state::RifrafState
    consensus_stages::Vector{Vector{DNASeq}}
    error_probs::EstimatedProbs = EstimatedProbs(Array{Prob, 2}(),
                                                 Array{Prob, 1}(),
                                                 Array{Prob, 2}())
    aln_error_probs::Vector{Float64} = Float64[]
end

function equal_ranges(a_range::Tuple{Int,Int},
                      b_range::Tuple{Int,Int})
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

function seq_score_deletion(A::BandedArray{Score}, B::BandedArray{Score},
                            acol::Int, bcol::Int)
    Acol = sparsecol(A, acol)
    Bcol = sparsecol(B, bcol)
    (amin, amax), (bmin, bmax) = equal_ranges(row_range(A, acol),
                                              row_range(B, bcol))
    asub = view(Acol, amin:amax)
    bsub = view(Bcol, bmin:bmax)
    return summax(asub, bsub)
end

const BOFFSETS = Dict(Substitution => 2,
                      Insertion => 1,
                      Deletion => 2)

function score_nocodon(proposal::Proposal,
                       A::BandedArray{Score}, B::BandedArray{Score},
                       pseq::RifrafSequence,
                       newcols::Array{Score, 2})
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
    amin, amax = row_range(A, min(new_acol, ncols))
    for i in amin:amax
        seq_base = i > 1 ? pseq.seq[i-1] : DNA_Gap
        newcols[i, 1], _ = update(A, i, new_acol,
                                  seq_base, proposal.base, pseq;
                                  newcols=newcols, acol=acol)
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

function get_consensus_substring(proposal::Proposal, seq::DNASeq,
                                 n_after::Int)
    # next valid position in sequence after this proposal
    prefix = if typeof(proposal) in (Substitution, Insertion)
        DNASeq([proposal.base])
    else
        DNASeq()
    end
    pos = proposal.pos
    next_posn = pos + 1
    stop = min(next_posn + n_after - 1, length(seq))
    suffix = seq[next_posn:stop]
    return DNASeq(prefix, suffix)
end

function score_proposal(proposal::Proposal,
                        A::BandedArray{Score}, B::BandedArray{Score},
                        consensus::DNASeq,
                        pseq::RifrafSequence,
                        newcols::Array{Score, 2})
    if !do_codon_moves(pseq)
        return score_nocodon(proposal, A, B, pseq, newcols)
    end
    t = typeof(proposal)
    # last valid column of A
    acol_offset = t == Insertion ? 0 : -1
    acol = proposal.pos + acol_offset + 1

    # first column of B to use
    first_bcol = acol + BOFFSETS[t]
    # last column of B to use
    last_bcol = first_bcol + CODON_LENGTH - 1

    if t == Deletion && acol == (size(A)[2] - 1)
        # suffix deletions do not need recomputation
        return A[end, end - 1]
    end

    nrows, ncols = size(A)

    # if we'd go to or past the last column of B, just recompute the
    # rest of A, including A[end, end]
    just_a = last_bcol >= ncols

    # number of columns after new columns to also recompute.
    n_after = !just_a ? CODON_LENGTH : length(consensus) - proposal.pos

    # number of bases changed/inserted
    n_new_bases = (t == Deletion ? 0 : 1)
    if n_new_bases == 0 && n_after == 0
        error("no new columns need to be recomputed.")
    end

    n_new = n_new_bases + n_after
    sub_consensus = get_consensus_substring(proposal, consensus, n_after)
    # compute new columns
    for j in 1:n_new
        range_col = min(acol + j, ncols)
        amin, amax = row_range(A, range_col)
        for i in amin:amax
            seq_base = i > 1 ? pseq.seq[i-1] : DNA_Gap
            newcols[i, j], _ = update(A, i, acol + j,
                                      seq_base, sub_consensus[j], pseq;
                                      newcols=newcols, acol=acol)
        end
    end

    if just_a
        return newcols[nrows, n_new]
    end

    # add up results
    best_score = -Inf
    for j in 1:CODON_LENGTH
        new_j = n_new - CODON_LENGTH + j
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

function score_proposal(m::Proposal,
                        state::RifrafState,
                        newcols::Array{Score, 2},
                        use_ref::Bool)
    score = 0.0
    for si in 1:length(state.batch_seqs)
        score += score_proposal(m, state.As[si], state.Bs[si], state.consensus,
                                state.batch_seqs[si], newcols)
    end
    if use_ref
        score += score_proposal(m, state.A_ref, state.B_ref, state.consensus,
                                state.reference, newcols)
    end
    return score
end

function all_proposals(stage::Stage,
                       consensus::DNASeq,
                       indel_correction_only::Bool,
                       indel_seeds::Vector{Proposal}=Proposal[],
                       seed_neighborhood::Int=CODON_LENGTH)
    # FIXME: look into reducing allocations by not producing new
    # instances of every possible proposal
    len = length(consensus)
    ins_positions = Set{Int}()
    del_positions = Set{Int}()

    for p in indel_seeds
        if typeof(p) == Insertion
            for idx in max(p.pos - seed_neighborhood, 0):min(p.pos + seed_neighborhood, len)
                push!(ins_positions, idx)
            end
        else
            for idx in max(p.pos - seed_neighborhood, 1):min(p.pos + seed_neighborhood, len)
                push!(del_positions, idx)
            end
        end
    end
    do_subs = stage != STAGE_FRAME || !indel_correction_only
    do_indels = stage == STAGE_INIT ||
                stage == STAGE_FRAME ||
                stage == STAGE_SCORE
    no_seeds = length(indel_seeds) == 0
    function _it()
        if do_indels
            for base in BASES
                produce(Insertion(0, base))
            end
        end
        for j in 1:len
            if do_subs
                # substitutions
                for base in BASES
                    if consensus[j] != base
                        produce(Substitution(j, base))
                    end
                end
            end
            if do_indels
                # single indels
                if no_seeds || j in del_positions
                    produce(Deletion(j))
                end
                if no_seeds || j in ins_positions
                    for base in BASES
                        produce(Insertion(j, base))
                    end
                end
            end
        end
    end
    Task(_it)
end

function moves_to_proposals(moves::Vector{Trace},
                            consensus::DNASeq, seq::RifrafSequence)
    proposals = Proposal[]
    i, j = (0, 0)
    for move in moves
        i, j = offset_forward(move, i, j)

        score = seq.match_scores[max(i, 1)]
        next_score = seq.match_scores[min(i + 1, length(seq))]
        del_score = min(score, next_score)

        if move == TRACE_MATCH
            if seq.seq[i] != consensus[j]
                push!(proposals, Substitution(j, seq.seq[i]))
            end
        elseif move == TRACE_INSERT
            push!(proposals, Insertion(j, seq.seq[i]))
        elseif move == TRACE_DELETE
            push!(proposals, Deletion(j))
        end
    end
    return proposals
end

"""Only get proposals that appear in at least one alignment"""
function alignment_proposals(Amoves::Vector{BandedArray{Trace}},
                             consensus::DNASeq,
                             sequences::Vector{RifrafSequence},
                             do_indels::Bool)
    result = Set{Proposal}()
    for (Amoves, seq) in zip(Amoves, sequences)
        moves = backtrace(Amoves)
        for proposal in moves_to_proposals(moves, consensus, seq)
            if do_indels || typeof(proposal) == Substitution
                push!(result, proposal)
            end
        end
    end
    return collect(result)
end

"""Inserted bases surrounding a substitution"""
function best_surrounding_ins_bases(moves::Vector{Trace},
                                    seq::RifrafSequence,
                                    aln_idx::Int, seq_idx::Int)
    ins_bases = Dict{DNANucleotide, Tuple{Score, Score}}()
    # search before
    search_idx = aln_idx - 1
    seq_search_idx = seq_idx - 1
    while search_idx >= 1 && moves[search_idx] == TRACE_INSERT
        search_base = seq.seq[seq_search_idx]
        if !haskey(ins_bases, search_base)
            ins_bases[search_base] = (-Inf, -Inf)
        end
        match_score = seq.match_scores[seq_search_idx]
        best_match_score, _ = ins_bases[search_base]
        if match_score > best_match_score
            ins_score = seq.ins_scores[seq_search_idx]
            ins_bases[search_base] = (match_score, ins_score)
        end
        search_idx -= 1
        seq_search_idx -= 1
    end
    # search after
    search_idx = aln_idx + 1
    seq_search_idx = seq_idx + 1
    while search_idx <= length(moves) && moves[search_idx] == TRACE_INSERT
        search_base = seq.seq[seq_search_idx]
        if !haskey(ins_bases, search_base)
            ins_bases[search_base] = (-Inf, -Inf)
        end
        match_score = seq.match_scores[seq_search_idx]
        best_match_score, _ = ins_bases[search_base]
        if match_score > best_match_score
            ins_score = seq.ins_scores[seq_search_idx]
            ins_bases[search_base] = (match_score, ins_score)
        end
        search_idx += 1
        seq_search_idx += 1
    end
    return ins_bases
end

"""Deleted bases surrounding a substitution"""
function surrounding_del_bases(moves::Vector{Trace},
                               consensus::DNASeq,
                               seq::RifrafSequence, aln_idx::Int,
                               cons_idx::Int, seq_idx::Int)
    result = Dict{DNANucleotide, Tuple{Score, Score}}()
    del_idx = seq_idx + 1

    # search backwards
    old_deletion_score = seq.del_scores[del_idx-1]
    new_deletion_score = seq.del_scores[del_idx]
    search_idx = aln_idx - 1
    cons_search_idx = cons_idx - 1
    while search_idx >= 1 && moves[search_idx] == TRACE_DELETE
        search_base = consensus[cons_search_idx]
        result[search_base] = (old_deletion_score, new_deletion_score)
        search_idx -= 1
        cons_search_idx -= 1
    end

    # search forwards
    search_idx = aln_idx + 1
    cons_search_idx = cons_idx + 1
    old_deletion_score = seq.del_scores[del_idx]
    new_deletion_score = seq.del_scores[del_idx-1]
    while search_idx <= length(seq) && moves[search_idx] == TRACE_DELETE
        search_base = consensus[cons_search_idx]
        if haskey(result, search_base)
            prev_old, prev_new = result[search_base]
            if new_deletion_score > prev_new
                result[search_base] = (old_deletion_score, new_deletion_score)
            end
        else
            result[search_base] = (old_deletion_score, new_deletion_score)
        end
        search_idx += 1
        cons_search_idx += 1
    end
    return result
end

function seq_posn_scores(seq::RifrafSequence, seq_idx::Int)
    if seq_idx > 0
        return (seq.match_scores[seq_idx],
                log10(1. - exp10(seq.mismatch_scores[seq_idx])),
                seq.mismatch_scores[seq_idx],
                seq.ins_scores[seq_idx],
                seq.del_scores[seq_idx + 1])
    else
        return (-Inf,
                -Inf,
                -Inf,
                -Inf,
                seq.del_scores[seq_idx + 1])
    end
end

function update_deltas(sub_deltas::Array{Score, 2},
                       del_deltas::Array{Score, 1},
                       ins_deltas::Array{Score, 2},
                       consensus::DNASeq,
                       seq::RifrafSequence, moves::Vector{Trace},
                       do_subs::Bool, do_indels::Bool)
    # FIXME: modularize this function and test each part

    # local alignment fixes:
    # - insertions next to a matching base get pushed to end of match
    # - deletions get pushed to the end of poly-base runs
    # - substitutions that convert matches to mismatches try to swap
    #   with neighboring insertions and deletions
    # - insertions next to an eligible (mis)match get swapped

    # push insertions and deletions to the end, on the fly
    # insertions that match consensus should *keep* getting pushed
    # deletions in a poly-base run should also get pushed
    deletion_score = seq.del_scores[1]
    insertion_bases = Dict{DNANucleotide, Score}(DNA_A => deletion_score,
                                                 DNA_C => deletion_score,
                                                 DNA_G => deletion_score,
                                                 DNA_T => deletion_score)

    del_base = DNA_Gap
    del_idx = 0
    pushed_del_score = -Inf

    cons_idx = 0
    seq_idx = 0
    n_moves = length(moves)
    prev_del_score = -Inf
    for (aln_idx, move) in enumerate(moves)
        cbase = DNA_Gap
        sbase = DNA_Gap
        if move in (TRACE_MATCH, TRACE_INSERT)
            seq_idx += 1
            sbase = seq.seq[seq_idx]
        end
        if move in (TRACE_MATCH, TRACE_DELETE)
            cons_idx += 1
            cbase = consensus[cons_idx]
        end

        all_scores = seq_posn_scores(seq, seq_idx)
        (match_score, error_score, mm_score, ins_score, del_score) = all_scores

        if do_indels
            if del_base != DNA_Gap && del_base != cbase
                # handle pushed deletions
                del_deltas[del_idx] += pushed_del_score
                del_base = DNA_Gap
                del_idx = 0
                pushed_del_score = -Inf
            end
            if move != TRACE_INSERT
                # handle insertions before this position.
                for (ins_base, delta) in insertion_bases
                    if cbase == ins_base
                        # do nothing; keep pushing this insertion
                    else
                        # update `insertion_bases` with candidate
                        # insertion after this position
                        # this move is a mismatch, consider
                        # swapping with prev and next insertions

                        # FIXME: this is causing duplicate insertions
                        # before and after the same base. Detect this
                        # case and account for it.
                        new_delta = del_score
                        if move == TRACE_MATCH && sbase == ins_base && cbase != sbase
                            prev_swap_delta = (del_score + match_score) - mm_score
                            next_swap_delta = (prev_del_score + match_score) - mm_score
                            delta = max(delta, prev_swap_delta)
                            new_delta = max(new_delta, next_swap_delta)
                        end
                        ins_deltas[cons_idx, BASEINTS[ins_base]] += delta
                        insertion_bases[ins_base] = new_delta
                    end
                end
            end
        end
        old_match_score = (sbase == cbase ? match_score : mm_score)
        if move == TRACE_INSERT && do_indels
            # consider insertion proposals.
            # push all insertion proposals to the end
            for (baseint, new_base) in enumerate(BASES)
                # try all insertions, converting TRACE_INSERT to (mis)match
                sbase = seq.seq[seq_idx]
                new_score = (sbase == new_base) ? match_score : mm_score
                insertion_bases[new_base] = max(insertion_bases[new_base],
                                                new_score - ins_score)
            end
        elseif move == TRACE_MATCH
            # consider all substitution proposals
            if do_subs
                # TODO: do this on the fly
                best_ins = best_surrounding_ins_bases(moves, seq, aln_idx, seq_idx)
                del_bases = surrounding_del_bases(moves, consensus, seq,
                                                  aln_idx, cons_idx, seq_idx)

                for (baseint, new_base) in enumerate(BASES)
                    if new_base == cbase
                        continue
                    end
                    delta = -Inf

                    if haskey(best_ins, new_base)
                        # do substitution and align with a different insertion.
                        # ins/match -> match/ins
                        other_match_score, other_ins_score = best_ins[new_base]
                        old_score = other_ins_score + match_score
                        new_score = other_match_score + ins_score
                        delta = max(delta, new_score - old_score)
                    end
                    if sbase == new_base
                        # proposal would help
                        delta = max(delta, (match_score - mm_score))
                    elseif sbase == cbase
                        if haskey(del_bases, sbase)
                            # del/match -> match/del
                            # do substitution and delete this base, moving seq
                            # base to different position
                            old_del_score, new_del_score = del_bases[sbase]
                            old_score = old_match_score + old_del_score
                            new_score = match_score + new_del_score
                            delta = max(delta, new_score - old_score)
                        end
                        # proposal would hurt
                        delta = max(delta, (mm_score - match_score))
                    else
                        # proposal would still be a mismatch. do nothing.
                    end
                    if delta > -Inf
                        sub_deltas[cons_idx, baseint] += delta
                    end
                end
            end
            # consider deleting the consensus base.
            # This would convert (mis)match to insertion
            del_base = cbase
            del_idx = cons_idx
            pushed_del_score = max(pushed_del_score,
                                   ins_score - old_match_score)
        elseif do_indels && move == TRACE_DELETE
            # no reason to consider substitutions. delta = 0.
            # consider deletion proposal. would convert deletion
            # to no-op.
            del_base = cbase
            del_idx = cons_idx
            pushed_del_score = max(pushed_del_score,
                                   0 - del_score)
        end
        prev_del_score = del_score
     end
    if do_indels
        # handle deletion at end
        if del_base != DNA_Gap
            del_deltas[end] += pushed_del_score
        end

        # handle insertions at the end
        for (base, delta) in insertion_bases
            ins_deltas[end, BASEINTS[base]] += delta
        end
    end
end

"""Use model surgery heuristic to get proposals with positive score deltas."""
function surgery_candidates(state::RifrafState, do_indels::Bool)
    sub_deltas = zeros(Score, (length(state.consensus), 4))
    del_deltas = zeros(Score, (length(state.consensus)))
    ins_deltas = zeros(Score, (length(state.consensus) + 1, 4))

    do_subs = true
    for (seq, Amoves) in zip(state.batch_seqs, state.Amoves)
        moves = backtrace(Amoves)
        update_deltas(sub_deltas, del_deltas, ins_deltas,
                      state.consensus, seq, moves,
                      do_subs, do_indels)
    end
    # only return proposals with positive deltas
    result = ScoredProposal[]
    nrows, ncols = size(sub_deltas)
    if do_subs
        for i in 1:nrows, j in 1:ncols
            if sub_deltas[i, j] > 0.0
                push!(result, ScoredProposal(Substitution(i, BASES[j]),
                                           state.score + sub_deltas[i, j]))
            end
        end
    end
    if do_indels
        for i in 1:nrows
            if del_deltas[i] > 0.0
                push!(result, ScoredProposal(Deletion(i),
                                           state.score + del_deltas[i]))
            end
        end
        for i in 1:(nrows+1), j in 1:ncols
            if ins_deltas[i, j] > 0.0
                push!(result, ScoredProposal(Insertion(i - 1, BASES[j]),
                                           state.score + ins_deltas[i, j]))
            end
        end
    end
    return result
end

function get_candidates(state::RifrafState, params::RifrafParams;
                        indel_seeds::Vector{Proposal}=Proposal[])
    if params.do_alignment_proposals && params.do_surgery_proposals
        # TODO: allow both and merge
        error("cannot do both alignment and surgery proposals")
    end
    use_ref = (state.stage == STAGE_FRAME)

    maxlen = maximum(length(s) for s in state.batch_seqs)
    nrows = max(maxlen, length(state.reference)) + 1
    newcols = zeros(Score, (nrows, CODON_LENGTH + 1))

    if (state.stage == STAGE_INIT ||
        state.stage == STAGE_REFINE) && params.do_surgery_proposals
        do_indels = state.stage == STAGE_INIT
        return surgery_candidates(state, do_indels)
    end

    # multi-line ternary
    proposals = if (state.stage == STAGE_INIT ||
        state.stage == STAGE_REFINE) && params.do_alignment_proposals
        do_indels = state.stage == STAGE_INIT
        alignment_proposals(state.Amoves, state.consensus,
                            state.batch_seqs, do_indels)
    else
        all_proposals(state.stage, state.consensus,
                      params.indel_correction_only, indel_seeds)
    end

    candidates = ScoredProposal[]
    for p in proposals
        score = score_proposal(p, state, newcols, use_ref)
        if score > state.score
            push!(candidates, ScoredProposal(p, score))
        end
    end
    return candidates
end

function base_consensus(d::Dict{DNANucleotide, Score})
    return minimum((v, k) for (k, v) in d)[2]
end

function has_single_indels(consensus::DNASeq,
                           reference::RifrafSequence)
    has_right_length = length(consensus) % CODON_LENGTH == 0
    moves = align_moves(consensus, reference)
    result = TRACE_INSERT in moves || TRACE_DELETE in moves
    if !result && !has_right_length
        error("consensus length is not a multiple of three")
    end
    return result
end

function single_indel_proposals(reference::DNASeq,
                                consensus::RifrafSequence)
    moves = align_moves(reference, consensus,
                        skew_matches=true)
    results = Proposal[]
    cons_idx = 0
    ref_idx = 0
    for move in moves
        if move == TRACE_MATCH
            cons_idx += 1
            ref_idx += 1
        elseif move == TRACE_INSERT
            cons_idx += 1
            push!(results, Deletion(cons_idx))
        elseif move == TRACE_DELETE
            ref_idx += 1
            push!(results, Insertion(cons_idx, reference[ref_idx]))
        elseif move == TRACE_CODON_INSERT
            cons_idx += 3
        elseif move == TRACE_CODON_DELETE
            ref_idx += 3
        end
    end
    return results
end

function initial_state(consensus::DNASeq,
                       sequences::Vector{RifrafSequence},
                       reference::DNASeq,
                       params::RifrafParams)
    # batch size calculations
    batch_size = params.batch_size
    if batch_size < 0 || batch_size > length(sequences)
        batch_size = length(sequences)
    end
    base_batch_size = batch_size

    if length(consensus) == 0
        # choose highest-quality sequence as initial consensus
        idx = indmax([logsumexp10(s.match_scores) for s in sequences])
        consensus = sequences[idx].seq
    end

    maxlen = maximum(map(length, sequences))
    rlen = length(reference)

    # initialize sequence banded arrays
    # use longest sequence to ensure data reallocation is rarely necessary
    seq_shape = (maxlen + 1, maxlen + 1)
    As = BandedArray{Score}[BandedArray(Score, seq_shape, params.bandwidth,
                                        default=-Inf)
                            for _ in 1:batch_size]
    Bs = BandedArray{Score}[BandedArray(Score, seq_shape, params.bandwidth,
                                        default=-Inf)
                            for _ in 1:batch_size]
    Amoves = BandedArray{Trace}[BandedArray(Trace, seq_shape, params.bandwidth)
                                for _ in 1:batch_size]


    # initialize reference alignment arrays
    ref_shape = (maxlen + 1, rlen + 1)
    ref_bandwidth = params.bandwidth
    if length(reference) == 0
        ref_shape = (1, 1)
        ref_bandwidth = 1
    end
    A_ref = BandedArray(Score, ref_shape, ref_bandwidth, default=-Inf)
    B_ref = BandedArray(Score, ref_shape, ref_bandwidth, default=-Inf)
    Amoves_ref = BandedArray(Trace, ref_shape, ref_bandwidth)

    # just a placeholder for now
    ref_error_log_p = fill(0.0, length(reference))
    refseq = RifrafSequence(reference, ref_error_log_p,
                            params.bandwidth, params.ref_scores)

    stage = minimum(params.enabled_stages)
    return RifrafState(consensus=consensus,
                       reference=refseq,
                       ref_scores=params.ref_scores,
                       batch_size=batch_size,
                       base_batch_size=base_batch_size,
                       maxlen=maxlen,
                       sequences=sequences,
                       As=As, Bs=Bs, Amoves=Amoves,
                       A_ref=A_ref, B_ref=B_ref, Amoves_ref=Amoves_ref,
                       stage=stage)
end

function use_ref(ref, stage, use_ref_for_qvs)
    if length(ref) == 0
        return false
    end
    if stage == STAGE_FRAME
        return true
    end
    if stage == STAGE_SCORE && use_ref_for_qvs
        return true
    end
    return false
end

function rescore!(state::RifrafState, use_ref_for_qvs::Bool)
    state.score = sum(A[end, end] for A in state.As)
    if use_ref(state.reference, state.stage, use_ref_for_qvs)
        state.score += state.A_ref[end, end]
    end
end

"""Do forward and backward alignments, and compute scores"""
function realign!(state::RifrafState, use_ref_for_qvs::Bool)
    seqs = state.batch_seqs
    shape = (state.maxlen + 1, state.maxlen + 1)
    if state.realign_As
        # add new alignments, in case batch size increased
        for i in (length(state.As) + 1):length(seqs)
            push!(state.As, BandedArray(Score, shape, seqs[i].bandwidth))
            push!(state.Amoves, BandedArray(Trace, shape, seqs[i].bandwidth))
        end
        for (s, A, Amoves) in zip(seqs, state.As, state.Amoves)
            forward_moves!(state.consensus, s, A, Amoves)
        end
        if use_ref(state.reference, state.stage, use_ref_for_qvs)
            forward_moves!(state.consensus, state.reference,
                           state.A_ref, state.Amoves_ref)
        end
    end
    if state.realign_Bs
        # add new alignments, in case batch size increased
        for i in (length(state.Bs) + 1):length(seqs)
            push!(state.Bs, BandedArray(Score, shape, seqs[i].bandwidth))
        end
        for (s, B) in zip(seqs, state.Bs)
            backward!(state.consensus, s, B)
        end
        if use_ref(state.reference, state.stage, use_ref_for_qvs)
            backward!(state.consensus, state.reference, state.B_ref)
        end
    end
end

function realign_rescore!(state::RifrafState, use_ref_for_qvs::Bool)
    realign!(state, use_ref_for_qvs)
    rescore!(state, use_ref_for_qvs)
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
    return EstimatedProbs(sub_probs, del_probs, ins_probs)
end

function estimate_probs(state::RifrafState,
                        use_ref_for_qvs::Bool)
    sequences = state.batch_seqs
    reference = state.reference
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
    newcols = zeros(Score, (nrows, CODON_LENGTH + 1))

    use_ref = (length(reference) > 0) && use_ref_for_qvs
    # do not want to penalize indels too harshly, otherwise consensus
    # appears too likely.
    # TODO: how to set scores correctly for estimating qvs?
    # qv_ref_scores = Scores(ErrorModel(1.0, 1.0, 1.0, 0.0, 0.0))
    for m in all_proposals(STAGE_SCORE, state.consensus, false)
        score = score_proposal(m, state, newcols, use_ref)
        if typeof(m) == Substitution
            sub_scores[m.pos, BASEINTS[m.base]] = score
        elseif typeof(m) == Deletion
            del_scores[m.pos] = score
        elseif typeof(m) == Insertion
            ins_scores[m.pos + 1, BASEINTS[m.base]] = score
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

function estimate_point_probs(probs::EstimatedProbs)
    pos_probs = hcat(probs.sub, probs.del)
    no_point_error_prob = maximum(pos_probs, 2)
    # multiple by 0.5 to avoid double counting.
    no_ins_error_prob = 1.0 - 0.5 * sum(probs.ins, 2)
    result = 1.0 - broadcast(*, no_point_error_prob,
                             no_ins_error_prob[1:end-1],
                             no_ins_error_prob[2:end])
    return reshape(result, length(result))
end

function base_distribution(base::DNANucleotide, ilp)
    lp = log10(1.0 - exp10(ilp))
    result = fill(lp - log10(3), 4)
    result[BASEINTS[base]] = ilp
    return result
end

"""Assign error probabilities to consensus bases.

Compute posterior error probability from all sequence bases
that aligned to the consensus base.

"""
function alignment_error_probs(tlen::Int,
                               seqs::Vector{RifrafSequence},
                               Amoves::Vector{BandedArray{Trace}})
    # FIXME: incorporate scores
    # FIXME: account for indels
    probs = zeros(tlen, 4)
    for (s, Am) in zip(seqs, Amoves)
        moves = backtrace(Am)
        result = zeros(Int, tlen + 1)
        i, j = (1, 1)
        for move in moves
            i, j = offset_forward(move, i, j)
            s_i = i - 1
            t_j = j - 1
            if move == TRACE_MATCH
                probs[t_j, 1:4] += base_distribution(s.seq[s_i],
                                                     s.match_scores[s_i])
            end
        end
    end
    probs = exp10(probs)
    probs = 1.0 - maximum(probs ./ sum(probs, 2), 2)
    return vec(probs)
end

function check_args(scores, reference, params)
    if (scores.mismatch >= 0.0 || scores.mismatch == -Inf ||
        scores.insertion >= 0.0 || scores.insertion == -Inf ||
        scores.deletion >= 0.0 || scores.deletion == -Inf)
        error("scores must be between -Inf and 0.0")
    end
    if scores.codon_insertion > -Inf || scores.codon_deletion > -Inf
        error("error model cannot allow codon indels")
    end

    if length(reference) > 0
        if params.ref_indel_penalty >= 0.0
            error("ref_indel_penalty must be < 0.0")
        end
        if params.min_ref_indel_score >= 0.0
            error("min_ref_indel_score must be < 0.0")
        end
        if (params.ref_scores.mismatch >= 0.0 ||
            params.ref_scores.insertion >= 0.0 ||
            params.ref_scores.deletion >= 0.0 ||
            params.ref_scores.codon_insertion >= 0.0 ||
            params.ref_scores.codon_deletion >= 0.0)
            error("ref scores cannot be >= 0")
        end
        if (params.ref_scores.mismatch == -Inf ||
            params.ref_scores.insertion == -Inf ||
            params.ref_scores.deletion == -Inf ||
            params.ref_scores.codon_insertion == -Inf ||
            params.ref_scores.codon_deletion == -Inf)
            error("ref scores cannot be -Inf")
        end
        if (params.ref_scores.insertion < params.min_ref_indel_score ||
            params.ref_scores.deletion < params.min_ref_indel_score)
            error("ref indel scores are less than specified minimum")
        end
    end

    if length(params.enabled_stages) == 0
        error("no stages enabled")
    end
    if params.max_iters < 1
        error("invalid max iters: $(params.max_iters)")
    end
end

function handle_candidates!(candidates::Vector{ScoredProposal},
                            state::RifrafState,
                            params::RifrafParams)
    old_consensus = state.consensus
    chosen_cands = choose_candidates(candidates, params.min_dist)
    if params.verbose > 1
        println(STDERR, "  found $(length(candidates)) candidates")
    end
    if params.verbose > 1
        println(STDERR, "  filtered to $(length(chosen_cands)) candidates")
    end
    state.consensus = apply_proposals(old_consensus,
                                      Proposal[c.proposal
                                               for c in chosen_cands])
    state.realign_As = true
    state.realign_Bs = false
    realign_rescore!(state, params.use_ref_for_qvs)
    # detect if a single proposal is better
    # note: this may not always be correct, because score_proposal() is not exact
    if length(chosen_cands) > 1 &&
        (state.score < chosen_cands[1].score ||
         isapprox(state.score, chosen_cands[1].score))
        if params.verbose > 1
            println(STDERR, "  rejecting multiple candidates in favor of best")
        end
        chosen_cands = ScoredProposal[chosen_cands[1]]
        state.consensus = apply_proposals(old_consensus,
                                          Proposal[c.proposal
                                                   for c in chosen_cands])
    else
        # no need to realign As unless batch changes
        state.realign_As = false
    end
    state.realign_Bs = true
    proposal_counts = [length(filter(c -> (typeof(c.proposal) == t),
                                     chosen_cands))
                       for t in [Substitution, Insertion, Deletion]]
end

function finish_stage!(state::RifrafState,
                       params::RifrafParams)
    if params.verbose > 1
        println(STDERR, "  no candidates found. finishing $(state.stage).")
    end
    if state.stage == STAGE_INIT
        if length(state.reference) == 0 || !(STAGE_FRAME in params.enabled_stages)
            state.converged = true
        else
            state.stage = STAGE_FRAME
            # estimate reference error rate
            # TODO: use consensus estimated error rate here too
            edit_dist = levenshtein(convert(String, state.reference.seq),
                                    convert(String, state.consensus))
            ref_error_rate = edit_dist / max(length(state.reference),
                                             length(state.consensus))
            # needs to be < 0.5, otherwise matches aren't rewarded at all
            state.ref_error_rate = min(max(ref_error_rate, 1e-10), 0.5)
            ref_error_log_p = fill(log10(state.ref_error_rate), length(state.reference))
            state.reference = RifrafSequence(state.reference.seq, ref_error_log_p,
                                             params.bandwidth, state.ref_scores)
        end
    elseif state.stage == STAGE_FRAME
        if !has_single_indels(state.consensus, state.reference)
            consensus_ref = state.consensus
            state.stage = STAGE_REFINE
        elseif (state.ref_scores.insertion == params.min_ref_indel_score ||
                state.ref_scores.deletion == params.min_ref_indel_score)
            if params.verbose > 1
                println(STDERR, "  NOTE: alignment had single indels but indel scores already minimized.")
            end
            state.stage = STAGE_REFINE
        else
            state.penalties_increased = true
            # TODO: this is not probabilistically correct
            state.ref_scores = Scores(state.ref_scores.mismatch,
                                      max(state.ref_scores.insertion + params.ref_indel_penalty,
                                          params.min_ref_indel_score),
                                      max(state.ref_scores.deletion + params.ref_indel_penalty,
                                          params.min_ref_indel_score),
                                      state.ref_scores.codon_insertion,
                                      state.ref_scores.codon_deletion)
            state.reference = RifrafSequence(state.reference, state.ref_scores)
            if params.verbose > 1
                println(STDERR, "  NOTE: alignment to reference had single indels. increasing penalty.")
            end
        end
    elseif state.stage == STAGE_REFINE
        state.converged = true
    else
        error("invalid stage: $(state.stage)")
    end
end

function get_indel_seeds(state::RifrafState, params::RifrafParams)
    indel_proposals = Proposal[]
    if state.stage == STAGE_FRAME && params.seed_indels
        cons_errors = alignment_error_probs(length(state.consensus),
                                            state.batch_seqs, state.Amoves)
        # ensure none are 0.0
        cons_errors = [max(p, 1e-10) for p in cons_errors]
        # flip insertions and deletions
        cons_scores = Scores(params.ref_scores.mismatch,
                             params.ref_scores.deletion,
                             params.ref_scores.insertion,
                             params.ref_scores.codon_deletion,
                             params.ref_scores.codon_insertion)

        cons = RifrafSequence(state.consensus, log10(cons_errors),
                                      params.bandwidth, cons_scores)
        indel_proposals = single_indel_proposals(state.reference.seq, cons)
    end
    return indel_proposals
end

"""Samble a subset of `data`"""
function resample(data, batch_size)
    if batch_size < length(data)
        indices = rand(1:length(data), batch_size)
        return data[indices]
    else
        return data
    end
end

function resample!(state::RifrafState)
    state.batch_seqs = resample(state.sequences, state.batch_size)
    if state.batch_size < length(state.sequences)
        state.realign_As = true
    end
end

"""Check that score did not decrease.

If it did dicrease too much, increase batch size, resample, and
realign again.

"""
function check_score!(state::RifrafState, params::RifrafParams,
                      old_score::Float64)
    if (state.score < old_score &&
        !state.penalties_increased &&
        state.batch_size == length(state.sequences) &&
        state.stage_iterations[Int(state.stage)] > 1)
        if params.verbose > 1
            println(STDERR, "  WARNING: not using batches, but score decreased.")
        end
    end
    if ((state.score - old_score) / old_score > params.batch_threshold &&
        !state.penalties_increased &&
        state.batch_size < length(state.sequences) &&
        state.stage_iterations[Int(state.stage)] > 1)
        # TODO: could speed this up by adding new sequences and
        # realigning only the new ones
        state.batch_size = min(state.batch_size + state.base_batch_size,
                               length(state.sequences))
        resample!(state)
        state.realign_As = true
        state.realign_Bs = true
        realign_rescore!(state, params.use_ref_for_qvs)
        if params.verbose > 1
            println(STDERR, "  NOTE: increased batch size to $(batch_size). new score: $(state.score)")
        end
    end

end

function rifraf(dnaseqs::Vector{DNASeq},
                error_log_ps::Vector{Vector{LogProb}},
                scores::Scores;
                consensus::DNASeq=DNASeq(),
                reference::DNASeq=DNASeq(),
                params::RifrafParams=RifrafParams())
    check_args(scores, reference, params)
    enabled_stages = Set(params.enabled_stages)
    sequences = RifrafSequence[RifrafSequence(s, p, params.bandwidth, scores)
                        for (s, p) in zip(dnaseqs, error_log_ps)]

    state = initial_state(consensus, sequences, reference, params)

    # keep track of every iteration's consensus
    consensus_stages = [[] for _ in 1:(Int(typemax(Stage)) - 1)]

    state.realign_As = true
    state.realign_Bs = true
    old_score = -Inf

    for i in 1:params.max_iters
        # skip to next valid stage
        while state.stage < STAGE_SCORE && !(state.stage in enabled_stages)
            state.stage = next_stage(state.stage)
        end
        if state.stage == STAGE_SCORE
            break
        end
        state.stage_iterations[Int(state.stage)] += 1
        push!(consensus_stages[Int(state.stage)], state.consensus)

        if params.verbose > 1
            println(STDERR, "iteration $i : $(state.stage)")
        end
        if params.verbose > 2
            println(STDERR, "  consensus: $(state.consensus)")
        elseif params.verbose > 1
            println(STDERR, "  consensus length: $(length(state.consensus))")
        end

        # resample sequences, realign, and possibly adjust batch size
        resample!(state)
        state.realign_Bs = true
        realign_rescore!(state, params.use_ref_for_qvs)
        check_score!(state, params, old_score)
        old_score = state.score
        if params.verbose > 1
            println(STDERR, "  score: $(state.score)")
        end

        # get candidate changes to consensus
        state.realign_As = false
        state.penalties_increased = false
        indel_seeds = get_indel_seeds(state, params)
        candidates = get_candidates(state, params; indel_seeds=indel_seeds)

        # handle candidates
        state.realign_As = true
        if length(candidates) > 0
            handle_candidates!(candidates, state, params)
        else
            finish_stage!(state, params)
        end
        if state.converged
            break
        end
    end
    state.stage = STAGE_SCORE
    if params.verbose > 0
        println(STDERR, "done. converged: $(state.converged)")
    end

    result = RifrafResult(consensus=state.consensus,
                          params=params,
                          state=state,
                          consensus_stages=consensus_stages)

    if STAGE_SCORE in enabled_stages
        # FIXME: recomputing for all sequences is costly, but using batch
        # is less accurate
        state.realign_As = true
        state.realign_Bs = true
        realign_rescore!(state, params.use_ref_for_qvs)
        result.error_probs = estimate_probs(state, params.use_ref_for_qvs)
        result.aln_error_probs = alignment_error_probs(length(state.consensus),
                                                     state.batch_seqs, state.Amoves)
    end
    return result
end

function rifraf(sequences::Vector{DNASeq},
                phreds::Vector{Vector{Phred}},
                scores::Scores;
                kwargs...)
    if any(minimum(p) < 0 for p in phreds)
        error("phred score cannot be negative")
    end
    error_log_ps = phred_to_log_p(phreds)
    return rifraf(sequences, error_log_ps, scores; kwargs...)
end
