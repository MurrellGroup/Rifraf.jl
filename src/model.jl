@enum(Stage,
      STAGE_INIT = 1,    # no reference; all proposals
      STAGE_FRAME = 2,   # reference; indel proposals
      STAGE_REFINE = 3,  # no reference; substitutions
      STAGE_SCORE = 4)

function next_stage(s::Stage)
    return Stage(Int(s) + 1)
end

"""
Estimated error probabilities of a consensus sequence.

# Fields
- `sub::Array{Prob,2}`: Probabilities for all four bases.
- `del::Array{Prob,1}`: Probabilities for a deletion
- `ins::Array{Prob,2}`: Probabilities for all insertions.

"""
struct EstimatedProbs
    sub::Array{Prob,2}
    del::Array{Prob,1}
    ins::Array{Prob,2}
end

"""
The parameters for a RIFRAF run.

# Fields
- `scores::Scores = Scores(ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))`

- `ref_scores::Scores = Scores(ErrorModel(10.0, 1e-1, 1e-1, 1.0, 1.0))`

- `ref_indel_mult::Score = 3.0`: multiplier for single indel penalties
  in alignment with the reference

- `max_ref_indel_mults::Int = 5`: maximum multiplier increases for
  single indel penalty

- `ref_error_mult::Float64 = 1.0`: multiplier for estimated reference
  error rate.

- `do_init::Bool = true`: enable initialization stage

- `do_frame::Bool = true`: enable frame correction stage

- `do_refine::Bool = true`: enable refinement stage

- `do_score::Bool = false`: enable scoring stage

- `do_alignment_proposals::Bool = true`: only propose changes that
  occur in pairwise alignments

- `seed_indels::Bool = true`: seed indel locations from the alignment
  to reference

- `indel_correction_only::Bool = true`: only propose indels during
  frame correction stage

- `use_ref_for_qvs::Bool = false`: use reference alignment when
  estimating quality scores

- `bandwidth::Int = (3 * CODON_LENGTH)`: alignment bandwidth

- `bandwidth_pvalue::Float64 = 0.1`: p-value for increasing bandwidth

- `min_dist::Int = (5 * CODON_LENGTH)`: distance between accepted
  candidate proposals

- `batch_fixed::Bool = true`: use top sequences for initial stage and
  frame correction

- `batch_fixed_size::Int = 5`: size of fixed batch

- `batch_size::Int = 20`: batch size; if <= 1, no batching is used

- `batch_randomness::Float64 = 0.9`: batch randomness
  -  `0`: top n get picked
  -  `0.5`: weight according to estimated errors
  -  `1`: completely random

- `batch_mult::Float64 = 0.7`: multiplier to reduce batch randomness

- `batch_threshold::Float64 = 0.1`: score threshold for increasing
  batch size

- `max_iters::Int = 100`: maximum total iterations across all stages
  before giving up

- `verbose::Int = 0`: verbosity level
  - `0`: nothing
  - `1`: print iteration and score
  - `2`: also print step within each iteration
  - `3`: also print full consensus sequence

"""
@with_kw mutable struct RifrafParams
    scores::Scores = Scores(ErrorModel(1.0, 2.0, 2.0, 0.0, 0.0))
    ref_scores::Scores = Scores(ErrorModel(10.0, 1e-1, 1e-1, 1.0, 1.0))

    # multiplier for single indel penalty
    ref_indel_mult::Score = 3

    # limit for single indel penalty
    max_ref_indel_mults::Score = 5

    # multiplier for estimated reference error rate
    ref_error_mult::Float64 = 1.0

    # enable or disable various stages of the run
    do_init::Bool = true
    do_frame::Bool = true
    do_refine::Bool = true
    do_score::Bool = false

    # only propose changes that occur in pairwise alignments
    do_alignment_proposals::Bool = true

    # seed indel locations from the alignment to reference
    seed_indels::Bool = true

    # only propose indels during frame correction stage
    indel_correction_only::Bool = true

    # use reference alignment when estimating quality scores
    use_ref_for_qvs::Bool = false

    # alignment bandwidth
    bandwidth::Int = (3 * CODON_LENGTH)

    # p-value for increasing bandwidth
    bandwidth_pvalue::Float64 = 0.1

    # distance between accepted candidate proposals
    min_dist::Int = (5 * CODON_LENGTH)

    # use top sequences for initial stage and frame correction
    batch_fixed::Bool = true
    batch_fixed_size::Int = 5

    # if <= 1, no batching is used
    batch_size::Int = 20

    # 0: top n get picked
    # 0.5: weight according to estimated errors
    # 1: completely random
    batch_randomness::Float64 = 0.9

    # multiplier to reduce batch randomness
    batch_mult::Float64 = 0.7

    # score threshold for increasing batch size
    batch_threshold::Float64 = 0.1

    # maximum total iterations across all stages before giving up
    max_iters::Int = 100

    # verbosity level
    # 0: nothing
    # 1: print iteration and score
    # 2: also print step within each iteration
    # 3: also print full consensus sequence
    verbose::Int = 0
end

"""The mutable state that keeps track of a RIFRAF run."""
@with_kw mutable struct RifrafState
    score::Score = -Inf
    consensus::DNASeq
    ref_scores::Scores
    ref_error_rate::Float64 = -Inf
    n_ref_indel_mults::Int = 0
    reference::RifrafSequence
    batch_fixed_size::Int
    batch_size::Int
    base_batch_size::Int
    batch_randomness::Float64 = 0.9
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

"""
    RifrafResult()

The result of a RIFRAF run.

# Fields
- `consensus::DNASeq`: the consensus found by RIFRAF.

- `params::RifrafParams`: the parameters used for this run.

- `state::RifrafState`: the final state of the run.

- `consensus_stages::Vector{Vector{DNASeq}}`:

- `error_probs::EstimatedProbs`: estimated per-base probabilities for
  each position. Only available if `params.do_score` is `true`.

- `aln_error_probs::Vector{Float64}`: combined per-base error
  probabilities. Only available if `params.do_score` is `true`.

"""
@with_kw mutable struct RifrafResult
    consensus::DNASeq
    params::RifrafParams
    state::RifrafState
    consensus_stages::Vector{Vector{DNASeq}}
    error_probs::EstimatedProbs = EstimatedProbs(Array{Prob,2}(0, 0),
                                                 Array{Prob,1}(0),
                                                 Array{Prob,2}(0, 0))
    aln_error_probs::Vector{Float64} = Float64[]
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
                       newcols::Array{Score,2})
    if size(A)[1] != length(pseq) + 1
        error("wrong size array")
    end
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
        seq_base = i > 1 ? pseq.seq[i - 1] : DNA_Gap
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
                        newcols::Array{Score,2})
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
            seq_base = i > 1 ? pseq.seq[i - 1] : DNA_Gap
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
                        newcols::Array{Score,2},
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
    results = Proposal[]
    if do_indels
        for base in BASES
            push!(results, Insertion(0, base))
        end
    end
    for j in 1:len
        if do_subs
            # substitutions
            for base in BASES
                if consensus[j] != base
                    push!(results, Substitution(j, base))
                end
            end
        end
        if do_indels
            # single indels
            if no_seeds || j in del_positions
                push!(results, Deletion(j))
            end
            if no_seeds || j in ins_positions
                for base in BASES
                    push!(results, Insertion(j, base))
                end
            end
        end
    end
    results
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

function get_candidates(state::RifrafState, params::RifrafParams;
                        indel_seeds::Vector{Proposal}=Proposal[])
    use_ref = (state.stage == STAGE_FRAME)

    maxlen = maximum(length(s) for s in state.batch_seqs)
    nrows = max(maxlen, length(state.reference)) + 1
    newcols = zeros(Score, (nrows, CODON_LENGTH + 1))

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

function base_consensus(d::Dict{DNA,Score})
    return minimum((v, k) for (k, v) in d)[2]
end

function has_single_indels(consensus::DNASeq,
                           reference::RifrafSequence)
    moves = align_moves(consensus, reference)
    return TRACE_INSERT in moves || TRACE_DELETE in moves
end

function single_indel_proposals(consensus::DNASeq,
                                reference::RifrafSequence)
    moves = align_moves(consensus, reference,
                        skew_matches=true)
    results = Proposal[]
    cons_idx = 0
    ref_idx = 0
    for move in moves
        if move == TRACE_MATCH
            cons_idx += 1
            ref_idx += 1
        elseif move == TRACE_INSERT
            ref_idx += 1
            push!(results, Insertion(cons_idx, reference.seq[ref_idx]))
        elseif move == TRACE_DELETE
            cons_idx += 1
            push!(results, Deletion(cons_idx))
        elseif move == TRACE_CODON_INSERT
            ref_idx += 3
        elseif move == TRACE_CODON_DELETE
            cons_idx += 3
        end
    end
    return results
end

function initial_state(consensus::DNASeq,
                       sequences::Vector{RifrafSequence},
                       reference::DNASeq,
                       params::RifrafParams)
    # batch size calculations
    batch_size = params.batch_size > 1 ? params.batch_size : length(sequences)
    batch_size = min(batch_size, length(sequences))
    batch_fixed_size = min(params.batch_fixed_size,
                           length(sequences))
    base_batch_size = batch_size

    if length(consensus) == 0
        # choose highest-quality sequence as initial consensus
        idx = indmax([logsumexp10(s.match_scores) for s in sequences])
        consensus = sequences[idx].seq
    end

    maxlen = maximum(map(length, sequences))
    rlen = length(reference)

    # initialize sequence banded arrays
    As = BandedArray{Score}[]
    Bs = BandedArray{Score}[]
    Amoves = BandedArray{Trace}[]

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

    return RifrafState(consensus=consensus,
                       reference=refseq,
                       ref_scores=params.ref_scores,
                       batch_fixed_size=batch_fixed_size,
                       batch_size=batch_size,
                       base_batch_size=base_batch_size,
                       maxlen=maxlen,
                       sequences=sequences,
                       As=As, Bs=Bs, Amoves=Amoves,
                       A_ref=A_ref, B_ref=B_ref, Amoves_ref=Amoves_ref)
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

"""Forward alignment with bandwidth checking.

If the number of errors is greater than the `1 - pvalue` right tail of
a Poisson with `μ = seq.error_rate`, increase bandwidth.

"""
function smart_forward_moves!(consensus::DNASeq, seq::RifrafSequence,
                              A::BandedArray{Score},
                              B::BandedArray{Score},
                              Amoves::BandedArray{Trace},
                              pvalue::Float64)
    max_bandwidth = min(seq.bandwidth * 2^5, length(consensus), length(seq))
    if seq.bandwidth_fixed
        max_bandwidth = seq.bandwidth
    end
    n_errors = typemax(Int)
    old_n_errors = typemax(Int)
    while true
        forward_moves!(consensus, seq, A, Amoves)
        if seq.bandwidth_fixed || seq.bandwidth >= max_bandwidth
            break
        end
        old_n_errors = n_errors
        n_errors = count_errors(Amoves, consensus, seq.seq)
        threshold = cquantile(Poisson(seq.est_n_errors), pvalue)
        if n_errors > threshold && n_errors < old_n_errors
            seq.bandwidth = min(seq.bandwidth * 2, max_bandwidth)
            newbandwidth!(A, seq.bandwidth)
            newbandwidth!(B, seq.bandwidth)
            newbandwidth!(Amoves, seq.bandwidth)
        else
            break
        end
    end
    seq.bandwidth_fixed = true
end

"""Do forward and backward alignments, and compute scores.

Do smart bandwidth updating.

"""
function realign!(state::RifrafState, params::RifrafParams)
    seqs = state.batch_seqs
    shape = (state.maxlen + 1, state.maxlen + 1)
    # add new arrays, in case batch size increased
    for i in (length(state.As) + 1):length(seqs)
        push!(state.As, BandedArray(Score, shape, seqs[i].bandwidth))
        push!(state.Bs, BandedArray(Score, shape, seqs[i].bandwidth))
        push!(state.Amoves, BandedArray(Trace, shape, seqs[i].bandwidth))
    end

    if state.realign_As
        if params.verbose >= 2
            println(STDERR, "    realigning As")
        end
        for (s, A, B, Amoves) in zip(seqs, state.As, state.Bs, state.Amoves)
            smart_forward_moves!(state.consensus, s, A, B, Amoves,
                                 params.bandwidth_pvalue)
        end
        if use_ref(state.reference, state.stage, params.use_ref_for_qvs)
            smart_forward_moves!(state.consensus, state.reference,
                                 state.A_ref, state.B_ref, state.Amoves_ref,
                                 params.bandwidth_pvalue)
        end
    end
    if state.realign_Bs
        if params.verbose >= 2
            println(STDERR, "    realigning Bs")
        end
        for (s, B) in zip(seqs, state.Bs)
            backward!(state.consensus, s, B)
        end
        if use_ref(state.reference, state.stage, params.use_ref_for_qvs)
            backward!(state.consensus, state.reference, state.B_ref)
        end
    end
end

function realign_rescore!(state::RifrafState, params::RifrafParams)
    realign!(state, params)
    rescore!(state, params.use_ref_for_qvs)
end

"""convert per-proposal log differences to a per-base error rate"""
function normalize_log_differences(sub_scores,
                                   del_scores,
                                   ins_scores,
                                   state_score)
    # per-base insertion score is mean of neighboring insertions
    pos_scores = hcat(sub_scores, del_scores)
    pos_exp = exp10.(pos_scores)
    pos_probs = broadcast(/, pos_exp, sum(pos_exp, 2))
    ins_exp = exp10.(ins_scores)
    ins_probs = broadcast(/, ins_exp, exp10.(state_score) + sum(ins_exp, 2))
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
                             no_ins_error_prob[1:end - 1],
                             no_ins_error_prob[2:end])
    return reshape(result, length(result))
end

function base_distribution(base::DNA, ilp)
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
    probs = exp10.(probs)
    probs = 1.0 - maximum(probs ./ sum(probs, 2), 2)
    return vec(probs)
end

function check_params(scores, reference, params)
    if (scores.mismatch >= 0.0 || scores.mismatch == -Inf ||
        scores.insertion >= 0.0 || scores.insertion == -Inf ||
        scores.deletion >= 0.0 || scores.deletion == -Inf)
        error("scores must be between -Inf and 0.0")
    end
    if scores.codon_insertion > -Inf || scores.codon_deletion > -Inf
        error("error model cannot allow codon indels")
    end

    if length(reference) > 0
        if params.ref_error_mult <= 0.0
            error("ref_error_mult must be > 0.0")
        end
        if params.ref_indel_mult <= 0.0
            error("ref_indel_mult must be > 0.0")
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
        if (params.max_ref_indel_mults < 0)
            error("ref_indel_increases must be >= 0")
        end
    end

    if !any([params.do_init, params.do_frame, params.do_refine, params.do_score])
        error("no stages enabled")
    end
    if params.max_iters < 1
        error("invalid max iters: $(params.max_iters)")
    end
    if params.batch_fixed && params.batch_fixed_size <= 1
        error("batch_fixed_size must be > 1")
    end
    if params.batch_randomness < 0.0 || params.batch_randomness > 1.0
        error("batch_randomness must be between 0.0 and 1.0")
    end
    if params.batch_mult < 0.0 || params.batch_mult > 1.0
        error("batch_mult must be between 0.0 and 1.0")
    end
    if params.batch_threshold < 0.0 || params.batch_mult > 1.0
        error("batch_threshold must be between 0.0 and 1.0")
    end
end

function handle_candidates!(candidates::Vector{ScoredProposal},
                            state::RifrafState,
                            params::RifrafParams)
    old_consensus = state.consensus
    chosen_cands = choose_candidates(candidates, params.min_dist)
    if params.verbose >= 2
        println(STDERR, "    found $(length(candidates)) candidates; filtered to $(length(chosen_cands))")
    end
    if params.verbose >= 3
        println(STDERR, "    chosen: $chosen_cands")
    end
    state.consensus = apply_proposals(old_consensus,
                                      Proposal[c.proposal
                                               for c in chosen_cands])
    state.realign_As = true
    state.realign_Bs = false
    realign_rescore!(state, params)
    # detect if a single proposal is better
    # note: this may not always be correct, because score_proposal() is not exact
    if length(chosen_cands) > 1 &&
        (state.score < chosen_cands[1].score ||
         isapprox(state.score, chosen_cands[1].score))
        if params.verbose >= 2
            println(STDERR, "    rejecting multiple candidates in favor of best")
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
    if params.verbose >= 2
        println(STDERR, "    no candidates found in $(state.stage).")
    end
    if state.stage == STAGE_INIT
        if length(state.reference) == 0 || !(params.do_frame)
            state.converged = true
        else
            state.stage = STAGE_FRAME
            # estimate reference error rate
            edit_dist = edit_distance(state.consensus, state.reference.seq)
            ref_error_rate = edit_dist / max(length(state.reference),
                                             length(state.consensus))
            ref_error_rate *= params.ref_error_mult
            # needs to be < 0.5, otherwise matches aren't rewarded at all
            state.ref_error_rate = min(max(ref_error_rate, 1e-10), 0.5)
            ref_error_log_p = fill(log10(state.ref_error_rate), length(state.reference))
            state.reference = RifrafSequence(state.reference.seq, ref_error_log_p,
                                             params.bandwidth, state.ref_scores)

            # check if we've already converged
            if !has_single_indels(state.consensus, state.reference)
                state.converged = true
            end
        end
    elseif state.stage == STAGE_FRAME
        if !has_single_indels(state.consensus, state.reference)
            state.stage = STAGE_REFINE
        elseif (state.n_ref_indel_mults == params.max_ref_indel_mults)
            if params.verbose >= 2
                println(STDERR, "    NOTE: alignment had single indels but reached penalty limit.")
            end
            state.stage = STAGE_REFINE
        else
            state.penalties_increased = true
            if state.n_ref_indel_mults < params.max_ref_indel_mults
                state.n_ref_indel_mults += 1
            else
                error("Tried to illegally increase n_ref_indel_mults")
            end
            # TODO: this is not probabilistically correct
            mult = params.ref_indel_mult ^ state.n_ref_indel_mults
            state.ref_scores = Scores(state.ref_scores.mismatch,
                                      state.ref_scores.insertion * mult,
                                      state.ref_scores.deletion * mult,
                                      state.ref_scores.codon_insertion,
                                      state.ref_scores.codon_deletion)
            state.reference = RifrafSequence(state.reference, state.ref_scores)
            if params.verbose >= 2
                println(STDERR, "    NOTE: alignment to reference had single indels. increasing penalty.")
            end
        end
    elseif state.stage == STAGE_REFINE
        state.converged = true
    else
        error("  invalid stage: $(state.stage)")
    end
end

"""Samble a subset of `data`"""
function resample(data, n, wv)
    # TODO: optionally softmax the weight vector
    if n < length(data)
        return StatsBase.sample(data, wv, n, replace=false)
    else
        return data
    end
end

"""Reweight weight vector to increase or decrease the weights of the
top `n` values.

`randomness` can be between 0 and 1

0: top n get picked
0.5: weight according to estimated errors
1: completely random

"""
function reweight(wv, n, randomness)
    # TODO: there is probably a one liner that does all this
    if randomness < 0.0 || randomness > 1.0
        error("randomness must be between 0.0 and 1.0")
    end
    wv = wv ./ sum(wv)
    indices = reverse(sortperm(wv))[1:n]
    endpoint = wv
    weight = 0.0
    if randomness > 0.5
        weight = (randomness - 0.5) * 2.0
        endpoint = fill(1.0 / length(wv), length(wv))
    elseif randomness < 0.5
        weight = 1.0 - randomness * 2.0
        endpoint = zeros(length(wv))
        endpoint[indices] = 1.0 / n
    end
    result = weight * endpoint + (1.0 - weight) * wv
    return result
end

function resample!(state::RifrafState, params::RifrafParams;
                   verbose::Int=0)
    err_weights = [s.est_n_errors for s in state.sequences]
    if ((state.stage == STAGE_INIT || state.stage == STAGE_FRAME) &&
        params.batch_fixed)
        # always return the top `n`
        # TODO: when to switch to random batches?
        indices = sortperm(err_weights)[1:state.batch_fixed_size]
        state.batch_seqs = state.sequences[indices]
        if verbose >= 2
            println(STDERR, "    kept fixed batch")
        end
        return
    end
    wv = Weights(reweight(1.0 - err_weights ./ sum(err_weights),
                          state.batch_size, state.batch_randomness))
    state.batch_seqs = resample(state.sequences, state.batch_size, wv)
    did_sample = state.batch_size < length(state.sequences)
    if did_sample
        state.realign_As = true
    end
    if verbose >= 2
        if did_sample
            println(STDERR, "    sampled $(state.batch_size) new sequences")
        else
            println(STDERR, "    sampled all sequences")
        end
    end
end

"""Check that score did not decrease.

If it did dicrease too much, increase batch size, resample, and
realign again.

"""
function check_score!(state::RifrafState, params::RifrafParams,
                      old_score::Float64)
    if params.verbose >= 2
        println(STDERR, "    score: $(state.score)")
    end

    if (!state.penalties_increased &&
        state.batch_size == length(state.sequences) &&
        state.stage_iterations[Int(state.stage)] > 1)
        if state.score < old_score
            if params.verbose >= 2
                println(STDERR, "    WARNING: not using batches, but score decreased.")
            end
        elseif state.score == old_score
            if params.verbose >= 2
                println(STDERR, "    score did not change. ending stage.")
            end
            return false
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
        if params.verbose >= 2
            println(STDERR, "    NOTE: increased batch size to $(state.batch_size).")
        end
        resample!(state, params, verbose=params.verbose)
        state.realign_As = true
        state.realign_Bs = true
        realign_rescore!(state, params)
        if params.verbose >= 2
            println(STDERR, "    new score: $(state.score)")
        end
    end
    return true
end

function rifraf(dnaseqs::Vector{DNASeq},
                error_log_ps::Vector{Vector{LogProb}};
                consensus::DNASeq=DNASeq(),
                reference::DNASeq=DNASeq(),
                params::RifrafParams=RifrafParams())
    check_params(params.scores, reference, params)

    sequences = RifrafSequence[RifrafSequence(s, p, params.bandwidth, params.scores)
                               for (s, p) in zip(dnaseqs, error_log_ps)]

    state = initial_state(consensus, sequences, reference, params)

    # TODO: this is ugly
    enabled_stages = Set{Stage}()
    if params.do_init
        push!(enabled_stages, STAGE_INIT)
    end
    if params.do_frame
        push!(enabled_stages, STAGE_FRAME)
    end
    if params.do_refine
        push!(enabled_stages, STAGE_REFINE)
    end
    if params.do_score
        push!(enabled_stages, STAGE_SCORE)
    end

    # keep track of every iteration's consensus
    consensus_stages = [[] for _ in 1:(Int(typemax(Stage)) - 1)]

    state.realign_As = true
    state.realign_Bs = true
    old_score = -Inf

    for iter in 1:params.max_iters
        # skip to next valid stage
        while state.stage < STAGE_SCORE && !(state.stage in enabled_stages)
            state.stage = next_stage(state.stage)
        end
        if state.stage == STAGE_SCORE
            break
        end
        state.stage_iterations[Int(state.stage)] += 1
        push!(consensus_stages[Int(state.stage)], state.consensus)

        if params.verbose >= 1
            println(STDERR, "iteration $iter : $(state.stage) : $(state.score)")
        end
        if params.verbose >= 3
            println(STDERR, "  consensus: $(state.consensus)")
        elseif params.verbose >= 2
            println(STDERR, "  consensus length: $(length(state.consensus))")
        end

        # resample sequences, realign, and possibly adjust batch size
        if params.verbose >= 2
            println(STDERR, "  step: resample")
        end
        resample!(state, params; verbose=params.verbose)

        if params.verbose >= 2
            println(STDERR, "  step: realign and rescore")
        end
        realign_rescore!(state, params)

        if params.verbose >= 2
            println(STDERR, "  step: check score")
        end
        if check_score!(state, params, old_score)
            old_score = state.score

            # get candidate changes to consensus
            state.penalties_increased = false

            indel_seeds = if state.stage == STAGE_FRAME && params.seed_indels
                single_indel_proposals(state.consensus, state.reference)
            else
                Proposal[]
            end
            candidates = get_candidates(state, params; indel_seeds=indel_seeds)

            # handle candidates
            state.realign_As = true
            if length(candidates) > 0
                if params.verbose >= 2
                    println(STDERR, "  step: handle candidates")
                end
                handle_candidates!(candidates, state, params)
            else
                if params.verbose >= 2
                    println(STDERR, "  step: finish_stage")
                end
                finish_stage!(state, params)
            end
        else
            finish_stage!(state, params)
        end
        if state.converged
            break
        end

        # update batch randomness
        if ((!params.batch_fixed ||
            (state.stage == STAGE_REFINE &&
             state.stage_iterations[Int(STAGE_REFINE)] > 1)) &&
            state.batch_size < length(state.sequences))
            state.batch_randomness *= params.batch_mult
            if params.verbose >= 2
                println(STDERR, "  batch randomness decreased to $(state.batch_randomness)")
            end
        end
    end
    state.stage = STAGE_SCORE

    result = RifrafResult(consensus=state.consensus,
                          params=params,
                          state=state,
                          consensus_stages=consensus_stages)

    if params.do_score
        if params.verbose >= 2
            println(STDERR, "computing consensus quality scores")
        end

        # FIXME: recomputing for all sequences is costly, but using batch
        # is less accurate
        # possibly use top n sequences here
        state.realign_As = true
        state.realign_Bs = true
        realign_rescore!(state, params)
        result.error_probs = estimate_probs(state, params.use_ref_for_qvs)
        result.aln_error_probs = alignment_error_probs(length(state.consensus),
                                                     state.batch_seqs, state.Amoves)
    end
    if params.verbose >= 1
        println(STDERR, "done. converged: $(state.converged)")
    end
    return result
end

"""
    rifraf(dnaseqs, phreds; kwargs...)

Find a consensus sequence for a set of DNA sequences.

Returns an instance of `RifrafResult`.

# Arguments
- `dnaseqs::Vector{DNASeq}`: reads for which to find a consensus

- `phreds::Vector{Vector{Phred}}`: Phred scores for `dnaseqs`

- `consensus::DNASeq=DNASeq()`: initial consensus; if not given,
  defaults to the sequence in `dnaseqs` with the lowest mean error
  rate

- `reference::DNASeq=DNASeq()`: reference for frame correction

- `params::RifrafParams=RifrafParams()`

"""
function rifraf(dnaseqs::Vector{DNASeq},
                phreds::Vector{Vector{Phred}};
                kwargs...)
    # TODO: do not require DNASeq only!
    if any(minimum(p) < 0 for p in phreds)
        error("phred score cannot be negative")
    end
    error_log_ps = phred_to_log_p(phreds)
    myseqs = DNASeq[DNASeq(s) for s in dnaseqs]
    return rifraf(myseqs, error_log_ps; kwargs...)
end

"""

Adjust error probabilities so that the expected number of errors
matches the edit distance to the consensus sequence.

"""
function calibrate_phreds(s::DNASeq, phred::Vector{Phred}, consensus::DNASeq)
    # TODO: should use our own alignments here
    n_errors = edit_distance(consensus, s)
    errors = phred_to_p(phred)
    return errors * float(n_errors) / sum(errors)
end

"""rifraf-style fast frameshift correction"""
function correct_shifts(consensus::DNASeq,
                        reference::DNASeq;
                        log_p::LogProb=-1.0,
                        bandwidth::Int=-1,
                        scores::Scores=Scores(ErrorModel(10.0, 1e-5, 1e-5, 1.0, 1.0)))
    log_ps = fill(log_p, length(reference))
    if bandwidth < 0
        bandwidth = Int(ceil(min(length(consensus), length(reference)) * 0.1))
    end
    refseq = RifrafSequence(reference, log_ps, bandwidth, scores)
    a, b = align(consensus, refseq)
    proposals = single_indel_proposals(consensus, refseq)
    result = apply_proposals(consensus, proposals)
end
