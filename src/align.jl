@enum(DPMove, dp_none=0, dp_match=1, dp_ins=2, dp_del=3,
      dp_codon_ins=4, dp_codon_del=5)


# all offsets are relative to the consensus.
# so an insertion is a base NOT in the consensus.
const OFFSETS = ([1, 1],  # sub
                 [1, 0],  # insertion
                 [0, 1],  # deletion
                 [3, 0],  # codon insertion
                 [0, 3])  # codon deletion

function move_scores(t_base::DNANucleotide,
                     s_base::DNANucleotide,
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
    offset = OFFSETS[Int(move)]
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
                s_base::DNANucleotide, t_base::DNANucleotide,
                pseq::RifrafSequence,
                scores::Scores;
                newcols::Array{Float64, 2}=Array(Float64, (0, 0)),
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
    result = update_helper(newcols, A, i, j, acol, dp_match, match_score,
                           result...)
    result = update_helper(newcols, A, i, j, acol, dp_ins, ins_score, result...)
    result = update_helper(newcols, A, i, j, acol, dp_del, del_score, result...)
    codon_ins_score, codon_del_score = codon_move_scores(i-1, pseq.error_log_p,
                                                         scores)

    if scores.codon_insertion > -Inf && i > CODON_LENGTH
        result = update_helper(newcols, A, i, j, acol, dp_codon_ins,
                               codon_ins_score, result...)
    end
    if scores.codon_deletion > -Inf && j > CODON_LENGTH
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
function forward_moves(t::DNASequence, s::RifrafSequence,
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
            sbase = i > 1 ? s.seq[i-1] : DNA_Gap
            tbase = j > 1 ? t[j-1] : DNA_Gap
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
function forward(t::DNASequence, s::RifrafSequence,
                 scores::Scores)
    result = BandedArray(Float64, (length(s) + 1, length(t) + 1), s.bandwidth)
    nrows, ncols = size(result)
    for j = 1:ncols
        start, stop = row_range(result, j)
        for i = start:stop
            if i == 1 && j == 1
                continue
            end
            sbase = i > 1 ? s.seq[i-1] : DNA_Gap
            tbase = j > 1 ? t[j-1] : DNA_Gap
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
function backward(t::DNASequence, s::RifrafSequence,
                  scores::Scores)
    result = forward(reverse(t), reverse(s), scores)
    return flip(result)
end


function backtrace(moves::BandedArray{Int})
    taken_moves = DPMove[]
    i, j = moves.shape
    while i > 1 || j > 1
        m = moves[i, j]
        push!(taken_moves, DPMove(m))
        ii, jj = OFFSETS[m]
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
        ii, jj = OFFSETS[m]
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

function moves_to_aligned_seqs(moves::Vector{DPMove},
                               t::DNASequence, s::DNASequence)
    aligned_t = DNASequence()
    aligned_s = DNASequence()
    i, j = (0, 0)
    for move in moves
        (ii, jj) = OFFSETS[Int(move)]
        i += ii
        j += jj
        if move == dp_match
            push!(aligned_t, t[j])
            push!(aligned_s, s[i])
        elseif move == dp_ins
            push!(aligned_t, DNA_Gap)
            push!(aligned_s, s[i])
        elseif move == dp_del
            push!(aligned_t, t[j])
            push!(aligned_s, DNA_Gap)
        elseif move == dp_codon_ins
            append!(aligned_t, [DNA_Gap, DNA_Gap, DNA_Gap])
            append!(aligned_s, [s[i-2], s[i-1], s[i]])
        elseif move == dp_codon_del
            append!(aligned_t, [t[j-2], t[j-1], t[j]])
            append!(aligned_s, [DNA_Gap, DNA_Gap, DNA_Gap])
        end
    end
    return aligned_t, aligned_s
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
        (ii, jj) = OFFSETS[Int(move)]
        i += ii
        j += jj
    end
    if j > last_j
        result[j] = i
        last_j = j
    end
    return result
end

function align_moves(t::DNASequence, s::RifrafSequence,
                     scores::Scores;
                     trim::Bool=false,
                     skew_matches::Bool=false)
    A, Amoves = forward_moves(t, s, scores;
                              trim=trim, skew_matches=skew_matches)
    return backtrace(Amoves)
end


function align(t::DNASequence, s::RifrafSequence,
               scores::Scores;
               trim::Bool=false,
               skew_matches::Bool=false)
    moves = align_moves(t, s, scores, trim=trim,
                        skew_matches=skew_matches)
    return moves_to_aligned_seqs(moves, t, s.seq)
end

function align(t::DNASequence, s::DNASequence, phreds::Vector{Int8},
               scores::Scores,
               bandwidth::Int;
               trim::Bool=false,
               skew_matches::Bool=false)
    moves = align_moves(t, RifrafSequence(s, phreds, bandwidth), scores,
                        trim=trim, skew_matches=skew_matches)
    return moves_to_aligned_seqs(moves, t, s)
end

function moves_to_proposals(moves::Vector{DPMove},
                            consensus::DNASequence, seq::RifrafSequence)
    proposals = Proposal[]
    i, j = (0, 0)
    for move in moves
        (ii, jj) = OFFSETS[Int(move)]
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
function alignment_proposals(Amoves::Vector{BandedArray{Int}},
                             consensus::DNASequence,
                             sequences::Vector{RifrafSequence},
                             do_subs::Bool,
                             do_indels::Bool)
    result = Set{Proposal}()
    for (Amoves, seq) in zip(Amoves, sequences)
        moves = backtrace(Amoves)
        for proposal in moves_to_proposals(moves, consensus, seq)
            if (typeof(proposal) == Substitution && do_subs) || do_indels
                push!(result, proposal)
            end
        end
    end
    return collect(result)
end
