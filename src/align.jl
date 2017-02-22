import Base.zero

# trace type for pairwise alignment
typealias Trace Int8

# trace values
const TRACE_NONE         = Trace(0)
const TRACE_MATCH        = Trace(1)
const TRACE_INSERT       = Trace(2)
const TRACE_DELETE       = Trace(3)
const TRACE_CODON_INSERT = Trace(4)
const TRACE_CODON_DELETE = Trace(5)


# all offsets are relative to the consensus.
# so an insertion is a base NOT in the consensus.

# matrix = use_acol ? :(newcols) : :(A)
# j_to_use = use_acol ? :(prev_j - acol) : :(prev_j)

# TODO: version without newcol
# TODO: get rid of inband() call
macro update_score(move, move_score, prev_cell)
    return quote
        score = $prev_cell + $move_score
        if score > final_score
            final_score = score
            final_move = $move
        end
    end
end


# generate 8 different versions of this function, for efficiency:
# - forward and backward
# - with and without codon moves
# - with and without acols
for direction in [:(forward), :(backward)]
    for use_codon in [true, false]
        for use_newcols in [true, false]
            op = direction == :forward ? :(-) : :(+)
            s1 = use_codon ? "_codon" : ""
            s2 = use_newcols ? "_newcols" : ""
            funcname = parse("update_$direction$s1$s2")
            # only change op for when use_newcols is false
            match_prev = use_newcols ? :(newcols[i-1, j - acol - 1]) : parse("A[i $op 1, j $op 1]")
            ins_prev = use_newcols ? :(newcols[i-1, j - acol]) : parse("A[i $op 1, j]")
            del_prev = use_newcols ? :(newcols[i, j - acol - 1]) : parse("A[i, j $op 1]")
            codon_ins_prev = use_newcols ? :(newcols[i-CODON_LENGTH, j - acol]) : parse("A[i $op CODON_LENGTH, j]")
            codon_del_prev = use_newcols ? :(newcols[i, j - acol - CODON_LENGTH]) : parse("A[i, j $op CODON_LENGTH]")

            eval(quote
                function $funcname(A::BandedArray{Score},
                                   i::Int, j::Int,
                                   s_base::DNANucleotide, t_base::DNANucleotide,
                                   pseq::RifrafSequence;
                                   newcols::Array{Score, 2}=Array(Score, 0, 0),
                                   acol::Int=0,
                                   match_mult::Float64=0.0)
                    final_score = Score(-Inf)
                    final_move = TRACE_NONE

                    match_score = s_base == t_base ? pseq.match_scores[i] : pseq.mismatch_scores[i]
                    # TODO: version without match mult
                    if match_mult > 0.0
                        # TODO: implement this
                    end

                    @update_score(TRACE_MATCH, match_score, $match_prev)
                    @update_score(TRACE_INSERT, pseq.ins_scores[i], $ins_prev)
                    @update_score(TRACE_DELETE, pseq.del_scores[i], $del_prev)

                    @includeif($use_codon, begin
                               if i > CODON_LENGTH
                               @update_score(TRACE_CODON_INSERT, pseq.ins_scores[i], $codon_ins_prev)
                               end
                               if j > CODON_LENGTH
                               @update_score(TRACE_CODON_DELETE, pseq.del_scores[i], $codon_del_prev)
                               end
                               end)

                    @myassert(final_score == Score(-Inf), "new score is invalid")
                    @myassert(final_move == TRACE_NONE, "failed to find a move")

                    return final_score, final_move
                end
            end)
        end
    end
end


function forward_moves!(t::DNASeq, s::RifrafSequence,
                        result::BandedArray{Score},
                        moves::BandedArray{Trace};
                        trim::Bool=false,
                        skew_matches::Bool=false)
    new_shape = (length(s) + 1, length(t) + 1)
    resize!(result, new_shape)
    resize!(moves, new_shape)
    result[1, 1] = Score(0.0)
    moves[1, 1] = TRACE_NONE

    # fill first row
    # TODO: handle trim
    do_codon_ins = length(s.codon_ins_scores) > 0
    for j in 2:(length(t) + 1)
        result[1, j] = result[1, j-1] + s.ins_scores[j]
        moves[1, j] = TRACE_INSERT
        if do_codon && j > CODON_LENGTH
            cand_score = result[1, j-CODON_LENGTH] + s.codon_ins_scores[j]
            if result[1, j] < cand_score
                result[1, j] = cand_score
                moves[1, j] = TRACE_CODON_INSERT
            end
        end
    end
    # fill first column
    do_codon_del = length(s.codon_del_scores) > 0
    for i in 2:length(s) + 1
        result[i, 1] = result[i-1, 1] + s.del_scores[j]
        moves[i, 1] = TRACE_DELETE
        if do_codon && i > CODON_LENGTH
            cand_score = result[i-CODON_LENGTH, 1] + s.codon_del_scores[j]
            if result[i, 1] < cand_score
                result[i, 1] = cand_score
                moves[i, 1] = TRACE_CODON_DELETE
            end
        end
    end

    # TODO: handle this
    # TODO: this cannot make mismatches preferable to codon indels
    match_mult = skew_matches ? 0.1 : 0.0

    update_function = update_forward
    if length(s.codon_ins_scores) > 0 || length(s.codon_del_scores) > 0
        update_function = update_forward_codon
    end

    for j = 2:new_shape[2]
        start, stop = row_range(result, j)
        for i = start:stop
            sbase = s.seq[i-i]
            tbase = t[j-1]
            result[i, j], moves[i, j] = update_function(result, i, j, sbase, tbase, s;
                                                        match_mult=match_mult)
        end
    end
end


"""Does backtracing to find best alignment."""
function forward_moves(t::DNASeq, s::RifrafSequence;
                       padding::Int=0,
                       trim::Bool=false,
                       skew_matches::Bool=false)
    result = BandedArray(Score, (length(s) + 1, length(t) + 1), s.bandwidth;
                         padding=padding,
                         initialize=false,
                         default=-Inf)
    moves = BandedArray(Trace, size(result), s.bandwidth,
                        padding=padding,
                        initialize=false)
    forward_moves!(t, s, result, moves,
                   trim=trim, skew_matches=skew_matches)
    return result, moves
end

# TODO: code duplication with forward
# use the eval trick here too
function backward!(t::DNASeq, s::RifrafSequence,
                   result::BandedArray{Score};
                   padding::Int=0,
                   trim::Bool=false,
                   skew_matches::Bool=false)
    new_shape = (length(s) + 1, length(t) + 1)
    resize!(result, new_shape)
    resize!(moves, new_shape)
    result[end, end] = Score(0.0)
    moves[end, end] = TRACE_NONE

    # fill last row
    # TODO: handle trim
    do_codon_ins = length(s.codon_ins_scores) > 0
    for j in length(t):-1:1
        result[end, j] = result[end, j+1] + s.ins_scores[j]
        moves[end, j] = TRACE_INSERT
        if do_codon && j > CODON_LENGTH
            cand_score = result[end, j+CODON_LENGTH] + s.codon_ins_scores[j]
            if result[end, j] < cand_score
                result[end, j] = cand_score
                moves[end, j] = TRACE_CODON_INSERT
            end
        end
    end
    # fill last column
    do_codon_del = length(s.codon_del_scores) > 0
    for i in length(s):-1:1
        result[i, end] = result[i+1, end] + s.del_scores[j]
        moves[i, end] = TRACE_DELETE
        if do_codon && i > CODON_LENGTH
            cand_score = result[i+CODON_LENGTH, end] + s.codon_del_scores[j]
            if result[i, end] < cand_score
                result[i, end] = cand_score
                moves[i, end] = TRACE_CODON_DELETE
            end
        end
    end

    # TODO: handle this
    # TODO: this cannot make mismatches preferable to codon indels
    match_mult = skew_matches ? 0.1 : 0.0

    update_function = update_backward
    if length(s.codon_ins_scores) > 0 || length(s.codon_del_scores) > 0
        update_function = update_backward_codon
    end

    for j = (new_shape[2]-1):1
        start, stop = row_range(result, j)
        for i = stop:-1:start
            sbase = s.seq[i+1]
            tbase = t[j+1]
            result[i, j], moves[i, j] = update_function(result, i, j, sbase, tbase, s;
                                                        match_mult=match_mult)
        end
    end
end


"""
B[i, j] is the log probability of aligning s[i:end] to t[j:end].

"""
function backward(t::DNASeq, s::RifrafSequence;
                  padding::Int=0,
                  trim::Bool=false,
                  skew_matches::Bool=false)
    result = BandedArray(Score, (length(s) + 1, length(t) + 1), s.bandwidth;
                         padding=padding,
                         initialize=false,
                         default=-Inf)
    result = backward!(t, s, result,
                     padding=padding, trim=trim, skew_matches=skew_matches)
    return result
end


function backtrace(moves::BandedArray{Trace})
    taken_moves = Trace[]
    i, j = size(moves)
    while i > 1 || j > 1
        m = moves[i, j]
        push!(taken_moves, m)
        i, j = offset_backward(m, i, j)
    end
    return reverse(taken_moves)
end


function band_tolerance(Amoves::BandedArray{Trace})
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


function moves_to_aligned_seqs(moves::Vector{Trace},
                               t::DNASeq, s::DNASeq)
    aligned_t = DNASequence()
    aligned_s = DNASequence()
    i, j = (0, 0)
    for move in moves
        i, j = offset_forward(move, i, j)
        if move == TRACE_MATCH
            push!(aligned_t, t[j])
            push!(aligned_s, s[i])
        elseif move == TRACE_INSERT
            push!(aligned_t, DNA_Gap)
            push!(aligned_s, s[i])
        elseif move == TRACE_DELETE
            push!(aligned_t, t[j])
            push!(aligned_s, DNA_Gap)
        elseif move == TRACE_CODON_INSERT
            append!(aligned_t, [DNA_Gap, DNA_Gap, DNA_Gap])
            append!(aligned_s, [s[i-2], s[i-1], s[i]])
        elseif move == TRACE_CODON_DELETE
            append!(aligned_t, [t[j-2], t[j-1], t[j]])
            append!(aligned_s, [DNA_Gap, DNA_Gap, DNA_Gap])
        end
    end
    return aligned_t, aligned_s
end


""" Compute index vector mapping from position in `t` to position in
`s`.

"""
function moves_to_indices(moves::Vector{Trace},
                          tlen::Int, slen::Int)
    result = zeros(Int, tlen + 1)
    i, j = (1, 1)
    last_j = 0
    for move in moves
        if j > last_j
            result[j] = i
            last_j = j
        end
        i, j = offset_forward(move, i, j)
    end
    if j > last_j
        result[j] = i
        last_j = j
    end
    return result
end


function align_moves(t::DNASeq, s::RifrafSequence;
                     padding::Int=0,
                     trim::Bool=false,
                     skew_matches::Bool=false)
    A, Amoves = forward_moves(t, s;
                              padding=padding,
                              trim=trim,
                              skew_matches=skew_matches)
    return backtrace(Amoves)
end


function align(t::DNASeq, s::RifrafSequence;
               trim::Bool=false,
               skew_matches::Bool=false)
    moves = align_moves(t, s;
                        trim=trim,
                        skew_matches=skew_matches)
    return moves_to_aligned_seqs(moves, t, s.seq)
end


function align(t::DNASeq, s::DNASeq, phreds::Vector{Phred},
               scores::Scores,
               bandwidth::Int;
               trim::Bool=false,
               skew_matches::Bool=false)
    moves = align_moves(t, RifrafSequence(s, phreds, bandwidth, scores),
                        trim=trim, skew_matches=skew_matches)
    return moves_to_aligned_seqs(moves, t, s)
end


function moves_to_proposals(moves::Vector{Trace},
                            consensus::DNASeq, seq::RifrafSequence)
    proposals = Proposal[]
    i, j = (0, 0)
    for move in moves
        i, j = offset_forward(move, i, j)

        score = seq.match_log_p[max(i, 1)]
        next_score = seq.match_log_p[min(i + 1, length(seq))]
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
