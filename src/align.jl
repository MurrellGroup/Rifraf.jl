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

# generate 8 different versions of this function, for efficiency:
# - forward and backward
# - with and without codon moves
# - with and without acols
for isforward in [true, false]
    for use_codon in [true, false]
        for use_newcols in [true, false]
            if !isforward && use_newcols
                # we never compute new columns in the backwards case,
                # so no need to generate this function
                continue
            end
            op = isforward ? :(-) : :(+)
            direction = isforward ? "forward" : "backward"
            s1 = use_codon ? "_codon" : ""
            s2 = use_newcols ? "_newcols" : ""
            validrow = isforward ? :(i > 1) : :(i < size(A)[1])
            validcol = isforward ? :(j > 1) : :(j < size(A)[2])
            funcname = parse("update_$direction$s1$s2")
            seq_i = isforward ? :(i-1) : :(i)
            codon_ins_i = isforward ? :(i - CODON_LENGTH) : :(i)
            codon_ins_check = isforward ? :(i > CODON_LENGTH) : :(i <= nrows - CODON_LENGTH)
            codon_del_check = isforward ? :(j > CODON_LENGTH) : :(j <= ncols - CODON_LENGTH)

            # `use_newcols` is only used in the forward case, so no
            # need to do any other indexing in that case

            # TODO: write a view for a banded array column, so all this can be removed
            match_prev = use_newcols ? :(j - acol > 1 ? newcols[i-1, j-acol-1] : A[i-1, j-1]) : parse("A[i $op 1, j $op 1]")
            ins_prev = use_newcols ? :(newcols[i-1, j-acol]) : parse("A[i $op 1, j]")
            del_prev = use_newcols ? :(j - acol > 1 ? newcols[i, j-acol-1] : A[i, j-1]) : parse("A[i, j $op 1]")
            codon_ins_prev = use_newcols ? :(newcols[i-CODON_LENGTH, j - acol]) : parse("A[i $op CODON_LENGTH, j]")
            codon_del_prev = use_newcols ? :(j - acol - CODON_LENGTH > 0 ? newcols[i, j - acol - CODON_LENGTH] : A[i, j-CODON_LENGTH]) : parse("A[i, j $op CODON_LENGTH]")

            codon_block = use_codon ? quote
                if $codon_ins_check
                    cur_score = $codon_ins_prev + pseq.codon_ins_scores[$codon_ins_i]
                    if cur_score > final_score
                        final_score = cur_score
                        final_move = TRACE_CODON_INSERT
                    end
                end
                if $codon_del_check
                    cur_score = $codon_del_prev + pseq.codon_del_scores[i]
                    if cur_score > final_score
                        final_score = cur_score
                        final_move = TRACE_CODON_DELETE
                    end
                end
            end : :()

            eval(quote
                 function $funcname(A::BandedArray{Score},
                                    i::Int, j::Int,
                                    s_base::DNANucleotide, t_base::DNANucleotide,
                                    pseq::RifrafSequence;
                                    newcols::Array{Score, 2}=Array(Score, 0, 0),
                                    acol::Int=0,
                                    match_mult::Float64=0.0)
                 @myassert($validrow, "update called on gap row")
                 @myassert($validcol, "update called on gap column")

                 nrows, ncols = size(A)
                 match_score = s_base == t_base ? pseq.match_scores[$seq_i] : pseq.mismatch_scores[$seq_i]
                 # TODO: version without match mult
                 if match_mult > 0.0
                     # TODO: implement this
                 end

                 final_score = Score(-Inf)
                 final_move = TRACE_NONE

                 cur_score = $match_prev + match_score
                 if cur_score > final_score
                   final_score = cur_score
                   final_move = TRACE_MATCH
                 end

                 cur_score = $ins_prev + pseq.ins_scores[$seq_i]
                 if cur_score > final_score
                   final_score = cur_score
                   final_move = TRACE_INSERT
                 end

                 cur_score = $del_prev + pseq.del_scores[i]
                 if cur_score > final_score
                   final_score = cur_score
                   final_move = TRACE_DELETE
                 end

                 $codon_block

                 @myassert(final_score != Score(-Inf), "new score is invalid")
                 @myassert(final_move != TRACE_NONE, "failed to find a move")
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
    docodondel = do_codon_del(s)
    for j in 2 : min((result.h_offset + result.bandwidth + 1), length(t) + 1)
        result[1, j] = result[1, j-1] + s.del_scores[1]
        moves[1, j] = TRACE_DELETE
        if docodondel && j > CODON_LENGTH
            cand_score = result[1, j-CODON_LENGTH] + s.codon_del_scores[1]
            if result[1, j] < cand_score
                result[1, j] = cand_score
                moves[1, j] = TRACE_CODON_DELETE
            end
        end
    end
    # fill first column
    docodonins = do_codon_ins(s)
    for i in 2 : min(result.v_offset + result.bandwidth + 1, length(s) + 1)
        result[i, 1] = result[i-1, 1] + s.ins_scores[i-1]
        moves[i, 1] = TRACE_INSERT
        if docodonins && i > CODON_LENGTH
            cand_score = result[i-CODON_LENGTH, 1] + s.codon_ins_scores[i-CODON_LENGTH]
            if result[i, 1] < cand_score
                result[i, 1] = cand_score
                moves[i, 1] = TRACE_CODON_INSERT
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
        for i in max(start, 2):stop
            sbase = s.seq[i-1]
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
    result[end, end] = Score(0.0)

    # fill last row
    # TODO: handle trim
    docodondel = do_codon_ins(s)
    for j in length(t) : -1 :max(1, length(t) - result.h_offset - result.bandwidth + 1)
        result[end, j] = result[end, j+1] + s.del_scores[end]
        if docodondel && j <= new_shape[2]-CODON_LENGTH
            cand_score = result[end, j+CODON_LENGTH] + s.codon_del_scores[end]
            if result[end, j] < cand_score
                result[end, j] = cand_score
            end
        end
    end
    # fill last column
    docodonins = do_codon_ins(s)
    for i in length(s) : -1 : max(1, length(s) - result.v_offset - result.bandwidth + 1)
        result[i, end] = result[i+1, end] + s.ins_scores[i]
        if docodonins && i <= new_shape[1]-CODON_LENGTH
            cand_score = result[i+CODON_LENGTH, end] + s.codon_ins_scores[i]
            if result[i, end] < cand_score
                result[i, end] = cand_score
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

    second_last_row = new_shape[1] - 1
    second_last_col = new_shape[2] - 1
    for j = second_last_col : -1 : 1
        start, stop = row_range(result, j)
        for i in min(stop, second_last_row) : -1 : start
            sbase = s.seq[i]
            tbase = t[j]
            result[i, j], _ = update_function(result, i, j, sbase, tbase, s;
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
    backward!(t, s, result,
              padding=padding, trim=trim, skew_matches=skew_matches)
    return result
end


const OFFSETS = ([1, 1],  # sub
                 [1, 0],  # insertion
                 [0, 1],  # deletion
                 [3, 0],  # codon insertion
                 [0, 3])  # codon deletion


function offset_forward(move::Trace, i::Int, j::Int)
    (a, b) = OFFSETS[move]
    return (i + a, j + b)
end


function offset_backward(move::Trace, i::Int, j::Int)
    (a, b) = OFFSETS[move]
    return (i - a, j - b)
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
