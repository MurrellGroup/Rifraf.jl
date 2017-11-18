import Base.zero

# trace type for pairwise alignment
const Trace = Int8

# trace values
const TRACE_NONE         = Trace(0)
const TRACE_MATCH        = Trace(1)
const TRACE_INSERT       = Trace(2)
const TRACE_DELETE       = Trace(3)
const TRACE_CODON_INSERT = Trace(4)
const TRACE_CODON_DELETE = Trace(5)

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

function update_helper(final_score::Score, final_move::Trace,
                       move_score::Score, move::Trace,
                       newcols::Array{Score,2},
                       A::BandedArray{Score},
                       i::Int, j::Int, acol::Int)
    prev_i, prev_j = offset_backward(move, i, j)
    rangecol = min(prev_j, size(A)[2])
    if inband(A, prev_i, rangecol)
        score = if acol < 1 || prev_j <= acol
            A[prev_i, prev_j] + move_score
        else
            newcols[prev_i, prev_j - acol] + move_score
        end
        if score > final_score
            return score, move
        end
    end
    return final_score, final_move
end

function update(A::BandedArray{Score},
                i::Int, j::Int,
                s_base::DNA, t_base::DNA,
                pseq::RifrafSequence;
                newcols::Array{Score,2}=Array{Score}((0, 0)),
                doreverse::Bool=false,
                acol::Int=-1, trim::Bool=false,
                skew_matches::Bool=false)
    final_score = Score(-Inf)
    final_move = TRACE_NONE

    nrows, ncols = size(A)
    seqlen = length(pseq)
    # TODO: this cannot make mismatches preferable to codon indels
    seq_i = doreverse ? min(seqlen, seqlen - (i - 1) + 1) : max(i - 1, 1)
    del_i = doreverse ? nrows - i + 1 : i
    match_score = (s_base == t_base) ? pseq.match_scores[seq_i] : pseq.mismatch_scores[seq_i]
    ins_score = pseq.ins_scores[seq_i]
    del_score = pseq.del_scores[del_i]

    if skew_matches && s_base != t_base
        match_score *= 0.99
    end
    # allow terminal insertions for free
    if trim && ((j == 1) || (j == ncols))
        ins_score = 0.0
    end
    final_score, final_move = update_helper(final_score, final_move,
                                            match_score, TRACE_MATCH,
                                            newcols, A, i, j, acol)
    final_score, final_move = update_helper(final_score, final_move,
                                            ins_score, TRACE_INSERT,
                                            newcols, A, i, j, acol)
    final_score, final_move = update_helper(final_score, final_move,
                                            del_score, TRACE_DELETE,
                                            newcols, A, i, j, acol)

    if do_codon_moves(pseq)
        if do_codon_ins(pseq) && i > CODON_LENGTH
            codon_i = i - CODON_LENGTH
            if doreverse
                codon_i = length(pseq.codon_ins_scores) - codon_i + 1
            end
            codon_ins_score = pseq.codon_ins_scores[codon_i]
            final_score, final_move = update_helper(final_score, final_move,
                                                    codon_ins_score, TRACE_CODON_INSERT,
                                                    newcols, A, i, j, acol)
        end
        if do_codon_del(pseq) && j > CODON_LENGTH
            codon_del_score = pseq.codon_del_scores[del_i]
            final_score, final_move = update_helper(final_score, final_move,
                                                    codon_del_score, TRACE_CODON_DELETE,
                                                    newcols, A, i, j, acol,)
        end
    end
    if final_score == Score(-Inf)
        error("new score is invalid")
    end
    if final_move == TRACE_NONE
        error("failed to find a move")
    end
    return final_score, final_move
end

function forward_moves!(t::DNASeq, s::RifrafSequence,
                        result::BandedArray{Score},
                        moves::BandedArray{Trace};
                        trim::Bool=false,
                        skew_matches::Bool=false)
    new_shape = (length(s) + 1, length(t) + 1)
    newbandwidth!(result, s.bandwidth)
    newbandwidth!(moves, s.bandwidth)
    resize!(result, new_shape)
    resize!(moves, new_shape)
    result[1, 1] = Score(0.0)
    moves[1, 1] = TRACE_NONE
    nrows, ncols = new_shape
    for j = 1:ncols
        start, stop = row_range(result, j)
        for i = start:stop
            if i == 1 && j == 1
                continue
            end
            sbase = i > 1 ? s.seq[i - 1] : DNA_Gap
            tbase = j > 1 ? t[j - 1] : DNA_Gap
            result[i, j], moves[i, j] = update(result, i, j, sbase, tbase, s;
                                               doreverse=false,
                                               trim=trim,
                                               skew_matches=skew_matches)
        end
    end
end

"""Does backtracing to find best alignment."""
function forward_moves(t::DNASeq, s::RifrafSequence;
                       trim::Bool=false,
                       skew_matches::Bool=false)
    result = BandedArray(Score, (length(s) + 1, length(t) + 1), s.bandwidth;
                         default=-Inf, initialize=false)
    moves = BandedArray(Trace, size(result), s.bandwidth, initialize=false)
    forward_moves!(t, s, result, moves,
                   trim=trim, skew_matches=skew_matches)
    return result, moves
end

function forward!(t::DNASeq, s::RifrafSequence,
                  result::BandedArray{Score};
                  doreverse=false,
                  trim::Bool=false,
                  skew_matches::Bool=false)
    new_shape = (length(s) + 1, length(t) + 1)
    newbandwidth!(result, s.bandwidth)
    resize!(result, new_shape)
    result[1, 1] = Score(0.0)
    nrows, ncols = new_shape
    for j = 1:ncols
        start, stop = row_range(result, j)
        for i = start:stop
            if i == 1 && j == 1
                continue
            end
            sbase = i > 1 ? s.seq[doreverse ? length(s) - (i - 1) + 1 : i - 1] : DNA_Gap
            tbase = j > 1 ? t[doreverse ? length(t) - (j - 1) + 1 : j - 1] : DNA_Gap
            result[i, j], _ = update(result, i, j, sbase, tbase, s;
                                     doreverse=doreverse,
                                     trim=trim,
                                     skew_matches=skew_matches)
        end
    end
end

"""
F[i, j] is the log probability of aligning s[1:i-1] to t[1:j-1].

"""
function forward(t::DNASeq, s::RifrafSequence;
                 doreverse::Bool=false,
                 trim::Bool=false,
                 skew_matches::Bool=false)
    result = BandedArray(Score, (length(s) + 1, length(t) + 1), s.bandwidth;
                         initialize=false, default=-Inf)
    forward!(t, s, result; doreverse=doreverse, trim=trim,
             skew_matches=skew_matches)
    return result
end

function backward!(t::DNASeq, s::RifrafSequence,
                   result::BandedArray{Score})
    # this actually seems faster than a dedicated backwards alignment
    # algorithm. possibly because of cache effects.
    forward!(t, s, result, doreverse=true)
    flip!(result)
end

"""
B[i, j] is the log probability of aligning s[i:end] to t[j:end].

"""
function backward(t::DNASeq, s::RifrafSequence)
    result = forward(t, s, doreverse=true)
    flip!(result)
    return result
end

function backtrace_indices(moves::BandedArray{Trace};
                           start::Tuple{Int,Int}=(0, 0))
    result = Tuple{Int,Int}[]
    i, j = start
    if i == 0 || j == 0
        i, j = size(moves)
    end
    while i > 1 || j > 1
        m = moves[i, j]
        i, j = offset_backward(m, i, j)
        push!(result, (i, j))
    end
    return reverse(result)
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

function count_errors(arr::BandedArray{Trace},
                      t::DNASeq, s::DNASeq)
    moves = backtrace(arr)
    a, b = moves_to_aligned_seqs(moves, t, s)
    return sum([1 for (x, y) in zip(a, b) if x != y])
end

function count_errors(t::DNASeq, s::RifrafSequence)
    _, Amoves = forward_moves(t, s, skew_matches=true)
    return count_errors(Amoves, t, s.seq)
end

"""should be close to Levenshtein"""
function edit_distance(t::DNASeq, s::DNASeq)
    log_ps = fill(-1.0, length(s))
    bandwidth = Int(ceil(min(length(t), length(s)) * 0.5))
    scores = Scores(ErrorModel(1.0, 1.0, 1.0))
    seq = RifrafSequence(s, log_ps, bandwidth, scores)
    _, Amoves = forward_moves(t, seq, skew_matches=true)
    return count_errors(Amoves, t, s)
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
        i, j = offset_backward(Amoves[i, j], i, j)
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
            append!(aligned_t, DNASequence([DNA_Gap, DNA_Gap, DNA_Gap]))
            append!(aligned_s, DNASequence([s[i - 2], s[i - 1], s[i]]))
        elseif move == TRACE_CODON_DELETE
            append!(aligned_t, DNASequence([t[j - 2], t[j - 1], t[j]]))
            append!(aligned_s, DNASequence([DNA_Gap, DNA_Gap, DNA_Gap]))
        end
    end
    return aligned_t, aligned_s
end

"""
    moves_to_indices(moves, tlen, slen)

Compute index vector mapping positions in `t` to positions in `s`.

If a position in `t` is an insertion relative to `s`, it maps to the
last mapped base in `s` (or to `0`, if it is the first base in `t`).

"""
function moves_to_indices(moves::Vector{Trace},
                          tlen::Int, slen::Int)
    result = Int[]
    i, j = (0, 0)
    last_j = 0
    for move in moves
        i, j = offset_forward(move, i, j)
        if j > last_j
            push!(result, i)
            last_j = j
        end
    end
    return result
end

function align_moves(t::DNASeq, s::RifrafSequence;
                     trim::Bool=false,
                     skew_matches::Bool=false)
    A, Amoves = forward_moves(t, s;
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

function align_score(t::DNASeq, s::DNASeq, phreds::Vector{Phred},
                     scores::Scores,
                     bandwidth::Int;
                     trim::Bool=false,
                     skew_matches::Bool=false)
    A, Amoves = forward_moves(t, RifrafSequence(s, phreds, bandwidth, scores);
                              trim=trim,
                              skew_matches=skew_matches)
    moves = backtrace(Amoves)
    return moves_to_aligned_seqs(moves, t, s), A[end, end]
end
