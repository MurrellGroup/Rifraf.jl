using DataFrames
using Levenshtein
using BioSymbols
using BioSequences
using IterTools

using Rifraf


"""indices of reads with right length and error cutoff"""
function valid_read_indices(error_range, sequences, phreds, mean_errors; length_cutoff=LENGTH_CUTOFF)
    (min_error, max_error) = error_range
    median_length = median(map(length, sequences))
    return [i for i in 1:length(phreds) if ((min_error <= mean_errors[i] <= max_error) &&
                                            (abs(length(sequences[i]) - median_length) < length_cutoff))]
end

"""random sample of `nseqs` reads"""
function sample_reads(nseqs::Int,
                      cand_indices::Vector{Int},
                      sequences::Vector{DNASeq},
                      phreds::Vector{Vector{Rifraf.Phred}},
                      names::Vector{String},
                      mean_errors::Vector{Float64})
    indices = rand(cand_indices, nseqs)
    sampled_seqs = DNASeq[sequences[i] for i in indices]
    sampled_phreds = Vector{Rifraf.Phred}[phreds[i] for i in indices]
    sampled_names = String[names[i] for i in indices]
    if any([length(s) != length(p) for (s, p) in zip(sampled_seqs, sampled_phreds)])
        for (ss, sp) in zip(sampled_seqs, sampled_phreds)
            println("$(length(ss)), $(length(sp))")
            println("")
        end
    end
    return sampled_seqs, sampled_phreds, sampled_names, indices
end

"""alignment moves"""
function get_moves(seq, reference, align_errors)
    align_scores = Scores(align_errors)
    align_bandwidth = 100
    pseq = Rifraf.RifrafSequence(DNASeq(seq), fill(-1.0, length(seq)),
                                 align_bandwidth, align_scores)
    moves = Rifraf.align_moves(DNASeq(reference), pseq; trim=true)
    moves
end

"""get [start, stop] indices to trim read ends by aligning to reference"""
function trim_ends_indices(seq, reference)
    align_errors = ErrorModel(1e5, 1e-3, 1e-3, 0.0, 0.0)
    moves = get_moves(seq, reference, align_errors)
    x = 1
    while moves[x] == Rifraf.TRACE_INSERT
        x += 1
    end
    y = length(moves)
    while moves[y] == Rifraf.TRACE_INSERT
        y -= 1
    end
    n_end = length(moves) - y
    y = length(seq) - n_end
    return x, y
end

"""align to reference and trim, to remove primers"""
function trim_reads(seqs, phreds, reference)
    refstring = String(reference)
    final_seqs = DNASeq[]
    final_phreds = Vector{Rifraf.Phred}[]
    for (seq, phred) in zip(seqs, phreds)
        if length(seq) != length(phred)
            error("mismatched lengths")
        end
        seqstring = String(seq)

        # take reverse complement, if necessary
        d1 = Levenshtein.levenshtein(seqstring, refstring)
        d2 = Levenshtein.levenshtein(String(reverse_complement(seq)),
                                     refstring)
        if d1 > d2
            seq = reverse_complement(seq)
            phred = reverse(phred)
        end

        # trim
        x, y = trim_ends_indices(seq, reference)
        push!(final_seqs, seq[x:y])
        push!(final_phreds, phred[x:y])
    end
    return final_seqs, final_phreds
end

function shorten_seq_index(seq, reference, target_len; codon_moves=false)
    align_errors = if codon_moves
        ErrorModel(10, 1e-5, 1e-5, 1.0, 1.0)
    else
        ErrorModel(10., 1., 1., 0.0, 0.0)
    end
    moves = get_moves(seq, reference, align_errors)
    indices = Rifraf.moves_to_indices(moves, length(reference), length(seq))
    stop = indices[target_len]
    if stop % 3 > 0
        stop += (3 - (stop % 3))
        stop = min(length(reference), stop)
    end
    return stop
end

function insert_base_helper(template, seq, phreds, index, base, align_errors)
    moves = get_moves(seq, template, align_errors)
    indices = Rifraf.moves_to_indices(moves, length(template), length(seq))
    read_index = indices[index]

    error_rate = mean(Rifraf.phred_to_p(phreds))
    final_seq, final_phreds = if rand(Bernoulli(error_rate)) == 1
        # small chance of sequencing error causing a deletion in the read
        # chance of missing base depends on sequence error rate
        seq, phreds
    else
        final_seq = DNASeq(seq[1:(read_index-1)], DNASeq([base]), seq[read_index:end])
        # uniform between neighbor phreds.
        phred1 = phreds[max(read_index-1, 1)]
        phred2 = phreds[min(read_index, length(template))]
        new_phred = rand(min(phred1,phred2):max(phred1,phred2))
        final_phreds = vcat(phreds[1:(read_index-1)], [new_phred], phreds[read_index:end])
        final_seq, final_phreds
    end
    return final_seq, final_phreds
end

function delete_base_helper(template, seq, phreds, index, align_errors)
    moves = get_moves(seq, template, align_errors)
    indices = Rifraf.moves_to_indices(moves, length(template), length(seq))
    read_index = indices[index]
    final_seq = copy(seq)
    final_phreds = copy(phreds)
    deleteat!(final_seq, read_index)
    deleteat!(final_phreds, read_index)
    return final_seq, final_phreds
end

"""

Get a random position in `sequence`. If `in_hompolymer`, it must be in
a homopolymer region at least `length` bp long. Otherwise it must be
in a non-homopolymer region.

"""
function sample_position(sequence, in_homopolymer; hplength=4)
    # run-length encoding
    groups = IterTools.groupby(x -> x, sequence)
    lengths = map(length, groups)

    # the index of each group in the overall sequence
    indices = [i + 1 for i in cumsum(lengths)[1:(end-1)]]
    insert!(indices, 1, 1)

    pairs = zip(lengths, indices)

    pairs = if in_homopolymer
        # filter out illegal positions
        [p for p in pairs if p[1] >= hplength]
    else
        [p for p in pairs if p[1] < hplength]
    end
    len, start = rand(pairs)
    return start + rand(1:len) - 1
end


function insert_base(template, seqs, phreds, in_homopolymer)
    # sample an existing base
    base_index = sample_position(template, in_homopolymer)
    # equal probability of insertion on either side of that base
    index = rand([base_index, base_index + 1])
    base = if in_homopolymer
        # match homopolymer bases
        template[base_index]
    else
        # do not match neighbors
        rand([DNA_A, DNA_C, DNA_G, DNA_T])
    end

    final_template = DNASeq(DNASeq(template[1:(index-1)]), DNASeq([base]), DNASeq(template[index:end]))
    align_errors = ErrorModel(10., 1., 1., 0.0, 0.0)
    final_seqs = DNASeq[]
    final_phreds = Vector{Rifraf.Phred}[]
    for (seq, phred) in zip(seqs, phreds)
        seq, phred = insert_base_helper(template, seq, phred, index, base, align_errors)
        push!(final_seqs, seq)
        push!(final_phreds, phred)
    end
    return (final_template, final_seqs, final_phreds, index)
end


function delete_base(template, seqs, phreds, in_homopolymer)
    index = sample_position(template, in_homopolymer)
    final_template = DNASeq(DNASeq(template[1:(index-1)]), DNASeq(template[(index+1):end]))
    align_errors = ErrorModel(10., 1., 1., 0.0, 0.0)

    final_seqs = DNASeq[]
    final_phreds = Vector{Rifraf.Phred}[]
    for (seq, phred) in zip(seqs, phreds)
        seq, phred = delete_base_helper(template, seq, phred, index, align_errors)
        push!(final_seqs, seq)
        push!(final_phreds, phred)
    end
    return (final_template, final_seqs, final_phreds, index)
end
