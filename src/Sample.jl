module Sample

using Bio.Seq
using Distributions

using Quiver2.IO

export rbase, mutate_base, random_seq, sample_from_template, sample

function rbase()
    bases = [DNA_A, DNA_C, DNA_G, DNA_T]
    return bases[rand(1:4)]
end

function mutate_base(base::DNANucleotide)
    result = rbase()
    while result == base
        result = rbase()
    end
    return result
end

function random_seq(n)
    return DNASequence([rbase() for i in 1:n])
end

function sample_from_template(template, error_rate,
                              sub_ratio, ins_ratio, del_ratio;
                              beta=-1.0)
    # do deletions
    del_rate = error_rate * del_ratio
    del_d = Bernoulli(del_rate)
    seq = DNANucleotide[]
    for j = 1:length(template)
        if !convert(Bool, rand(del_d))
            push!(seq, template[j])
        end
    end

    # generate per-base error probs
    if beta < 0
        beta = 1 / error_rate
    end
    alpha = -beta * error_rate / (error_rate - 1)
    prob_d = Beta(alpha, beta)
    probs = rand(prob_d, length(seq))

    # do substitutions
    for j = 1:length(seq)
        if rand(Bernoulli(probs[j] * sub_ratio)) == 1
            seq[j] = mutate_base(seq[j])
        end
    end

    final_seq = DNANucleotide[]
    final_probs = Float64[]
    # do insertions, drawing new quality scores
    for j = 1:length(seq)
        if rand(Bernoulli(probs[j] * ins_ratio)) == 1
            push!(final_seq, rbase())
            push!(final_probs, rand(prob_d))
        end
        push!(final_seq, seq[j])
        push!(final_probs, probs[j])
    end
    # insertions at end
    if rand(Bernoulli(probs[end] * ins_ratio)) == 1
        push!(final_seq, rbase())
        push!(final_probs, rand(prob_d))
    end

    # round log_ps and cap value so PHRED scores fit in a UInt8, which
    # is used by BioJulia's FASTQRecord.
    minval = Float64(typemax(UInt8)) / (-10.0)
    return DNASequence(final_seq), max(round(log10(final_probs), 1), minval)
end

function sample(nseqs::Int, len::Int,
                max_error_rate::Float64,
                sub_part::Float64, ins_part::Float64, del_part::Float64;
                error_rate_alpha::Float64=5.0, error_rate_beta::Float64=1.0,
                quality_beta::Float64=-1.0)
    error_d = Beta(error_rate_alpha, error_rate_beta)
    template = random_seq(len)

    sub_ratio = (sub_part / (sub_part + ins_part + del_part))
    ins_ratio = (ins_part / (sub_part + ins_part + del_part))
    del_ratio = (del_part / (sub_part + ins_part + del_part))

    seqs = DNASequence[]
    log_ps = Vector{Float64}[]

    for i = 1:nseqs
        error_rate = rand(error_d) * max_error_rate
        seq, log_p = sample_from_template(template, error_rate,
                                          sub_ratio, ins_ratio, del_ratio,
                                          beta=quality_beta)
        push!(seqs, seq)
        push!(log_ps, log_p)
    end
    return DNASequence(template), seqs, log_ps
end

"""Write template into FASTA and sequences into FASTQ."""
function write_samples(filename, template, seqs, log_ps)
    write_template(string(filename, "-template.fasta"), template)
    write_sequences(string(filename, "-sequences.fastq"), seqs, log_ps)
end

"""Read template from FASTA and sequences from FASTQ."""
function read_samples(filename)
    template = read_template(string(filename, "-template.fasta"))
    seqs, log_ps = read_sequences(string(filename, "-sequences.fastq"))
    return template, seqs, log_ps
end

end
