module Sample

using Bio.Seq
using Distributions

using Quiver2.QIO

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
                              beta=-1.0, codon=true)
    offset = codon ? 2 : 0
    # do deletions
    del_rate = error_rate * del_ratio / (offset + 1)
    del_d = Bernoulli(del_rate)
    seq = DNANucleotide[]
    skip = 0
    for j = 1:length(template)
        if skip > 0
            skip -= 1
        elseif j < (length(template) - offset) && rand(del_d) == 0
            if codon
                push!(seq, template[j])
                push!(seq, template[j+1])
                push!(seq, template[j+2])
                skip = 2
            else
                push!(seq, template[j])
            end
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
    if codon
        ins_ratio /= 3
    end
    # do insertions, drawing new quality scores
    for j = 1:length(seq)
        while rand(Bernoulli(probs[j] * ins_ratio)) == 1
            for k=1:(codon ? 3 : 1)
                push!(final_seq, rbase())
                push!(final_probs, rand(prob_d))
            end
        end
        push!(final_seq, seq[j])
        push!(final_probs, probs[j])
    end
    while rand(Bernoulli(probs[end] * ins_ratio)) == 1
        for k=1:(codon ? 3 : 1)
            push!(final_seq, rbase())
            push!(final_probs, rand(prob_d))
        end
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
    if len % 3 != 0
        error("Reference length must be a multiple of three")
    end
    reference = random_seq(len)

    sub_ratio = (sub_part / (sub_part + ins_part + del_part))
    ins_ratio = (ins_part / (sub_part + ins_part + del_part))
    del_ratio = (del_part / (sub_part + ins_part + del_part))

    error_d = Beta(error_rate_alpha, error_rate_beta)

    template_error_rate = rand(error_d) * max_error_rate
    template, template_log_p = sample_from_template(reference, template_error_rate,
                                                    sub_ratio, ins_ratio, del_ratio,
                                                    beta=quality_beta, codon=true)

    seqs = DNASequence[]
    log_ps = Vector{Float64}[]

    for i = 1:nseqs
        error_rate = rand(error_d) * max_error_rate
        seq, log_p = sample_from_template(template, error_rate,
                                          sub_ratio, ins_ratio, del_ratio,
                                          beta=quality_beta, codon=false)
        push!(seqs, seq)
        push!(log_ps, log_p)
    end
    return DNASequence(reference), DNASequence(template), template_log_p, seqs, log_ps
end

"""Write template into FASTA and sequences into FASTQ."""
function write_samples(filename, template, seqs, log_ps)
    write_template(string(filename, "-template.fasta"), template)
    write_fastq(string(filename, "-sequences.fastq"), seqs, log_ps)
end

"""Read template from FASTA and sequences from FASTQ."""
function read_samples(filename)
    template = read_template(string(filename, "-template.fasta"))
    seqs, log_ps = read_fastq(string(filename, "-sequences.fastq"))
    return template, seqs, log_ps
end

end
