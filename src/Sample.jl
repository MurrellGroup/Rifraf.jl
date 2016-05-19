module Sample

using Bio.Seq
using Distributions

using Quiver2.QIO

export rbase, mutate_base, random_seq, sample_from_template, sample, BetaAlt

function rbase()
    bases = [DNA_A, DNA_C, DNA_G, DNA_T]
    return bases[Base.rand(1:4)]
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

function to_phred(x)
    return -10.0 * log10(x)
end

function to_prob(x)
    return exp10(-x / 10.0)
end

function normalize(x, lower, upper)
    return (x - lower) / (upper - lower)
end

function restore(y, lower, upper)
    return y * (upper - lower) + lower
end

"""A beta distribution parametrized by mean and standar deviation."""
type BetaAlt <: Sampleable{Univariate,Continuous}
    mean::Float64
    std::Float64

    minval::Float64
    maxval::Float64

    alpha::Float64
    beta::Float64
    d::Sampleable{Univariate,Continuous}

    function BetaAlt(mean::Float64, std::Float64; minval::Float64=0.0, maxval::Float64=1.0)
        if mean <= minval || mean >= maxval
            error("mean must be in ($minval, $maxval)\n")
        end
        mean_scaled = normalize(mean, minval, maxval)
        scale = (maxval - minval)
        std_scaled = std / scale
        if std_scaled <= 0
            error("std must be > 0")
        end
        variance = std_scaled ^ 2
        max_variance = mean_scaled * (1 - mean_scaled)
        if variance >= max_variance
            variance = max_variance * 0.9999
        end
        alpha = (((1 - mean_scaled) / variance) - 1 / mean_scaled) * (mean_scaled ^ 2)
        beta = alpha * (1 / mean_scaled - 1)
        d = Beta(alpha, beta)
        new(mean, std, minval, maxval, alpha, beta, d)
    end
end

function rand(s::BetaAlt)
    return restore(Distributions.rand(s.d), s.minval, s.maxval)
end

"""
Draw a vector of samples.

This should be done automatically by Distributions.jl, but it results
in a stack overflow.

"""
function myrand(s::BetaAlt, n::Int)
    return [rand(s) for i in 1:n]
end

function ratios(sub_part, ins_part, del_part)
    denom = sub_part + ins_part + del_part
    sub_ratio = sub_part / denom
    ins_ratio = ins_part / denom
    del_ratio = del_part / denom
    return sub_ratio, ins_ratio, del_ratio
end


function sample_from_template(template, error_rate,
                              sub_part, ins_part, del_part,
                              error_std;
                              codon=true)
    sub_ratio, ins_ratio, del_ratio = ratios(sub_part, ins_part, del_part)
    # do deletions
    del_rate = error_rate * del_ratio
    if codon
        del_rate /= 3
    end
    del_d = Bernoulli(del_rate)
    seq = DNANucleotide[]
    skip = 0
    for j = 1:length(template)
        if skip > 0
            skip -= 1
            continue
        end
        if Distributions.rand(del_d) == 0 || (codon && ((length(template) - j) < 2))
            push!(seq, template[j])
        else
            # do deletion
            if codon
                skip = 2
            end
        end
    end

    # generate per-base error probs
    min_prob = to_prob(Float64(typemax(UInt8)))
    error_d = BetaAlt(error_rate, error_std, minval=min_prob)
    probs = myrand(error_d, length(seq))

    # do substitutions
    if sub_ratio > 0
        sub_d = Bernoulli(error_rate * sub_ratio)
        for j = 1:length(seq)
            s_rate = probs[j] * sub_ratio
            if s_rate > 0 && Distributions.rand(Bernoulli(s_rate)) == 1
                seq[j] = mutate_base(seq[j])
            end
        end
    end

    final_seq = DNANucleotide[]
    final_probs = Float64[]
    if codon
        ins_ratio /= 3
    end
    # do insertions, drawing new phred scores
    if ins_ratio > 0
        for j = 1:length(seq)
            push!(final_seq, seq[j])
            push!(final_probs, probs[j])
            ins_rate = probs[j] * ins_ratio
            while Distributions.rand(Bernoulli(ins_rate)) == 1
                for k=1:(codon ? 3 : 1)
                    push!(final_seq, rbase())
                    push!(final_probs, rand(error_d))
                end
            end
        end
    else
        final_seq = seq
        final_probs = probs
    end
    # round phreds and cap value so they fit in a UInt8, which
    # is used by BioJulia's FASTQRecord.
    log_ps = log10(final_probs)
    return DNASequence(final_seq), convert(Vector{Float64}, log_ps)
end

function sample(nseqs::Int, len::Int,
                template_error_rate::Float64,
                t_sub_part::Float64, t_ins_part::Float64, t_del_part::Float64,
                max_error_rate::Float64,
                sub_part::Float64, ins_part::Float64, del_part::Float64;
                error_rate_alpha::Float64=5.0, error_rate_beta::Float64=1.0,
                error_std::Float64=0.01)
    if len % 3 != 0
        error("Reference length must be a multiple of three")
    end

    t_sub_ratio, t_ins_ratio, t_del_ratio = ratios(t_sub_part, t_ins_part, t_del_part)
    reference = random_seq(len)
    template, template_log_p = sample_from_template(reference, template_error_rate,
                                                    t_sub_ratio, t_ins_ratio, t_del_ratio,
                                                    error_std, codon=true)


    error_d = Beta(error_rate_alpha, error_rate_beta)
    sub_ratio, ins_ratio, del_ratio = ratios(sub_part, ins_part, del_part)

    seqs = DNASequence[]
    log_ps = Vector{Float64}[]
    error_rates = Float64[]

    for i = 1:nseqs
        error_rate = Distributions.rand(error_d) * max_error_rate
        seq, log_p = sample_from_template(template, error_rate,
                                          sub_ratio, ins_ratio, del_ratio,
                                          error_std, codon=false)
        push!(seqs, seq)
        push!(log_ps, log_p)
        push!(error_rates, error_rate)
    end
    return DNASequence(reference), DNASequence(template), template_log_p, seqs, log_ps, error_rates
end

"""Write template into FASTA and sequences into FASTQ."""
function write_samples(filename, reference, template, seqs, log_ps)
    write_template(string(filename, "-reference.fasta"), template)
    write_template(string(filename, "-template.fasta"), template)
    write_fastq(string(filename, "-sequences.fastq"), seqs, log_ps)
end

"""Read template from FASTA and sequences from FASTQ."""
function read_samples(filename)
    reference = read_template(string(filename, "-reference.fasta"))
    template = read_template(string(filename, "-template.fasta"))
    seqs, log_ps = read_fastq(string(filename, "-sequences.fastq"))
    return reference, template, seqs, log_ps
end

end
