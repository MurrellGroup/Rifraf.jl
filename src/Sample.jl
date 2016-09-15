module Sample

using Bio.Seq
using Distributions
using Iterators

using Quiver2.QIO
using Quiver2.Util

export rbase, mutate_base, random_codon, random_seq, sample_from_reference, sample_from_template, sample, BetaAlt

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

function random_codon()
    return (rbase(), rbase(), rbase())
end

function random_seq(n)
    return DNASequence([rbase() for i in 1:n])
end

function normalize(x, lower, upper)
    return (x - lower) / (upper - lower)
end

function restore(y, lower, upper)
    return y * (upper - lower) + lower
end


const MIN_PROB = 1e-10
const MAX_PROB = 0.5


"""A beta distribution parametrized by mean and standard deviation."""
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
        std_scaled = std / (maxval - minval)
        if std_scaled <= 0
            error("std must be > 0")
        end
        variance = std_scaled ^ 2
        max_variance = mean_scaled * (1 - mean_scaled)
        if variance >= max_variance
            error("variance cannot be larger than mean * (1 - mean)")
        end
        nu = mean_scaled * (1 - mean_scaled) / variance - 1.0
        alpha = mean_scaled * nu
        beta = (1.0 - mean_scaled) * nu
        d = Beta(alpha, beta)
        new(mean, std, minval, maxval, alpha, beta, d)
    end
end

function rand(s::BetaAlt)
    return restore(Distributions.rand(s.d), s.minval, s.maxval)
end

"""
Add independent log-gaussian noise to each position in vector.

Keep within range [0, 1].
"""
function jitter_vector(x::Vector{Float64},
                       log_std::Float64)
    error = randn(length(x)) * log_std
    result = exp10(log10(x) + error)
    result[map(a -> a < MIN_PROB, result)] = MIN_PROB
    result[map(a -> a > MAX_PROB, result)] = MAX_PROB
    return result
end


function ratios(sub_part, ins_part, del_part)
    denom = sub_part + ins_part + del_part
    sub_ratio = sub_part / denom
    ins_ratio = ins_part / denom
    del_ratio = del_part / denom
    return sub_ratio, ins_ratio, del_ratio
end


function substitute(base, p::Float64)
    if Base.rand(Bernoulli(x)) == 1
        return mutate_base(base)
    end
    return base
end


function hmm_sample(sequence::DNASequence,
                    error_p::Vector{Float64},
                    error_ratios::Tuple{Float64, Float64, Float64};
                    codon=false)
    if codon && length(sequence) % 3 != 0
        error("sequence length is not multiple of 3")
    end
    sub_ratio, ins_ratio, del_ratio = ratios(error_ratios...)
    final_seq = []
    final_error_p = Float64[]
    skip = 0
    for i = 1:(length(sequence) + 1)
        p = (i > length(sequence) ? error_p[i - 1] : error_p[i])
        prev_p = (i == 1 ? error_p[1] : error_p[i - 1])
        # insertion between i-1 and i
        mean_p = mean([p, prev_p])
        ins_p = mean_p * ins_ratio
        if codon
            ins_p /= 3.0
        end
        if Base.rand(Bernoulli(ins_p)) == 1
            if codon
                push!(final_seq, random_codon()...)
                push!(final_error_p, collect(repeated(mean_p / ins_ratio, 3))...)
            else
                push!(final_seq, rbase())
                push!(final_error_p, mean_p)
            end
        end
        if i > length(sequence)
            break
        end

        # only skip after insertions, to ensure equal probability of
        # insertion and deletions
        if skip > 0
            skip -= 1
            continue
        end
        # deletion of i
        del_p = p * del_ratio
        if codon
            if i > length(sequence) - 2
                del_p = 0.0
            else
                del_p = mean(error_p[i:i+2]) * del_ratio / 3
            end
        end
        if Base.rand(Bernoulli(del_p)) == 1
            # skip position i, and possibly the entire codon starting at i
            skip = codon ? 2 : 0
        else
            # mutation of position i
            if Base.rand(Bernoulli(p * sub_ratio)) == 1
                push!(final_seq, mutate_base(sequence[i]))
            else
                push!(final_seq, sequence[i])
            end
            push!(final_error_p, p)
        end
    end
    return DNASequence(final_seq), final_error_p
end


function sample_from_reference(reference::DNASequence,
                               error_rate::Float64,
                               error_ratios::Tuple{Float64, Float64, Float64})
    error_p = error_rate *  ones(length(reference))
    template, e = hmm_sample(reference, error_p, error_ratios, codon=true)
    return template
end


"""
error_ratios: (sub, ins, del)
template_error_p: vector of per-base error rates
actual_error_std: standard deviation of Beta distribution
    for actual errors
reported_error_std: standard deviation of Beta distribution
    for reported erros
codon: if true, only do codon indels

"""
function sample_from_template(template::DNASequence,
                              template_error_p::Vector{Float64},
                              error_ratios::Tuple{Float64, Float64, Float64},
                              log_actual_error_std::Float64,
                              log_reported_error_std::Float64)
    # add noise to simulate measurement error
    jittered_error_p = jitter_vector(template_error_p,
                                     log_actual_error_std)

    seq, actual_error_p = hmm_sample(template, jittered_error_p,
                                     error_ratios, codon=false)

    # add noise to simulate quality score estimation error
    reported_error_p = jitter_vector(actual_error_p,
                                     log_reported_error_std)
    phreds = p_to_phred(reported_error_p)
    return DNASequence(seq), actual_error_p, phreds
end

"""
nseqs: number of sequences to sample
len: length of the reference
template_error_rate: overall error rate sampling template
    from reference
template_error_ratios: (sub, ins, del) error ratios for
    sampling template from reference
template_error_mean: mean of Beta distribution for drawing
    per-base template sequencing error rates
template_error_std: standard deviation of same Beta
log_seq_actual_std: standard deviation for jittering actual
    sequence per-base error rate (measurement noise)
log_seq_quality_std: standard deviation for jittering reported
    sequence per-base error rate (quality score estimation noise)
seq_error_ratios: (sub, ins, del) sequence error ratios

"""
function sample(nseqs::Int, len::Int,
                template_error_rate::Float64,
                template_error_ratios::Tuple{Float64, Float64, Float64},
                template_error_mean::Float64,
                template_error_std::Float64,
                log_seq_actual_std::Float64,
                log_seq_reported_std::Float64,
                seq_error_ratios::Tuple{Float64, Float64, Float64})
    if len % 3 != 0
        error("Reference length must be a multiple of three")
    end

    reference = random_seq(len)
    template = sample_from_reference(reference,
                                     template_error_rate,
                                     template_error_ratios)
    dist = BetaAlt(template_error_mean, template_error_std,
                   minval=MIN_PROB, maxval=MAX_PROB)
    template_error_p = Float64[rand(dist) for i = 1:length(template)]

    # left off here
    seqs = DNASequence[]
    actual_error_ps = Vector{Float64}[]
    phreds = Vector{Int8}[]

    for i = 1:nseqs
        (seq, actual_error_p, phred) = sample_from_template(template,
                                                            template_error_p,
                                                            seq_error_ratios,
                                                            log_seq_actual_std,
                                                            log_seq_reported_std)
        push!(seqs, seq)
        push!(actual_error_ps, actual_error_p)
        push!(phreds, phred)
    end
    return (DNASequence(reference), DNASequence(template),
            template_error_p, seqs, actual_error_ps, phreds)
end

"""Write template into FASTA and sequences into FASTQ."""
function write_samples(filename, reference, template, template_error, seqs, phreds)
    template_phred = p_to_phred(template_error)
    write_fasta(string(filename, "-reference.fasta"), [reference])
    write_fastq(string(filename, "-template.fastq"), [template], Vector{Int8}[template_phred])
    write_fastq(string(filename, "-sequences.fastq"), seqs, phreds)
end

"""Read template from FASTA and sequences from FASTQ."""
function read_samples(filename)
    reference = read_fasta(string(filename, "-reference.fasta"))[1]
    template_seqs, template_phreds = read_fastq(string(filename, "-template.fastq"))
    template = template_seqs[1]
    template_error = phred_to_p(template_phreds[1])
    seqs, phreds = read_fastq(string(filename, "-sequences.fastq"))
    return reference, template, template_error, seqs, phreds
end

end
