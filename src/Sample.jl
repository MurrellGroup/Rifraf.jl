module Sample

using Bio.Seq
using Distributions
using Iterators

using Quiver2.QIO
using Quiver2.Util
using Quiver2.Model

export rbase, mutate_base, random_codon
export random_seq, sample_reference
export sample_from_template, sample, BetaAlt
export sample_mixture

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

function mutate_seq(seq::DNASequence, n_diffs::Int)
    seq = copy(seq)
    positions = Base.rand(1:length(seq), n_diffs)
    for i in positions
        seq[i] = mutate_base(seq[i])
    end
    return seq
end


function random_codon()
    return (rbase(), rbase(), rbase())
end

function random_seq(n)
    return DNASequence([rbase() for i in 1:n])
end

function rescale(x, lower, upper)
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
        mean_scaled = rescale(mean, minval, maxval)
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


function substitute(base, p::Float64)
    if Base.rand(Bernoulli(x)) == 1
        return mutate_base(base)
    end
    return base
end


function hmm_sample(sequence::DNASequence,
                    error_p::Vector{Float64},
                    errors::ErrorModel)
    errors = Model.normalize(errors)
    codon = errors.codon_insertion > 0.0 || errors.codon_deletion > 0.0
    if codon && (errors.insertion > 0.0 || errors.deletion > 0.0)
        error("codon and non-codon indels are not both allowed")
    end
    if codon && length(sequence) % 3 != 0
        error("sequence length is not multiple of 3")
    end
    sub_ratio = errors.mismatch
    ins_ratio = errors.insertion
    del_ratio = errors.deletion
    if codon
        ins_ratio = errors.codon_insertion
        del_ratio = errors.codon_deletion
    end
    final_seq = []
    final_error_p = Float64[]
    seqbools = Bool[]
    tbools = Bool[]
    skip = 0
    for i = 1:(length(sequence) + 1)
        p = (i > length(sequence) ? error_p[i - 1] : error_p[i])
        prev_p = (i == 1 ? error_p[1] : error_p[i - 1])
        # insertion between i-1 and i
        max_p = max(p, prev_p)
        ins_p = max_p * ins_ratio
        if codon
            ins_p /= 3.0
        end
        if Base.rand(Bernoulli(ins_p)) == 1
            if codon
                push!(final_seq, random_codon()...)
                push!(final_error_p, collect(repeated(max_p, 3))...)
                push!(seqbools, false, false, false)
            else
                push!(final_seq, rbase())
                push!(final_error_p, max_p)
                push!(seqbools, false)
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
        if codon
            if i > length(sequence) - 2
                del_p = 0.0
            else
                del_p = maximum(error_p[i:i+2]) * del_ratio / 3.0
            end
        else
            del_p = p * del_ratio
        end
        if Base.rand(Bernoulli(del_p)) == 1
            # skip position i, and possibly the entire codon starting at i
            skip = codon ? 2 : 0
            append!(tbools, fill(false, skip + 1))
        else
            # mutation of position i
            if Base.rand(Bernoulli(p * sub_ratio)) == 1
                push!(final_seq, mutate_base(sequence[i]))
                push!(seqbools, false)
                push!(tbools, false)
            else
                push!(final_seq, sequence[i])
                push!(seqbools, true)
                push!(tbools, true)
            end
            push!(final_error_p, p)
        end
    end
    return DNASequence(final_seq), final_error_p, seqbools, tbools
end


function sample_reference(reference::DNASequence,
                          error_rate::Float64,
                          errors::ErrorModel)
    errors = Model.normalize(errors)
    if errors.insertion > 0.0 || errors.deletion > 0.0
        error("non-codon indels are not allowed in template")
    end
    error_p = error_rate * ones(length(reference))
    template, e, sbools, tbools = hmm_sample(reference, error_p, errors)
    return template
end


"""
template_error_p: vector of per-base error rates
actual_error_std: standard deviation of Beta distribution
    for actual errors
reported_error_std: standard deviation of Beta distribution
    for reported erros
codon: if true, only do codon indels

"""
function sample_from_template(template::DNASequence,
                              template_error_p::Vector{Float64},
                              errors::ErrorModel,
                              log_actual_error_std::Float64,
                              log_reported_error_std::Float64)
    errors = Model.normalize(errors)
    if errors.codon_insertion > 0.0 || errors.codon_deletion > 0.0
        error("codon indels are not allowed in sequences")
    end
    # add noise to simulate measurement error
    jittered_error_p = jitter_vector(template_error_p,
                                     log_actual_error_std)

    seq, actual_error_p, sbools, tbools = hmm_sample(template, jittered_error_p,
                                                     errors)

    # add noise to simulate quality score estimation error
    reported_error_p = jitter_vector(actual_error_p,
                                     log_reported_error_std)
    phreds = p_to_phred(reported_error_p)
    return DNASequence(seq), actual_error_p, phreds, sbools, tbools
end


function sample_mixture(nseqs::Tuple{Int, Int}, len::Int,
                        n_diffs::Int,
                        ref_error_rate::Float64,
                        ref_errors::ErrorModel,
                        template_error_mean::Float64,
                        template_error_std::Float64,
                        log_seq_actual_std::Float64,
                        log_seq_reported_std::Float64,
                        seq_errors::ErrorModel)
    if len % 3 != 0
        error("Reference length must be a multiple of three")
    end

    template1 = random_seq(len)
    template2 = mutate_seq(template1, n_diffs)
    templates = DNASequence[template1, template2]

    reference = sample_reference(template1,
                                 ref_error_rate,
                                 ref_errors)
    dist = BetaAlt(template_error_mean, template_error_std,
                   minval=MIN_PROB, maxval=MAX_PROB)
    template_error_p = Float64[rand(dist) for i = 1:len]

    seqs = DNASequence[]
    actual_error_ps = Vector{Float64}[]
    phreds = Vector{Int8}[]
    seqbools = Vector{Bool}[]
    tbools = Vector{Bool}[]
    for (t, n) in zip(templates, nseqs)
        for i = 1:n
            (seq, actual_error_p, phred, cb,
             db) = sample_from_template(t,
                                        template_error_p,
                                        seq_errors,
                                        log_seq_actual_std,
                                        log_seq_reported_std)
            push!(seqs, seq)
            push!(actual_error_ps, actual_error_p)
            push!(phreds, phred)
            push!(seqbools, cb)
            push!(tbools, db)
        end
    end
    return (DNASequence(reference),
            DNASequence[DNASequence(t) for t in templates],
            template_error_p, seqs, actual_error_ps, phreds,
            seqbools, tbools)

end

"""
nseqs: number of sequences to sample
len: length of the reference

ref_error_rate: overall error rate sampling template
    from reference
ref_errors: error model for sampling template from reference

template_error_mean: mean of Beta distribution for drawing
    per-base template sequencing error rates
template_error_std: standard deviation of same Beta

log_seq_actual_std: standard deviation for jittering actual
    sequence per-base log error rate (measurement noise)
log_seq_quality_std: standard deviation for jittering reported
    sequence per-base log error rate (quality score estimation noise)
seq_error_ratios: sequence error model

"""
function sample(nseqs::Int, len::Int,
                ref_error_rate::Float64,
                ref_errors::ErrorModel,
                template_error_mean::Float64,
                template_error_std::Float64,
                log_seq_actual_std::Float64,
                log_seq_reported_std::Float64,
                seq_errors::ErrorModel)
    (ref, templates, t_p, seqs, actual,
     phreds, cb, db) = sample_mixture((nseqs, 0), len, 0,
                              ref_error_rate,
                              ref_errors,
                              template_error_mean,
                              template_error_std,
                              log_seq_actual_std,
                              log_seq_reported_std,
                              seq_errors)
    return ref, templates[1], t_p, seqs, actual, phreds, cb, db
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
