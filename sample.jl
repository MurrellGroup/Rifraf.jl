module Sample

using Bio.Seq

using Distributions

export rbase, mutate_base, random_seq, sample_from_template, phred, sample

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

function phred(p)
    return Int8(min(typemax(Int8), floor(-10 * log10(p))))
end

function BetaAlt(mean::Float64, sample_size::Float64)
    alpha = sample_size * mean
    beta = alpha / mean - alpha
    return Beta(alpha, beta)
end

function sample_from_template(template, error_rate,
                              sub_ratio, ins_ratio, del_ratio)
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
    prob_d = BetaAlt(error_rate, convert(Float64, length(seq)))
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

    return Seq.FASTQSeqRecord("",
                              DNASequence(final_seq),
                              Seq.FASTQMetadata("", map(phred, final_probs)))
end

function sample(nseqs::Int, len::Int, max_error_rate::Float64, mean_error_rate::Float64,
                sub_part::Float64, ins_part::Float64, del_part::Float64)
    if mean_error_rate > max_error_rate
        error("max_error_rate must be greater than mean_error_rate")
    end
    error_d = BetaAlt(mean_error_rate / max_error_rate, convert(Float64, nseqs))
    template = random_seq(len)

    sub_ratio = (sub_part / (sub_part + ins_part + del_part))
    ins_ratio = (ins_part / (sub_part + ins_part + del_part))
    del_ratio = (del_part / (sub_part + ins_part + del_part))

    result = Seq.FASTQSeqRecord[]

    for i = 1:nseqs
        error_rate = rand(error_d) * max_error_rate
        record = sample_from_template(template, error_rate,
                                      sub_ratio, ins_ratio, del_ratio)
        push!(result, record)
    end
    return DNASequence(template), result
end

end
