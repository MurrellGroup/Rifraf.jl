module Sample

export rbase, mutate_base, random_seq, coinflip, sample_from_template, phred, sample

function rbase()
    bases = ['A', 'C', 'G', 'T']
    return bases[rand(1:4)]
end

function mutate_base(base::Char)
    result = rbase()
    while result == base
        result = rbase()
    end
    return result
end

function random_seq(n)
    return join([rbase() for i in 1:n])
end

function coinflip(p)
    return rand() < p
end

function sample_from_template(template, sub_rate, insertion_rate, deletion_rate)
    result = []
    for base in template
        while coinflip(insertion_rate)
            push!(result, rbase())
        end
        if coinflip(deletion_rate)
            continue
        end
        if coinflip(sub_rate)
            push!(result, mutate_base(base))
        else
            push!(result, base)
        end
    end
    while coinflip(insertion_rate)
        push!(result, rbase())
    end
    return join(result)
end

function phred(p)
    return -10 * log10(p)
end

function sample(n::Int, len::Int, sub_rate::Float64, ins_rate::Float64, del_rate::Float64)
    reads = ASCIIString[]
    phreds = Vector{Float64}[]
    template = random_seq(len)
    error_rate = sub_rate + ins_rate + del_rate
    for i in 1:n
        r = sample_from_template(template, sub_rate, ins_rate, del_rate)
        push!(phreds, phred(fill(error_rate, length(r))))
        push!(reads, r)
    end
    return template, reads, phreds
end

function tofasta(filename::ASCIIString, sequences::Vector{ASCIIString})
    outfile = open(filename, "w")
    for i in 1:length(sequences)
        s = sequences[i]
        write(outfile, ">seq_$(i)\n")
        write(outfile, "$s\n")
    end
    close(outfile)
end

end
