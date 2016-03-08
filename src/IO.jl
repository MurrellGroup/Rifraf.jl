__precompile__()
module IO

using Bio.Seq

export read_sequences, write_sequences, read_template, write_template

function read_sequences(filename)
    stream = open(filename, FASTQ)
    seqs = DNASequence[]
    log_ps = Vector{Float64}[]
    for entry in stream
        push!(seqs, entry.seq)
        lps = Float64[Float64(q) / (-10.0) for q in entry.metadata.quality]
        push!(log_ps, lps)
    end
    return seqs, log_ps
end

function write_sequences(filename, seqs, log_ps)
    stream = open(filename, "w")
    i = 0
    for (s, lps) in zip(seqs, log_ps)
        i += 1
        q = UInt8[UInt8(min(round(-10 * p), 80)) for p in lps]
        write(stream, Seq.FASTQSeqRecord(string("seq_", i), s, Seq.FASTQMetadata("", q)))
    end
    close(stream)
end

function read_template(filename)
    stream = open(filename, FASTA)
    for entry in stream
        return entry.seq
    end
end

function write_template(filename, template)
    template_record = Seq.FASTASeqRecord("template", template, Seq.FASTAMetadata(""))
    fasta_stream = open(filename, "w")
    write(fasta_stream, template_record)
    close(fasta_stream)
end

end
