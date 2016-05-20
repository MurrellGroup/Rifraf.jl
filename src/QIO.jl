module QIO

using Bio.Seq

export write_fasta, read_fastq, write_fastq, read_single, write_single

function write_fasta(filename, seqs)
    stream = open(filename, "w")
    for i = 1:length(seqs)
        s = seqs[i]
        write(stream, Seq.FASTASeqRecord(string("seq_", i), s, Seq.FASTAMetadata("")))
    end
    close(stream)
end

function read_fastq(filename)
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

function write_fastq(filename, seqs, log_ps)
    stream = open(filename, "w")
    i = 0
    for (s, lps) in zip(seqs, log_ps)
        i += 1
        q = UInt8[UInt8(min(round(-10 * p), 80)) for p in lps]
        write(stream, Seq.FASTQSeqRecord(string("seq_", i), s, Seq.FASTQMetadata("", q)))
    end
    close(stream)
end

function read_single(filename)
    stream = open(filename, FASTA)
    entry = Seq.FASTADNASeqRecord()
    while read!(stream, entry)
        return entry.seq
    end
end

function write_single(filename, seq; name::AbstractString="seq")
    record = Seq.FASTASeqRecord(name, seq, Seq.FASTAMetadata(""))
    fasta_stream = open(filename, "w")
    write(fasta_stream, record)
    close(fasta_stream)
end

end
