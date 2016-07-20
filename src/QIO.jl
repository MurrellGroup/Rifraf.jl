module QIO

using Bio.Seq

export write_fasta, read_fastq, write_fastq, read_fasta, write_single

function read_fasta(filename)
    stream = open(filename, FASTA)
    seqs = DNASequence[]
    for entry in stream
        push!(seqs, entry.seq)
    end
    return seqs
end

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
    phreds = Vector{UInt8}[]
    for entry in stream
        push!(seqs, entry.seq)
        push!(phreds, entry.metadata.quality)
    end
    return seqs, phreds
end

function write_fastq(filename, seqs, phreds::Vector{Vector{UInt8}})
    stream = open(filename, "w")
    i = 0
    for (s, q) in zip(seqs, phreds)
        i += 1
        write(stream, Seq.FASTQSeqRecord(string("seq_", i), s, Seq.FASTQMetadata("", q)))
    end
    close(stream)
end

function write_single(filename, seq; name::AbstractString="seq")
    record = Seq.FASTASeqRecord(name, seq, Seq.FASTAMetadata(""))
    fasta_stream = open(filename, "w")
    write(fasta_stream, record)
    close(fasta_stream)
end

end
