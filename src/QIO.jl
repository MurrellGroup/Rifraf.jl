module QIO

using Bio.Seq

export read_fasta, write_fasta, read_fastq, write_fastq

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
        write(stream, Seq.FASTASeqRecord(string("seq_", i), s))
    end
    close(stream)
end

function read_fastq(filename)
    stream = open(filename, FASTQ)
    seqs = DNASequence[]
    phreds = Vector{Int8}[]
    for entry in stream
        push!(seqs, entry.seq)
        push!(phreds, entry.metadata.quality)
    end
    return seqs, phreds
end

function write_fastq(filename, seqs, phreds::Vector{Vector{Int8}})
    stream = open(filename, "w")
    i = 0
    for (s, q) in zip(seqs, phreds)
        i += 1
        write(stream, Seq.FASTQSeqRecord(string("seq_", i), s, q))
    end
    close(stream)
end

end
