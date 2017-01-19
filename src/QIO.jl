module QIO

using Bio.Seq

export read_fasta_records, read_fasta, write_fasta, read_fastq_records, read_fastq, write_fastq

function read_fasta_records(filename)
    stream = open(FASTAReader, filename)
    records = FASTASeqRecord[]
    for entry in stream
        push!(records, entry)
    end
    return records
end

function read_fasta(filename)
    records = read_fasta_records(filename)
    return DNASequence[r.seq for r in records]
end

function write_fasta(filename, seqs)
    stream = open(FASTAWriter, filename)
    for i = 1:length(seqs)
        s = seqs[i]
        write(stream, Seq.FASTASeqRecord(string("seq_", i), s))
    end
    close(stream)
end

function read_fastq_records(filename)
    stream = open(FASTQReader, filename, quality_encoding=:sanger)
    records = FASTQSeqRecord[]
    for record in stream
        if any([q < 0 for q in record.metadata.quality])
            error("$(record.name) in $filename contains negative phred values")
        end
        push!(records, record)
    end
    return records
end

function read_fastq(filename)
    records = read_fastq_records(filename)
    seqs = DNASequence[]
    phreds = Vector{Int8}[]
    names = String[]
    for record in records
        push!(seqs, record.seq)
        push!(phreds, record.metadata.quality)
        push!(names, record.name)
    end
    return seqs, phreds, names
end

function write_fastq(filename, seqs, phreds::Vector{Vector{Int8}};
                     names=String[])
    stream = open(FASTQWriter, filename, quality_encoding=:sanger)
    i = 0
    if length(names) != length(seqs)
        names = [string("seq_", i) for i in 1:length(seqs)]
    end
    for (s, q, n) in zip(seqs, phreds, names)
        i += 1
        write(stream, Seq.FASTQSeqRecord(n, s, q))
    end
    close(stream)
end

end
