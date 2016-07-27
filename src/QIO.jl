module QIO

using Bio.Seq

export read_fasta_records, read_fasta, write_fasta, read_fastq_records, read_fastq, write_fastq

function read_fasta_records(filename)
    stream = open(filename, FASTA)
    records = FASTASeqRecord[]
    for entry in stream
        push!(records, entry)
    end
    return records
end

function read_fasta(filename)
    records = read_fastq_records(filename)
    return DNASequence[r.seq for r in records]
end

function write_fasta(filename, seqs)
    stream = open(filename, "w", FASTA)
    for i = 1:length(seqs)
        s = seqs[i]
        write(stream, Seq.FASTASeqRecord(string("seq_", i), s))
    end
    close(stream)
end

function read_fastq_records(filename)
    stream = open(filename, "r", FASTQ, quality_encoding=Seq.SANGER_QUAL_ENCODING)
    records = FASTQSeqRecord[]
    for entry in stream
        # BioJulia gets the offset wrong
        phred = entry.metadata.quality + 64 - 33
        if minimum(phred) < Int('!') - 33 || maximum(phred) > Int('~') - 33
            error("quality scores  not in range")
        end
        entry.metadata.quality = phred
        push!(records, entry)
    end
    return records
end

function read_fastq(filename)
    records = read_fastq_records(filename)
    seqs = DNASequence[]
    phreds = Vector{Int8}[]
    for record in records
        push!(seqs, record.seq)
        push!(phreds, record.metadata.quality)
    end
    return seqs, phreds
end

function write_fastq(filename, seqs, phreds::Vector{Vector{Int8}})
    stream = open(filename, "w", FASTQ, ascii_offset=33)
    i = 0
    for (s, q) in zip(seqs, phreds)
        i += 1
        write(stream, Seq.FASTQSeqRecord(string("seq_", i), s, q))
    end
    close(stream)
end

end
