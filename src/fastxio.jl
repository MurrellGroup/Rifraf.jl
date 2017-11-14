"""
    read_fasta_records(filename)

Read a FASTA file and return records.

# Returns:
- `records::Vector{FASTA.Record}`

"""
function read_fasta_records(filename)
    stream = open(FASTA.Reader, filename)
    records = FASTA.Record[]
    for entry in stream
        push!(records, entry)
    end
    return records
end

"""
    read_fasta(filename)

Read a FASTA file and convert to a given sequence type.

# Returns:
- `seqs::Vector{T}`

"""
function read_fasta(filename; seqtype=DNASequence)
    records = read_fasta_records(filename)
    return seqtype[sequence(r) for r in records]
end

"""
    write_fasta(filename, seqs; names)

Write sequences to a FASTA file.

# Arguments:
- `filename`: file into which to write
- `seqs`: sequences to write
- `names::Vector{String}`: optional list of corresponding names

"""
function write_fasta(filename, seqs; names=String[])
    if length(names) == 0
        names = ["seq_$i" for i in 1:length(seqs)]
    end
    stream = open(FASTA.Writer, filename)
    for i = 1:length(seqs)
        write(stream, FASTA.Record(names[i], seqs[i]))
    end
    close(stream)
end

"""
    read_fastq_records(filename)

Read a FASTQ file and return records.

# Returns:
- `records::Vector{FASTQ.Record}`

"""
function read_fastq_records(filename)
    stream = open(FASTQ.Reader, filename)
    records = FASTQ.Record[]
    for record in stream
        if any([q < 0 for q in FASTQ.quality(record)])
            error("$(FASTQ.identifier(record)) in $filename contains negative phred values")
        end
        push!(records, record)
    end
    return records
end

"""
    read_fastq(filename)

Read a FASTQ file and convert to a given sequence type.

# Returns:
- `seqs::Vector{T}`:
- `phreds::Vector{Vector{Phred}}`: Phred values
- `names::Vector{String}`: sequence names

"""
function read_fastq(filename; seqtype=DNASequence)
    records = read_fastq_records(filename)
    seqs = seqtype[]
    phreds = Vector{Phred}[]
    names = String[]
    for record in records
        push!(seqs, sequence(record))
        push!(phreds, FASTQ.quality(record))
        push!(names, FASTQ.identifier(record))
    end
    return seqs, phreds, names
end

"""
    write_fastq(filename, seqs, phreds; names)

Write sequences to a FASTA file.

# Arguments:
- `filename`: file into which to write
- `seqs`: sequences to write
- `phreds`: corresponding Phred scores
- `names::Vector{String}`: optional list of corresponding names

"""
function write_fastq(filename, seqs, phreds::Vector{Vector{Phred}};
                     names=String[])
    stream = open(FASTQ.Writer, filename)
    i = 0
    if length(names) != length(seqs)
        names = [string("seq_", i) for i in 1:length(seqs)]
    end
    for (s, q, n) in zip(seqs, phreds, names)
        i += 1
        write(stream, FASTQ.Record(n, s, q))
    end
    close(stream)
end
