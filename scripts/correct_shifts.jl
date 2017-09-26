using BioSymbols
using BioSequences
using ArgParse
using Glob

using Rifraf

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--multi-reference"
        help = "each sequence is followed by its reference"
        action = :store_true

        "--log-p"
        help = "log error probability"
        arg_type = Float64
        default = -1.0

        "--bandwidth"
        help = "alignment bandwidth. If < 0, choose dynamically."
        arg_type = Int
        default = -1

        "--errors"
        help = "comma-seperated reference error ratios - mm, codon ins, codon del"
        arg_type = String
        default = "10,0.00001,0.00001,1,1"

        "--verbose", "-v"
        help = "print progress"
        arg_type = Int
        default = 0

        # required arguments
        "input"
        help = "input fasta file, with sequence/reference alternating pairs"
        required = true

        "output"
        help = "output fasta file of corrected sequences"
        required = true
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()

    score_args = map(x -> parse(Float64, x), split(args["errors"], ","))
    scores = Scores(ErrorModel(score_args...))

    input = args["input"]
    input_sequences = Rifraf.read_fasta_records(input)

    multi = args["multi-reference"]
    sequences = multi ? input_sequences[1:2:end-1] : input_sequences[2:end]
    references = multi ? input_sequences[2:2:end] : repeated(input_sequences[1])

    stream = open(FASTAWriter, args["output"])
    for (seq, ref) in zip(sequences, references)
        result = Rifraf.correct_shifts(DNASeq(seq.seq), DNASeq(ref.seq);
                                       log_p=args["log-p"],
                                       bandwidth=args["bandwidth"],
                                       scores=scores)
        record = Seq.FASTASeqRecord(seq.name, result)
        write(stream, record)
    end
    close(stream)
end

main()
