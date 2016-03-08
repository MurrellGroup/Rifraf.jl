using Bio.Seq
using ArgParse

using Quiver2.Model
using Quiver2.IO

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--template"
        help = "FASTA file containing a template."

        "--name"
        help = "name to give consensus sequence"
        arg_type = ASCIIString
        default = "consensus"

        "--bandwidth"
        help = "alignment bandwidth"
        arg_type = Int
        default = 10

        "--min-dist"
        help = "minimum distance between mutations"
        arg_type = Int
        default = 9

        "--batch"
        help = "batch size; -1 for no batch iterations"
        arg_type = Int
        default = 10

        "--max-iters"
        help = "maximum iterations before giving up"
        arg_type = Int
        default = 100

        "--verbose", "-v"
        help = "print progress"
        action = :store_true

        "infile"
        help = "input FASTQ file"
        required = true

        "log_ins"
        help = "log10 insertion probability"
        arg_type = Float64
        required = true

        "log_del"
        help = "log10 deletion probability"
        arg_type = Float64
        required = true
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()

    infile = args["infile"]
    if args["verbose"]
        print(STDERR, "reading sequences from '$infile'\n")
    end
    sequences, log_ps = read_sequences(args["infile"])

    template = sequences[1]
    tfile = args["template"]
    if tfile != nothing
        if args["verbose"]
            print(STDERR, "reading template from '$tfile'\n")
        end
        template = read_template(args["template"])
    end

    if args["verbose"]
        print(STDERR, "starting run\n")
    end
    consensus, info = quiver2(template, sequences, log_ps,
                              args["log_ins"], args["log_del"],
                              bandwidth=args["bandwidth"], min_dist=args["min-dist"],
                              batch=args["batch"], max_iters=args["max-iters"],
                              verbose=args["verbose"])

    write(STDOUT, Seq.FASTASeqRecord(args["name"], consensus, Seq.FASTAMetadata("")))

end

main()
