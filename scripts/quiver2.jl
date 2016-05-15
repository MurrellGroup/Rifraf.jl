import Bio.Seq
@everywhere using Bio.Seq
using ArgParse
using Glob

import Quiver2.Model
import Quiver2.QIO
@everywhere using Quiver2.Model
@everywhere using Quiver2.QIO

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--prefix"
        help = "prependedto each filename to make label"
        arg_type = AbstractString
        default = ""

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

        "log_ins"
        help = "log10 insertion probability"
        arg_type = Float64
        required = true

        "log_del"
        help = "log10 deletion probability"
        arg_type = Float64
        required = true

        "input"
        help = "a single file or a glob. each filename should be unique."
        required = true
    end
    return parse_args(s)
end

@everywhere function dofile(file, args)
    if args["verbose"]
        print(STDERR, "reading sequences from '$(file)'\n")
    end
    sequences, log_ps = read_fastq(file)
    template = sequences[1]
    if args["verbose"]
        print(STDERR, "starting run\n")
    end
    consensus, info = quiver2(template, sequences, log_ps,
                              args["log_ins"], args["log_del"],
                              bandwidth=args["bandwidth"], min_dist=args["min-dist"],
                              batch=args["batch"], max_iters=args["max-iters"],
                              verbose=args["verbose"])
    return consensus
end

function main()
    args = parse_commandline()
    input = args["input"]

    dir, pattern = splitdir(input)
    infiles = glob(pattern, dir)
    names = [splitext(basename(f))[1] for f in infiles]
    if length(Set(names)) != length(names)
        error("Files do not have unique names")
    end

    @everywhere function f(x)
        return dofile(x, args)
    end
    results = pmap((f, a) -> dofile(f, a), infiles, [args for i in 1:length(infiles)])

    prefix = args["prefix"]
    for i=1:length(results)
        name = splitext(basename(infiles[i]))[1]
        label = string(prefix, name)
        write(STDOUT, Seq.FASTASeqRecord(label, results[i], Seq.FASTAMetadata("")))
    end
end

main()
