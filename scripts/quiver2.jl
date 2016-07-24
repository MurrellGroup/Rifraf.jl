import Bio.Seq
@everywhere using Bio.Seq
using ArgParse
using Glob

import Quiver2.Model
import Quiver2.QIO
import Quiver2.Util
@everywhere using Quiver2.Model
@everywhere using Quiver2.QIO
@everywhere using Quiver2.Util

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--reference"
        help = "single file, or a pattern to find reference files; must match inputs after sorting"
        arg_type = AbstractString
        default = ""

        "--prefix"
        help = "prepended to each filename to make label"
        arg_type = AbstractString
        default = ""

        "--keep-unique-name"
        help = "keep only unique middle part of filename"
        action = :store_true

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

        "--indel-fastq"
        help = "if given, a file to store indel base probs"
        arg_type = AbstractString
        default = ""

        "--verbose", "-v"
        help = "print progress"
        arg_type = Int
        default = 0

        "input"
        help = "a single file or a glob. each filename should be unique."
        required = true

        "output"
        help = "output fastq file"
        required = true

    end
    return parse_args(s)
end

@everywhere function dofile(file, reffile, args)
    if args["verbose"] > 0
        println(STDERR, "reading sequences from '$(file)'")
        if length(reffile) > 0
            println(STDERR, "reading reference from '$(reffile)'")
        end
    end
    reference = length(reffile) > 0 ? read_fasta(reffile)[1] : DNASequence("")

    sequences, phreds = read_fastq(file)
    template = sequences[1]
    if args["verbose"] > 0
        println(STDERR, "starting run")
    end
    return quiver2(template, sequences, phreds,
                   reference=reference,
                   bandwidth=args["bandwidth"],
                   min_dist=args["min-dist"],
                   batch=args["batch"],
                   max_iters=args["max-iters"],
                   verbose=args["verbose"])
end

function common_prefix(strings)
    minlen = minimum([length(s) for s in strings])
    cand = strings[1][1:minlen]
    x = 0
    for i = 1:minlen
        if all([s[i] == cand[i] for s in strings])
            x = i
	else
	    break
        end
    end
    return cand[1:x]
end

function common_suffix(strings)
    return reverse(common_prefix([reverse(s) for s in strings]))
end

function main()
    args = parse_commandline()
    input = args["input"]

    dir, pattern = splitdir(input)
    infiles = glob(pattern, dir)
    if length(infiles) == 0
       return
    end
    infiles = sort(infiles)
    names = [splitext(basename(f))[1] for f in infiles]
    if length(Set(names)) != length(names)
        error("Files do not have unique names")
    end

    refpattern = args["reference"]
    reffiles = []
    if length(refpattern) > 0
        if isfile(refpattern)
            reffiles = [refpattern for i in 1:length(infiles)]
        else
            dir, pattern = splitdir(refpattern)
            reffiles = glob(pattern, dir)
            if length(reffiles) != length(infiles)
                error("Wrong number of reference files")
            end
            reffiles = sort(reffiles)
        end
    else
        reffiles = ["" for i in 1:length(infiles)]
    end

    repeated_args = [args for i in 1:length(infiles)]
    results = pmap((f, r, a) -> dofile(f, r, a), infiles, reffiles, repeated_args)

    plen = 0
    slen = 0
    if args["keep-unique-name"]
        plen = length(common_prefix(names))
        snames = [n[plen+1:end] for n in names]
	slen = length(common_suffix(snames))
    end

    n_converged = 0
    prefix = args["prefix"]
    indel_handle = open("/dev/null", "w")
    if length(args["indel-fastq"]) > 0
        indel_handle = open(args["indel-fastq"], "w")
        write(indel_handle, "label",
              ",", "max_indel_p",
              ",", "sum_indel_p",
              ",", "indel_rate",
              "\n")
    end
    stream = open(args["output"], "w", FASTQ, ascii_offset=33)
    for i=1:length(results)
        if typeof(results[i]) == RemoteException
            throw(results[i])
        end
        consensus, base_probs, ins_probs, info = results[i]
        if info["converged"]
            n_converged += 1
            name = names[i]
            if args["keep-unique-name"]
                name = name[plen + 1:end - slen]
            end
            label = string(prefix, name)
            quality = Quiver2.Model.estimate_point_probs(base_probs, ins_probs)
            t_phred = p_to_phred(quality)
            record = Seq.FASTQSeqRecord(label, consensus, t_phred)
            write(stream, record)
            if length(args["indel-fastq"]) > 0
                indel_p = Quiver2.Model.estimate_indel_probs(base_probs, ins_probs)
                max_indel_p = maximum(indel_p)
                sum_indel_p = sum(indel_p)
                indel_rate = sum_indel_p / length(indel_p)
                write(indel_handle, label,
                      ",", string(max_indel_p),
                      ",", string(sum_indel_p),
                      ",", string(indel_rate),
                      "\n")
            end
        end
    end
    close(indel_handle)
    close(stream)
    if args["verbose"] > 0
        println(STDERR, "done. $n_converged / $(length(results)) converged.")
    end
end

main()
