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
        help = string("reference fasta file.",
                      " uses first sequence unless --reference-map is given.")
        arg_type = AbstractString
        default = ""

        "--reference-map"
        help = "file mapping input filename to reference id"
        arg_type = AbstractString
        default = ""

        "--prefix"
        help = "prepended to each filename to make name"
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

        "--indel-file"
        help = "if given, a file to store indel probabilities"
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

@everywhere function dofile(file, reffile, refid, args)
    if args["verbose"] > 0
        println(STDERR, "reading sequences from '$(file)'")
        if length(reffile) > 0
            println(STDERR, "reading reference from '$(reffile)'")
        end
    end
    reference = DNASequence("")
    ref_records = []
    if length(reffile) > 0
        ref_records = read_fasta_records(reffile)
    end
    if length(refid) > 0
        ref_record = collect(filter(r -> r.name == refid, ref_records))[1]
        reference = ref_record.seq
    elseif length(ref_records) > 0
        ref_record = ref_records[1]
        reference = ref_record.seq
    end

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
        if args["verbose"] > 0
            println(STDERR, "warning: no input files found.")
        end
        return
    end
    infiles = sort(infiles)
    basenames = [basename(f) for f in infiles]
    if length(Set(basenames)) != length(basenames)
        error("Files do not have unique names")
    end

    if length(args["reference"]) > 0
        if !isfile(args["reference"])
            error("reference file not found")
        end
        if length(args["reference-map"]) > 0
            if !isfile(args["reference-map"])
                error("reference map file not found")
            end
        end
    else
        if length(args["reference-map"]) > 0
            error("--reference-map is invalid without --reference")
        end
    end

    reffile = args["reference"]
    refmapfile = args["reference-map"]
    refids = ["" for i in 1:length(infiles)]
    if length(refmapfile) > 0
        refmap_handle = open(refmapfile)
        name_to_ref = Dict()
        for line in readlines(refmap_handle)
            name, refid = split(line)
            name_to_ref[name] = refid
        end
        close(refmap_handle)
        # only process seqs with a reference
        present_names = Set(keys(name_to_ref))
        infiles = sort(collect(filter(f -> basename(f) in present_names, infiles)))
        basenames = [basename(f) for f in infiles]
        refids = [name_to_ref[name] for name in basenames]
    end
    results = pmap((f, rid) -> dofile(f, reffile, rid, args),
                   infiles, refids, err_stop=true)

    plen = 0
    slen = 0
    if args["keep-unique-name"]
        plen = length(common_prefix(basenames))
        snames = [n[plen+1:end] for n in basenames]
        slen = length(common_suffix(snames))
    end

    n_converged = 0
    prefix = args["prefix"]
    indel_handle = open("/dev/null", "w")
    if length(args["indel-file"]) > 0
        indel_handle = open(args["indel-file"], "w")
        write(indel_handle, "name",
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
            name = basenames[i]
            if args["keep-unique-name"]
                name = name[plen + 1:end - slen]
            end
            seqname = string(prefix, name)
            quality = Quiver2.Model.estimate_point_probs(base_probs, ins_probs)
            t_phred = p_to_phred(quality)
            record = Seq.FASTQSeqRecord(seqname, consensus, t_phred)
            write(stream, record)
            if length(args["indel-file"]) > 0
                indel_p = Quiver2.Model.estimate_indel_probs(base_probs, ins_probs)
                max_indel_p = maximum(indel_p)
                sum_indel_p = sum(indel_p)
                indel_rate = sum_indel_p / length(indel_p)
                write(indel_handle, seqname,
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
