@everywhere using BioSymbols
@everywhere using BioSequences
using ArgParse
using Glob

@everywhere using Rifraf

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        # script-only arguments
        "--phred-cap"
        help = "maximum PHRED score"
        arg_type = Int8
        default = Int8(0)

        "--prefix"
        help = "prepended to each filename to make name"
        arg_type = String
        default = ""

        "--keep-unique-name"
        help = "keep only unique middle part of filename"
        action = :store_true

        # RIFRAF arguments
        "--reference"
        help = string("reference fasta file.",
                      " uses first sequence unless --reference-map is given.")
        arg_type = String
        default = ""

        "--reference-map"
        help = "file mapping input filename to reference id"
        arg_type = String
        default = ""

        "--ref-errors"
        help = "comma-seperated reference error ratios - mm, codon ins, codon del"
        arg_type = String
        default = "10,0.1,0.1,1,1"

        "--max-iters"
        help = "maximum iterations before giving up"
        arg_type = Int
        default = 100

        "--verbose", "-v"
        help = "print progress"
        arg_type = Int
        default = 0

        # required arguments
        "seq-errors"
        help = "comma-seperated sequence error ratios - mismatch, insertion, deletion"
        arg_type = String
        required = true

        "input"
        help = "a single file or a glob. each filename should be unique."
        required = true

        "output"
        help = "output fasta file"
        required = true

    end
    return parse_args(s)
end

@everywhere function dofile(file, reffile, refid, args)
    if args["verbose"] >= 1
        println(STDERR, "reading sequences from '$(file)'")
        if length(reffile) > 0
            println(STDERR, "reading reference from '$(reffile)'")
        end
    end
    reference = DNASeq()
    ref_records = []
    if length(reffile) > 0
        ref_records = Rifraf.read_fasta_records(reffile)
    end
    if length(refid) > 0
        filt_records = collect(filter(r -> FASTA.identifier(r) == refid, ref_records))
        if length(filt_records) == 0
            error("reference '$refid' not found in `$reffile`")
        end
        if length(filt_records) > 1
            error("multiple references with id '$refid' found in `$reffile`")
        end
        ref_record = filt_records[1]
        reference = DNASeq(sequence(ref_record))
    elseif length(ref_records) > 0
        ref_record = ref_records[1]
        reference = DNASeq(sequence(ref_record))
    end

    score_args = map(x -> parse(Float64, x), split(args["seq-errors"], ","))
    scores = Scores(ErrorModel(score_args...))

    ref_score_args = map(x -> parse(Float64, x), split(args["ref-errors"], ","))
    ref_scores = Scores(ErrorModel(ref_score_args...))

    sequences, phreds, _ = Rifraf.read_fastq(file, seqtype=DNASeq)
    if args["verbose"] >= 1
        println(STDERR, "starting run")
    end
    phred_cap = args["phred-cap"]
    if phred_cap > 0
        phreds = Vector{Phred}[Rifraf.cap_phreds(p, phred_cap)
                               for p in phreds]
    end
    params = RifrafParams(scores=scores,
                          ref_scores=ref_scores,
                          max_iters=args["max-iters"],
                          verbose=args["verbose"])
    return rifraf(sequences, phreds;
                  reference=reference,
                  params=params)
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
        if args["verbose"] >= 1
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
                   infiles, refids)

    plen = 0
    slen = 0
    if args["keep-unique-name"]
        plen = length(common_prefix(basenames))
        snames = [n[plen+1:end] for n in basenames]
        slen = length(common_suffix(snames))
    end

    n_converged = 0
    prefix = args["prefix"]
    stream = open(FASTA.Writer, args["output"])
    for (i, result) in enumerate(results)
        if typeof(result) == RemoteException
            throw(result)
        end
        if result.state.converged
            n_converged += 1
            name = basenames[i]
            if args["keep-unique-name"]
                name = name[plen + 1:end - slen]
            end
            seqname = string(prefix, name)
            record = FASTA.Record(seqname, result.consensus)
            write(stream, record)
        end
    end
    close(stream)
    if args["verbose"] >= 1
        println(STDERR, "done. $n_converged / $(length(results)) converged.")
    end
end

main()
