module Rifraf

using BioSequences
using Distributions
using Parameters
using StatsBase

import Base.length
import Base.reverse


export rifraf,
       RifrafParams,
       RifrafResult,
       ErrorModel,
       Scores,
       normalize,
       sample_sequences,
       write_samples,
       read_samples,
       sample_mixture,
       RifrafSequence,
       length,
       reverse,
       DNASeq,
       Phred,
       Prob,
       LogProb,
       Score,
       BandedArray,
       sparsecol,
       data_row,
       row_range,
       data_row_range,
       inband,
       flip!,
       newbandwidth!,
       align,
       calibrate_phreds,
       read_fastq_records,
       read_fastq,
       write_fastq,
       read_fasta_records,
       read_fasta,
       write_fasta

include("bandedarrays.jl")
include("types.jl")
include("phred.jl")
include("util.jl")
include("fastxio.jl")
include("proposals.jl")
include("errormodel.jl")
include("rifrafsequences.jl")
include("align.jl")
include("model.jl")
include("sample.jl")

end
