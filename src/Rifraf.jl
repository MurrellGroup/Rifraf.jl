module Rifraf

using Bio.Seq
using Distributions
using Levenshtein
using Parameters
using StatsBase

import Base.length
import Base.reverse


export rifraf,
       RifrafParams,
       ErrorModel,
       Scores,
       normalize,
       sample,
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
       align

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
