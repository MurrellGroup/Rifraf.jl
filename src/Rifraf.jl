module Rifraf

using Bio.Seq
using Distributions
using Levenshtein

import Base.length
import Base.reverse


export rifraf,
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
       BandedArray,
       sparsecol,
       data_row,
       row_range,
       data_row_range,
       inband,
       flip


# keep this as a submodule. might eventually break it out into its own
# project
include("bandedarrays.jl")
using Rifraf.BandedArrays

include("types.jl")
include("phred.jl")
include("util.jl")
include("fastxio.jl")
include("proposals.jl")
include("rifrafsequences.jl")
include("errormodel.jl")
include("align.jl")
include("model.jl")
include("sample.jl")

end

