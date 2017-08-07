using Base.Test

@testset "Rifraf" begin
    include("./test_rifrafsequences.jl")
    include("./test_bandedarrays.jl")
    include("./test_sample.jl")
    include("./test_proposals.jl")
    include("./test_align.jl")
    include("./test_model.jl")
    include("./test_correct_shifts.jl")
    include("./test_estimate_probs.jl")
end
