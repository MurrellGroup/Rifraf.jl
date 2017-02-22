using Base.Test

@testset "Rifraf" begin
    include("./test_rifrafsequences.jl")
    include("./test_bandedarrays.jl")
    include("./test_sample.jl")
    include("./test_proposals.jl")
    include("./test_rifraf.jl")
end
