using Base.Test

@testset "Quiver2" begin
    include("./test_bandedarrays.jl")
    include("./test_sample.jl")
    include("./test_proposals.jl")
    include("./test_align.jl")
    include("./test_rifraf.jl")
end
