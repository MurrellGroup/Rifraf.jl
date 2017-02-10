using Base.Test

@testset "Quiver2" begin
    include("./test_bandedarrays.jl")
    include("./test_sample.jl")
    include("./test_proposals.jl")
    include("./test_model.jl")
    include("./test_poa.jl")
end
