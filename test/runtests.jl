using AstroBase
using ERFA
using Test

@testset "AstroBase" begin
    include("earth_attitude.jl")
    include("bodies.jl")
    include("interfaces.jl")
end

