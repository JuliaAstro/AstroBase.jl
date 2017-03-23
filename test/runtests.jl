using AstroDynCoordinates
using AstroDynBase
using Base.Test

AstroDynBase.update()

@testset "AstroDynCoordinates" begin
    include("iau.jl")
    include("iers.jl")
end
