using AstroDynCoordinates
using AstroDynBase
using Base.Test

AstroDynBase.update()

@testset "AstroDynCoordinates" begin
    include("iau.jl")
    include("iers.jl")
    @testset "Common" begin
        ep = TDBEpoch(2000, 1, 1)
        rot = Rotation(IAUEarth, ITRF, ep)
        rot = Rotation(ITRF, IAUEarth, ep)
    end
end
