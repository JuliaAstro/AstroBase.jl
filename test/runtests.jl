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
        @test origin(rot) == IAUEarth
        @test target(rot) == ITRF
        rot = Rotation(ITRF, IAUEarth, ep)
        @test origin(rot) == ITRF
        @test target(rot) == IAUEarth
    end
    include("states.jl")
end
