using AstroDynBase
using Base.Test


@testset "AstroDynBase" begin
    @testset "Rotations" begin
        include("rotations.jl")
    end
    @testset "States" begin
    end
end
