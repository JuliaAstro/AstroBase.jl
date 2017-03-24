using AstroDynBase
using Base.Test

struct DefunctEphemeris <: Ephemeris end

@testset "AstroDynBase" begin
    include("elements.jl")
    include("rotations.jl")
    include("rotation_matrices.jl")
    @testset "States" begin
    end
    @testset "Ephemerides" begin
        load_ephemeris!(DefunctEphemeris())
        ep = TTEpoch(2000, 1, 1)
        @test_throws ErrorException state(ep, Earth)
        @test_throws ErrorException position(ep, Earth)
        @test_throws ErrorException velocity(ep, Earth)
    end
end
