using AstroBase
using Test

import SPICE

@testset "Utils" begin
    @testset "sec2rad" begin
        # Radians <-> arcseconds
        @test sec2rad(3600) == deg2rad(1)
        @test rad2sec(1) == rad2deg(1) * 3600
    end
    @testset "normalize_angle" begin
        @testset for a = -15.0:0.1:15.0, b = -15.0:0.2:15.0
            c = normalize_angle(a, b)
            @test (b - π) <= c
            @test c <= (b + π)
            two_k = round((a - c) / π)
            @test c ≈ a - two_k * π atol=1.0e-14
        end
    end
    @testset "frame" begin
        y = zeros(3)
        @test all(AstroBase.Util.frame(y) .== SPICE.frame(y))
        x = randn(3)
        @testset for i in 1:3
            exp = SPICE.frame(x)
            act = AstroBase.Util.frame(x)
            @test act[i] ≈ exp[i]
        end
    end
    @testset "plane_section" begin
        xy = plane_section([1.0, 2.0, 3.0], zeros(3), [0.0, 0.0, 1.0])
    end
end
