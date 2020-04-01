using AstroBase
using Test

import SPICE

@testset "Utils" begin
    @testset "Angles" begin
        x = randn()
        @test deg_to_rad(rad_to_deg(x)) ≈ x
        @test sec_to_rad(rad_to_sec(x)) ≈ x
        @test sec_to_deg(deg_to_sec(x)) ≈ x
        @test hms_to_days(days_to_hms(x)...) ≈ x
        @test dms_to_deg(deg_to_dms(x)...) ≈ x
        @test dms_to_rad(rad_to_dms(x)...) ≈ x
        @testset "normalize2pi" begin
            @testset for a = -15.0:0.5:15.0, b = -15.0:0.5:15.0
                c = normalize2pi(a, b)
                @test (b - π) <= c
                @test c <= (b + π)
                two_k = round((a - c) / π)
                @test c ≈ a - two_k * π atol=1.0e-14
            end
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
    @testset "spherical_to_cartesian" begin
        act = spherical_to_cartesian(3.0123, -0.999)
        exp = (-0.5366267667260523906, 0.0697711109765145365, -0.8409302618566214041)
        @testset for i in eachindex(exp, act)
            @test exp[i] ≈ act[i] atol=1e-12
        end

    end
end
