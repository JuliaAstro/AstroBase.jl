using AstroBase
using Test

import SPICE

@testset "Utils" begin
    @testset "Cubic Spline" begin
        x = [-1.45971551, -1.27401241, -1.24858571, -0.98854859, -0.85262672,
             -0.36269993,  0.12312452,  0.20519359,  0.82698527,  1.18579896]
        x_rev = reverse(x)
        y = [-1.12016875,  0.53793553, -0.32205336, -0.73225522,  0.53240318,
             -0.42591753,  0.96449187,  0.2450982 , -0.68313154,  0.52273895];
        spl = CubicSpline(x, y)
        spl_rev = CubicSpline(x_rev, y)
        @test spl(1.0) ≈ -0.07321713407025687
        @test spl_rev(1.0) ≈ -0.0765278000823969
        @test spl(-Inf) == y[1]
        @test spl(Inf) == y[end]
        @test_throws ArgumentError CubicSpline(x, y[1:end-1])
        @test_throws ArgumentError CubicSpline(x[1:3], y[1:3])
    end
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
