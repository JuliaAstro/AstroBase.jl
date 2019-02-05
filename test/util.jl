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
end
