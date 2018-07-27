using AstroBase
using ERFA
using Test

@testset "AstroBase" begin
    @testset "IERS" begin
        @test AstroBase.earth_rotation_angle(2.4578265e6, 0.30434616919175345) ≈ ERFA.era00(2.4578265e6, 0.30434616919175345)

        # Celestial to intermediate frame of date
        @test AstroBase.celestial_to_intermediate(0.2, 0.2, 0.2) ≈ ERFA.c2ixys(0.2, 0.2, 0.2)

        # Polar motion
        @test AstroBase.polar_motion(30, 30, 30) ≈ ERFA.pom00(30, 30, 30)
        @test AstroBase.polar_motion(20, 30, 50) ≈ ERFA.pom00(20, 30, 50)

        @test AstroBase.tio_locator(2.4578265e6, 0.30434616919175345) ≈ ERFA.sp00(2.4578265e6, 0.30434616919175345)

        # Radians <-> arcseconds
        @test AstroBase.sec2rad(3600) == deg2rad(1)
        @test AstroBase.rad2sec(1) == rad2deg(1) * 3600

        @test AstroBase.mean_anomaly(luna, 1)  ≈ ERFA.fal03(1)
        @test AstroBase.mean_anomaly(sun, 1) ≈ ERFA.falp03(1)
        @test AstroBase.mean_longitude_minus_lan(luna, 1) ≈ ERFA.faf03(1)
        @test AstroBase.mean_elongation(luna, 1) ≈ ERFA.fad03(1)
        @test AstroBase.mean_longitude_ascending_node(luna, 1) ≈ ERFA.faom03(1)
        @test AstroBase.mean_longitude(mercury, 1) ≈ ERFA.fame03(1)
        @test AstroBase.mean_longitude(venus, 1) ≈ ERFA.fave03(1)
        @test AstroBase.mean_longitude(earth, 1) ≈ ERFA.fae03(1)
        @test AstroBase.mean_longitude(mars, 1) ≈ ERFA.fama03(1)
        @test AstroBase.mean_longitude(jupiter, 1) ≈ ERFA.faju03(1)
        @test AstroBase.mean_longitude(saturn, 1) ≈ ERFA.fasa03(1)
        @test AstroBase.mean_longitude(uranus, 1) ≈ ERFA.faur03(1)
        @test AstroBase.mean_longitude(neptune, 1) ≈ ERFA.fane03(1)
        @test AstroBase.general_precession_in_longitude(1) ≈ ERFA.fapa03(1)
        let (x1, y1) = AstroBase.xy06(2.4578265e6, 0.30440190993249416)
            (x2, y2) = ERFA.xy06(2.4578265e6, 0.30440190993249416)
            @test x1 ≈ x2
            @test y1 ≈ y2
        end
        a, b = rand(3, 3), rand(3, 3)
        @test AstroBase.celestial_to_terrestrial_matrix(a, 0.2, b) ≈ ERFA.c2tcio(a, 0.2, b)
    end

    include("bodies.jl")
end
