using AstroBase
using ERFA
using Test

@testset "AstroBase" begin
    @testset "IERS" begin
        @test obliquity_of_ecliptic_06(2.4578265e6, 0.30434616919175345) ≈ ERFA.obl06(2.4578265e6, 0.30434616919175345)
        @test mean_obliquity_of_ecliptic(2.4578265e6, 0.30434616919175345) ≈ ERFA.obl80(2.4578265e6, 0.30434616919175345)
        @test earth_rotation_angle(2.4578265e6, 0.30434616919175345) ≈ ERFA.era00(2.4578265e6, 0.30434616919175345)

        let (w1, x1, y1, z1) = precession_fukushima_williams06(2.4578265e6, 0.30434616919175345)
            (w2, x2, y2, z2) = ERFA.pfw06(2.4578265e6, 0.30434616919175345)
            @test w1 ≈ w2
            @test x1 ≈ x2
            @test y1 ≈ y2
            @test z1 ≈ z2
        end

        # Celestial to intermediate frame of date
        @test celestial_to_intermediate(0.2, 0.2, 0.2) ≈ ERFA.c2ixys(0.2, 0.2, 0.2)

        # Polar motion
        @test polar_motion(30, 30, 30) ≈ ERFA.pom00(30, 30, 30)
        @test polar_motion(20, 30, 50) ≈ ERFA.pom00(20, 30, 50)

        @test tio_locator(2.4578265e6, 0.30434616919175345) ≈ ERFA.sp00(2.4578265e6, 0.30434616919175345)

        # Radians <-> arcseconds
        @test sec2rad(3600) == deg2rad(1)
        @test rad2sec(1) == rad2deg(1) * 3600

        @test mean_anomaly(luna, 1)  ≈ ERFA.fal03(1)
        @test mean_anomaly(sun, 1) ≈ ERFA.falp03(1)
        @test mean_longitude_minus_lan(luna, 1) ≈ ERFA.faf03(1)
        @test mean_elongation(luna, 1) ≈ ERFA.fad03(1)
        @test mean_longitude_ascending_node(luna, 1) ≈ ERFA.faom03(1)
        @test mean_longitude(mercury, 1) ≈ ERFA.fame03(1)
        @test mean_longitude(venus, 1) ≈ ERFA.fave03(1)
        @test mean_longitude(earth, 1) ≈ ERFA.fae03(1)
        @test mean_longitude(mars, 1) ≈ ERFA.fama03(1)
        @test mean_longitude(jupiter, 1) ≈ ERFA.faju03(1)
        @test mean_longitude(saturn, 1) ≈ ERFA.fasa03(1)
        @test mean_longitude(uranus, 1) ≈ ERFA.faur03(1)
        @test mean_longitude(neptune, 1) ≈ ERFA.fane03(1)
        @test general_precession_in_longitude(1) ≈ ERFA.fapa03(1)
        let (x1, y1) = xy06(2.4578265e6, 0.30440190993249416)
            (x2, y2) = ERFA.xy06(2.4578265e6, 0.30440190993249416)
            @test x1 ≈ x2
            @test y1 ≈ y2
        end

        let (x1, y1) = nutation_00b(2.4578265e6, 0.30440190993249416)
            (x2, y2) = ERFA.nut00b(2.4578265e6, 0.30440190993249416)
            @test x1 ≈ x2
            @test y1 ≈ y2
        end

        @test fukushima_williams_matrix(0.2,0.3,0.5,0.6) ≈ ERFA.fw2m(0.2,0.3,0.5,0.6)
        @test numat(0.7, 1.4, 1.3) ≈ ERFA.numat(0.7, 1.4, 1.3)

        let (ap, ao) = precession_rate_part_of_nutation(2.4578265e6, 0.30434616919175345)
            (ep, eo) = ERFA.pr00(2.4578265e6, 0.30434616919175345)
            @test ap ≈ ep
            @test ao ≈ eo
        end

        @test greenwich_mean_sidereal_time00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈  ERFA.gmst00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
        @test greenwich_mean_sidereal_time06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈  ERFA.gmst06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
        @test greenwich_mean_sidereal_time82(2.4578265e6, 0.30434616919175345) ≈ ERFA.gmst82(2.4578265e6, 0.30434616919175345)
        @test greenwich_mean_sidereal_time82(0.30434616919175345, 2.4578265e6) ≈ ERFA.gmst82(0.30434616919175345, 2.4578265e6)
        let (b1, p1, bp1) = bias_precession_matrix_00(2.4578265e6, 0.30434616919175345)
            (b2, p2, bp2) = ERFA.bp00(2.4578265e6, 0.30434616919175345)
            @test b1 ≈ b2
            @test p1 ≈ p2
            @test bp1 ≈ bp2
        end

        let a = rand(3,3)
            @test equation_of_origins(a, 0.3) ≈ ERFA.eors(a, 0.3)
        end

        let (x1, y1) = nutation(2.4578265e6, 0.30434616919175345),
            (x2, y2) = ERFA.nut80(2.4578265e6, 0.30434616919175345)
            @test x1 ≈ x2
            @test y1 ≈ y2
        end

        @test equation_of_equinoxes_complementary_terms(2.4578265e6, 0.30434616919175345) ≈ ERFA.eect00(2.4578265e6, 0.30434616919175345)
        @test equation_of_equinoxes_00(2.4578265e6, 0.30440190993249416, 1.5, 1.7) ≈ ERFA.ee00(2.4578265e6, 0.30440190993249416, 1.5, 1.7)

        @test s00(2.4578265e6, 0.30434616919175345, 0, 0) ≈ ERFA.s00(2.4578265e6, 0.30434616919175345, 0, 0)
        @test s06(2.4578265e6, 0.30434616919175345, 0, 0) ≈ ERFA.s06(2.4578265e6, 0.30434616919175345, 0, 0)

        let (x1, y1) = nutation_00a(2.4578265e6, 0.30440190993249416)
            (x2, y2) = ERFA.nut00a(2.4578265e6, 0.30440190993249416)
            @test x1 ≈ x2
            @test y1 ≈ y2
        end
    end

    @testset "Dependencies" begin
        @test nutation_matrix80(2.4578265e6, 0.30434616919175345) ≈ ERFA.nutm80(2.4578265e6, 0.30434616919175345)

    end

    include("bodies.jl")
end
