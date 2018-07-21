using AstroBase
using ERFA
using Base.Test

@testset "IERS" begin
    # Earth rotation angle
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

    let (ap, ao) = AstroBase.precession_rate_part_of_nutation(2.4578265e6, 0.30434616919175345)
        (ep, eo) = ERFA.pr00(2.4578265e6, 0.30434616919175345)
        @test ap ≈ ep
        @test ao ≈ eo
    end

    let (b1, p1, bp1) = AstroBase.bias_precession_matrix_00(2.4578265e6, 0.30434616919175345)
        (b2, p2, bp2) = ERFA.bp00(2.4578265e6, 0.30434616919175345)
        @test b1 ≈ b2
        @test p1 ≈ p2
        @test bp1 ≈ bp2
    end
end
