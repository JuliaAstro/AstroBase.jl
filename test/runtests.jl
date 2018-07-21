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

    @test AstroBase.obliquity_of_ecliptic_06(2.4578265e6, 0.30434616919175345) ≈ ERFA.obl06(2.4578265e6, 0.30434616919175345)

    let (w1, x1, y1, z1) = AstroBase.precession_fukushima_williams06(2.4578265e6, 0.30434616919175345)
        (w2, x2, y2, z2) = ERFA.pfw06(2.4578265e6, 0.30434616919175345)
        @test w1 ≈ w2
        @test x1 ≈ x2
        @test y1 ≈ y2
        @test z1 ≈ z2
    end
end
