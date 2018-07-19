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
    @test AstroBase.greenwich_mean_sidereal_time00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈  ERFA.gmst00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
    @test AstroBase.greenwich_mean_sidereal_time06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851) ≈  ERFA.gmst06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
end
