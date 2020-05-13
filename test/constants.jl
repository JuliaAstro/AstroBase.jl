using AstroBase.Constants
using AstroBase.Time: TCB, TCG, TDB, TT

@testset "Constants" begin
    @test astronomical_unit()[1] ≈ 1.49597870700e8
    @test astronomical_unit()[2] ≈ 3e-3
    @test drift_rate(TT, TCG) ≈ 6.969290134e-10
    @test drift_rate(TDB, TCB) ≈ 1.550519768e-8
    @test drift_rate(TCG, TCB)[1] ≈ 1.48082686741e-8
    @test drift_rate(TCG, TCB)[2] ≈ 2e-17
    @test earth_rotation_angle_j2000() ≈ 0.7790572732640
    @test earth_rotation_rate() ≈ 1.00273781191135448
    @test equatorial_radius_earth()[1] ≈ 6.3781366e3
    @test equatorial_radius_earth()[2] ≈ 1e-4
    @test gaussian_gravitational_const() ≈ 1.720209895e-2
    @test geocentric_gravitational_const()[1] ≈ 3.986004356e5
    @test geocentric_gravitational_const()[2] ≈ 8e-4
    @test geocentric_gravitational_const(TT)[1] ≈ 3.986004415e5
    @test geocentric_gravitational_const(TT)[2] ≈ 8e-4
    @test geocentric_gravitational_const(TCB)[1] ≈ 3.986004418e5
    @test geocentric_gravitational_const(TCB)[2] ≈ 8e-4
    @test geoid_potential()[1] ≈ 6.26368560e1
    @test geoid_potential()[2] ≈ 5e-7
    @test gravitational_const()[1] ≈ 6.67428e-20
    @test gravitational_const()[2] ≈ 6.7e-24
    @test heliocentric_gravitational_const()[1] ≈ 1.32712440041e11
    @test heliocentric_gravitational_const()[2] ≈ 1e1
    @test heliocentric_gravitational_const(TCB)[1] ≈ 1.32712442099e11
    @test heliocentric_gravitational_const(TCB)[2] ≈ 1e1
    @test j2_factor_earth()[1] ≈ 1.0826359e-3
    @test j2_factor_earth()[2] ≈ 1e-10
    @test j2_rate_earth()[1] ≈ -3e-9
    @test j2_rate_earth()[2] ≈ 6e-10
    @test mean_angular_velocity_earth() ≈ 7.292115e-5
    @test obliquity_j2000()[1] ≈ 8.4381406e4
    @test obliquity_j2000()[2] ≈ 1e-3
    @test speed_of_light() ≈ 2.99792458e5
    @test tdb0() ≈ -6.55e-5
    @test schwarzschild_radius_sun() ≈ 1.97412574336e-8

    @test mass_ratio_moon()[1] ≈ 1.23000371e-2
    @test mass_ratio_moon()[2] ≈ 4e-10
    @test mass_ratio_mercury()[1] ≈ 6.0236e6
    @test mass_ratio_mercury()[2] ≈ 3e2
    @test mass_ratio_venus()[1] ≈ 4.08523719e5
    @test mass_ratio_venus()[2] ≈ 8e-3
    @test mass_ratio_mars()[1] ≈ 3.09870359e6
    @test mass_ratio_mars()[2] ≈ 2e-2
    @test mass_ratio_jupiter()[1] ≈ 1.047348644e3
    @test mass_ratio_jupiter()[2] ≈ 1.7e-5
    @test mass_ratio_saturn()[1] ≈ 3.4979018e3
    @test mass_ratio_saturn()[2] ≈ 1e-4
    @test mass_ratio_uranus()[1] ≈ 2.290298e4
    @test mass_ratio_uranus()[2] ≈ 3e-2
    @test mass_ratio_neptune()[1] ≈ 1.941226e4
    @test mass_ratio_neptune()[2] ≈ 3e-2
    @test mass_ratio_pluto()[1] ≈ 1.36566e8
    @test mass_ratio_pluto()[2] ≈ 2.8e4
    @test mass_ratio_eris()[1] ≈ 1.191e8
    @test mass_ratio_eris()[2] ≈ 1.4e6
    @test mass_ratio_ceres()[1] ≈ 4.72e-10
    @test mass_ratio_ceres()[2] ≈ 3e-12
    @test mass_ratio_pallas()[1] ≈ 1.03e-10
    @test mass_ratio_pallas()[2] ≈ 3e-12
    @test mass_ratio_vesta()[1] ≈ 1.35e-10
    @test mass_ratio_vesta()[2] ≈ 3e-12
end

