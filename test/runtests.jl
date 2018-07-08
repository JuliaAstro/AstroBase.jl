using ERFA
using Base.Test


# write your own tests here

@testset "Conversions" begin
    @test AstroBase.earth_rotation_angle(2.4578265e6, 0.30434616919175345) ≈ ERFA.era00(2.4578265e6, 0.30434616919175345)
    @test AstroBase.mean_anomaly_of_moon(1)  ≈ ERFA.fal03(1)
    @test AstroBase.mean_anomaly_of_sun(1) ≈ ERFA.falp03(1)
    @test AstroBase.mean_longitude_of_moon_minus_mean_longitude_of_ascending_node(1)  ≈ ERFA.faf03(1)
    @test AstroBase.mean_elongation_of_moon_from_sun(1)  ≈ ERFA.fad03(1)
    @test AstroBase.mean_longitude_ascending_node_moon(1) ≈ mod2pi(ERFA.faom03(1))
    @test AstroBase.mean_longitude_of_mercury(1) ≈ ERFA.fame03(1)
    @test AstroBase.mean_longitude_of_venus(1) ≈ ERFA.fave03(1)
    @test AstroBase.mean_longitude_of_earth(1)  ≈ ERFA.fae03(1)
    @test AstroBase.mean_longitude_of_mars(1) ≈ ERFA.fama03(1)
    @test AstroBase.mean_longitude_of_jupiter(1) ≈ ERFA.faju03(1)
    @test AstroBase.mean_longitude_of_saturn(1) ≈ ERFA.fasa03(1)
    @test AstroBase.mean_longitude_of_uranus(1) ≈ ERFA.faur03(1)
    @test AstroBase.mean_longitude_of_neptune(1) ≈ ERFA.fane03(1)
    @test AstroBase.general_precession_in_longitude(1) ≈ ERFA.fapa03(1)
    let (x1, y1) = AstroBase.xy06(2.4578265e6, 0.30440190993249416)
        (x2, y2) = ERFA.xy06(2.4578265e6, 0.30440190993249416)
        @test x1 ≈ x2
        @test y1 ≈ y2
    end
end
