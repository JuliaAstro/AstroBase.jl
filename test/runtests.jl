using AstroBase
using ERFA
using Base.Test

@testset "IERS" begin
    # Earth rotation angle
    @test AstroBase.earth_rotation_angle(2.4578265e6, 0.30434616919175345) ≈ ERFA.era00(2.4578265e6, 0.30434616919175345)
    
    # Polar motion
    @test AstroBase.polar_motion(30, 30, 30) == ERFA.pom00(30, 30, 30)
    @test AstroBase.polar_motion(20, 30, 50) == ERFA.pom00(20, 30, 50)
end
