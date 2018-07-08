using AstroBase
using ERFA
using Base.Test


# write your own tests here
@testset "Conversions" begin
    @test AstroBase.earth_rotation_angle(2.4578265e6, 0.30434616919175345) ≈ ERFA.era00(2.4578265e6, 0.30434616919175345)
end
