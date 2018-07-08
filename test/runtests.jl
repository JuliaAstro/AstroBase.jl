using AstroBase
using ERFA
using Base.Test

# write your own tests here
@test AstroBase.tio_locator(2.4578265e6, 0.30434616919175345) == ERFA.sp00(2.4578265e6, 0.30434616919175345)
@test AstroBase.sec2rad(3600) == deg2rad(1)
@test AstroBase.rad2sec(1) == rad2deg(1) * 3600
