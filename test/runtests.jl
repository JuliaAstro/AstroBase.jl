using AstroBase
using ERFA
using Base.Test

# write your own tests here
@test AstroBase.sp00(2.4578265e6, 0.30434616919175345) == ERFA.sp00(2.4578265e6, 0.30434616919175345)
