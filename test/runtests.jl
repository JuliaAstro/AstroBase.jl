using AstroBase
using ERFA
using Base.Test


# write your own tests here
@test AstroBase.anp(-10) == ERFA.anp(-10)
@test AstroBase.anp(10) == ERFA.anp(10)
@test AstroBase.a00(2.4578265e6, 0.30434616919175345) ≈ ERFA.era00(2.4578265e6, 0.30434616919175345)
