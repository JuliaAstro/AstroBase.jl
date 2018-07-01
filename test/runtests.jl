using AstroBase
using ERFA
using Base.Test

# write your own tests here
@test AstroBase.fal03(1)  == ERFA.fal03(1)
@test AstroBase.falp03(1) == ERFA.falp03(1)
@test AstroBase.faf03(1)  == ERFA.faf03(1)
@test AstroBase.fad03(1)  == ERFA.fad03(1)
#@test AstroBase.faom03(1) == ERFA.faom03(1)
@test AstroBase.fame03(1) == ERFA.fame03(1)
@test AstroBase.fave03(1) == ERFA.fave03(1)
@test AstroBase.fae03(1)  == ERFA.fae03(1)
@test AstroBase.fama03(1) == ERFA.fama03(1)
@test AstroBase.faju03(1) == ERFA.faju03(1)
@test AstroBase.fasa03(1) == ERFA.fasa03(1)
@test AstroBase.faur03(1) == ERFA.faur03(1)
@test AstroBase.fane03(1) == ERFA.fane03(1)
@test AstroBase.fapa03(1) == ERFA.fapa03(1)
