using AstroBase
using ERFA
using Base.Test

# write your own tests here
@test AstroBase.polar_motion(30, 30, 30) == ERFA.pom00(30, 30, 30)
@test AstroBase.polar_motion(20, 30, 50) == ERFA.pom00(20, 30, 50)
