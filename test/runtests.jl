using AstroBase
using Base.Test

# write your own tests here
@test AstroBase.c2ixys(0.2, 0.2, 0.2) == ERFA.c2ixys(0.2, 0.2, 0.2)
