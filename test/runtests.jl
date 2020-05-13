using AstroBase
import AstroTime
using ERFA
using Test

AstroTime.load_test_eop()

@testset "AstroBase" begin
    # include("n_body.jl")
    include("astrometry.jl")
    include("bodies.jl")
    include("constants.jl")
    include("coords.jl")
    include("earth_attitude.jl")
    include("ephemerides.jl")
    include("frames.jl")
    include("two_body.jl")
    include("util.jl")
end

