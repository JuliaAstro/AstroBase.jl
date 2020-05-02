using AstroBase
import AstroTime
using ERFA
using Test

AstroTime.load_test_eop()

@testset "AstroBase" begin
    include("constants.jl")
    include("util.jl")
    include("earth_attitude.jl")
    include("bodies.jl")
    include("ephemerides.jl")
    include("two_body.jl")
    # include("n_body.jl")
    include("frames.jl")
    include("coords.jl")
end

