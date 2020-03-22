using AstroBase
import EarthOrientation
using ERFA
using Test

EarthOrientation.update()

@testset "AstroBase" begin
    include("util.jl")
    include("earth_attitude.jl")
    include("bodies.jl")
    include("ephemerides.jl")
    include("two_body.jl")
    include("n_body.jl")
    include("frames.jl")
    include("coords.jl")
end

