module AstroBase

using Reexport

@reexport using AstroTime

struct AstroException <: Exception
    msg::String
end

Base.showerror(io::IO, err::AstroException) = print(io, err.msg)

include(joinpath("constants", "Constants.jl"))
include("stubs.jl")
include(joinpath("util", "Util.jl"))
include(joinpath("bodies", "Bodies.jl"))
include(joinpath("earth_attitude", "EarthAttitude.jl"))
include(joinpath("two_body", "TwoBody.jl"))
include(joinpath("ephemerides", "Ephemerides.jl"))
include(joinpath("n_body", "NBody.jl"))
include(joinpath("frames", "Frames.jl"))
include(joinpath("coords", "Coords.jl"))

@reexport using .Util
@reexport using .Bodies
@reexport using .Ephemerides
@reexport using .EarthAttitude
@reexport using .TwoBody
@reexport using .Frames
@reexport using .Coords

end # module
