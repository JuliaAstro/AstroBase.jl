module AstroBase

using Reexport
using Rotations

@reexport using AstroTime

struct AstroException <: Exception
    msg::String
end

Base.showerror(io::IO, err::AstroException) = print(io, err.msg)

include("util.jl")
include(joinpath("bodies", "Bodies.jl"))
include("Interfaces.jl")
include("EarthAttitude.jl")
include("TwoBody.jl")
include("ephemerides.jl")

@reexport using .Bodies
@reexport using .EarthAttitude
@reexport using .TwoBody

end # module
