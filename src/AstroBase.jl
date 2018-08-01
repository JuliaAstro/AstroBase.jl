module AstroBase

using Reexport

include("util.jl")
include(joinpath("bodies", "Bodies.jl"))
include("EarthAttitude.jl")

@reexport using .Bodies
@reexport using .EarthAttitude

end # module
