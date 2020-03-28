module EarthAttitude

using AstroTime: days
using ..Util: sec2rad
using ReferenceFrameRotations: angle_to_dcm, angleaxis_to_dcm, compose_rotation

include("iau_models.jl")
include("fundamental.jl")
include("obliquity.jl")
include("nutation.jl")
include("precession.jl")
include("precession_nutation.jl")
include("icrs.jl")
include("rotation.jl")

const ARCSECONDS_IN_CIRCLE = 1296000.0


end
