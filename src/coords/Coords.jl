module Coords

import Base: ==

using AstroTime: Epoch, TimeScale

import ..position,
    ..state,
    ..velocity

using ..Bodies:
    CelestialBody,
    earth,
    grav_param

using ..Ephemerides: AbstractEphemeris

using ..Frames:
    AbstractFrame,
    icrf

import ..TwoBody:
    cartesian,
    keplerian,
    period

using StaticArrays: SVector

include("states.jl")
include("conversions.jl")
include("bodyfixed.jl")
include("events.jl")
include("time_series.jl")
include("trajectories.jl")

end
