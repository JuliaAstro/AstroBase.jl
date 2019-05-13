module Ephemerides

using AstroTime: Epoch
using ..Bodies:
    CelestialBody,
    NAIFId,
    SolarSystemBarycenter,
    naifid,
    path_ids,
    ssb

#####################
# AbstractEphemeris #
#####################

import ..position,
    ..position!,
    ..velocity,
    ..velocity!,
    ..position_velocity,
    ..position_velocity!

export AbstractEphemeris,
    position,
    position!,
    velocity,
    velocity!,
    position_velocity,
    position_velocity!

abstract type AbstractEphemeris end

function position end
function position! end
function velocity end
function velocity! end
function position_velocity end
function position_velocity end

# for f in (:position, :velocity)
#     fmut = Symbol(f, "!")
#     @eval begin
#         function $fmut(arr,
#                        eph::AbstractEphemeris,
#                        ep::Epoch,
#                        from::T,
#                        to::T;
#                        kwargs...
#                       ) where T<:CelestialBody
#             arr .+= 0.0
#         end
#
#         function $fmut(arr,
#                        eph::AbstractEphemeris,
#                        ep::Epoch,
#                        from::CelestialBody,
#                        to::CelestialBody;
#                        kwargs...
#                       )
#             path = path_ids(from, to)
#             for (origin, target) in zip(path[1:end-1], path[2:end])
#                 $fmut(arr, eph, ep, origin, target; kwargs...)
#             end
#             arr
#         end
#
#         function $f(eph::AbstractEphemeris,
#                     ep::Epoch,
#                     from::CelestialBody,
#                     to::CelestialBody;
#                     kwargs...
#                    )
#             $fmut(zeros(3), eph, ep, from, to; kwargs...)
#         end
#
#         function $fmut(arr, eph::AbstractEphemeris, ep::Epoch, to::CelestialBody; kwargs...)
#             $fmut(arr, eph, ep, ssb, to; kwargs...)
#         end
#
#         function $f(eph::AbstractEphemeris, ep::Epoch, to::CelestialBody; kwargs...)
#             $fmut(zeros(3), eph, ep, ssb, to; kwargs...)
#         end
#     end
# end
#
# function position_velocity!(pos, vel, eph::AbstractEphemeris, ep::Epoch, from::NAIFId, to::NAIFId; kwargs...)
#     position!(pos, eph, ep, from, to; kwargs...), velocity!(vel, eph, ep, from, to; kwargs...)
# end
#
# function position_velocity!(pos,
#                             vel,
#                             eph::AbstractEphemeris,
#                             ep::Epoch,
#                             from::CelestialBody,
#                             to::CelestialBody;
#                             kwargs...)
#     path = path_ids(from, to)
#     for (origin, target) in zip(path[1:end-1], path[2:end])
#         position_velocity!(pos, vel, eph, ep, origin, target; kwargs...)
#     end
#     pos, vel
# end
#
# position_velocity!(pos,
#                    vel,
#                    ::AbstractEphemeris,
#                    ::Epoch,
#                    ::T, ::T; kwargs...) where {T<:CelestialBody} = pos, vel
#
# function position_velocity(eph::AbstractEphemeris, ep::Epoch, from::CelestialBody, to::CelestialBody; kwargs...)
#     position_velocity!(zeros(3), zeros(3), eph, ep, from, to; kwargs...)
# end
#
# function position_velocity!(pos, vel, eph::AbstractEphemeris, ep::Epoch, to::CelestialBody; kwargs...)
#     position_velocity!(pos, vel, eph, ep, ssb, to; kwargs...)
# end
#
# function position_velocity(eph::AbstractEphemeris, ep::Epoch, to::CelestialBody; kwargs...)
#     position_velocity!(zeros(3), zeros(3), eph, ep, ssb, to; kwargs...)
# end

include("vsop87.jl")

end
