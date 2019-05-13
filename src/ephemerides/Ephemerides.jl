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
    ..state,
    ..state!

export AbstractEphemeris,
    position,
    position!,
    velocity,
    velocity!,
    state,
    state!

abstract type AbstractEphemeris end

function position end
function position! end
function velocity end
function velocity! end
function state end
function state end

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
# function state!(pos, vel, eph::AbstractEphemeris, ep::Epoch, from::NAIFId, to::NAIFId; kwargs...)
#     position!(pos, eph, ep, from, to; kwargs...), velocity!(vel, eph, ep, from, to; kwargs...)
# end
#
# function state!(pos,
#                             vel,
#                             eph::AbstractEphemeris,
#                             ep::Epoch,
#                             from::CelestialBody,
#                             to::CelestialBody;
#                             kwargs...)
#     path = path_ids(from, to)
#     for (origin, target) in zip(path[1:end-1], path[2:end])
#         state!(pos, vel, eph, ep, origin, target; kwargs...)
#     end
#     pos, vel
# end
#
# state!(pos,
#                    vel,
#                    ::AbstractEphemeris,
#                    ::Epoch,
#                    ::T, ::T; kwargs...) where {T<:CelestialBody} = pos, vel
#
# function state(eph::AbstractEphemeris, ep::Epoch, from::CelestialBody, to::CelestialBody; kwargs...)
#     state!(zeros(3), zeros(3), eph, ep, from, to; kwargs...)
# end
#
# function state!(pos, vel, eph::AbstractEphemeris, ep::Epoch, to::CelestialBody; kwargs...)
#     state!(pos, vel, eph, ep, ssb, to; kwargs...)
# end
#
# function state(eph::AbstractEphemeris, ep::Epoch, to::CelestialBody; kwargs...)
#     state!(zeros(3), zeros(3), eph, ep, ssb, to; kwargs...)
# end

include("vsop87.jl")

end
