module Interfaces

using AstroTime: Epoch
using ..Bodies:
    CelestialBody,
    SolarSystemBarycenter,
    path_ids,
    ssb

# AbstractEphemeris interface

import Base: position

export AbstractEphemeris,
    position,
    position!,
    velocity,
    velocity!,
    position_velocity,
    position_velocity!

abstract type AbstractEphemeris end

for f in (:position, :velocity)
    fmut = Symbol(f, "!")
    @eval begin
        function $fmut(arr,
                       eph::AbstractEphemeris,
                       ep::Epoch,
                       from::T,
                       to::T;
                       kwargs...
                      ) where T<:CelestialBody
            arr .+= 0.0
        end

        function $fmut(arr,
                       eph::AbstractEphemeris,
                       ep::Epoch,
                       from::CelestialBody,
                       to::CelestialBody;
                       kwargs...
                      )
            path = path_ids(from, to)
            for (origin, target) in zip(path[1:end-1], path[2:end])
                $fmut(arr, eph, ep, origin, target; kwargs...)
            end
            arr
        end

        function $f(eph::AbstractEphemeris,
                    ep::Epoch,
                    from::CelestialBody,
                    to::CelestialBody;
                    kwargs...
                   )
            $fmut(zeros(3), eph, ep, from, to; kwargs...)
        end

        function $fmut(arr, eph::AbstractEphemeris, ep::Epoch, to::CelestialBody; kwargs...)
            $fmut(arr, eph, ep, ssb, to; kwargs...)
        end

        function $f(eph::AbstractEphemeris, ep::Epoch, to::CelestialBody; kwargs...)
            $fmut(zeros(3), eph, ep, ssb, to; kwargs...)
        end
    end
end

function position_velocity!(r, v, eph::AbstractEphemeris, ep::Epoch, from::CelestialBody, to::CelestialBody)
    position!(r, eph, ep, from, to), velocity!(v, eph, ep, from, to)
end

function position_velocity(eph::AbstractEphemeris, ep::Epoch, from::CelestialBody, to::CelestialBody)
    position!(zeros(3), eph, ep, from, to), velocity!(zeros(3), eph, ep, from, to)
end

function position_velocity!(r, v, eph::AbstractEphemeris, ep::Epoch, to::CelestialBody)
    position!(r, eph, ep, ssb, to), velocity!(v, eph, ep, ssb, to)
end

function position_velocity(eph::AbstractEphemeris, ep::Epoch, to::CelestialBody)
    position!(zeros(3), eph, ep, ssb, to), velocity!(zeros(3), eph, ep, ssb, to)
end

end
