using AstroTime: Epoch
using ..Bodies: CelestialBody, ssb, SolarSystemBarycenter, ParentOfSSB

import Base: position

export AbstractEphemeris, state, position, velocity

abstract type AbstractEphemeris end

for (f, n) in zip((:state, :position, :velocity), (6, 3, 3))
    fmut = Symbol(f, "!")
    @eval begin
        $fmut(arr, ep::Epoch, to::CelestialBody) = $fmut(arr, ep, ssb, to)
        $f(ep::Epoch, to::CelestialBody) = $fmut(zeros($n), ep, ssb, to)

        function $fmut(arr, eph::E, ep::Epoch, to::CelestialBody,
            from::CelestialBody) where E<:AbstractEphemeris
            error($fmut, " not implemented for ephemeris of type '$E'.")
        end

        function $f(eph::E, ep::Epoch, to::CelestialBody,
            from::CelestialBody) where E<:AbstractEphemeris
            error($f, " not implemented for ephemeris of type '$E'.")
        end

        $fmut(arr, ep::Epoch, to::SolarSystemBarycenter) = zeros($n)
        $fmut(arr, ep::Epoch, from::ParentOfSSB, to::SolarSystemBarycenter) = zeros($n)
        $f(ep::Epoch, to::SolarSystemBarycenter) = zeros($n)
        $f(ep::Epoch, from::ParentOfSSB, to::SolarSystemBarycenter) = zeros($n)

        export $f, $fmut
    end
end
