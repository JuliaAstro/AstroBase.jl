using OptionalData
import Base: position

export Ephemeris, load_ephemeris!

abstract type Ephemeris end

@OptionalData EPHEMERIS Ephemeris

for f in (:state, :position, :velocity)
    @eval begin
        function $f(ep::Epoch, from::Type{<:CelestialBody},
            to::Type{<:CelestialBody})
            $f(get(EPHEMERIS), TDBEpoch(ep), to, from)
        end
        $f(ep::Epoch, to::Type{<:CelestialBody}) = $f(ep, SSB, to)

        function $f(eph::E, ep::Epoch, to::Type{<:CelestialBody},
            from::Type{<:CelestialBody}) where E<:Ephemeris
            error($f, " not implemented for ephemeris $E.")
        end

        export $f
    end
end

function load_ephemeris!(eph::E) where E<:Ephemeris
    push!(EPHEMERIS, eph)
end

function load_ephemeris!(::Type{E}, file) where E<:Ephemeris
    push!(EPHEMERIS, E, file)
end
