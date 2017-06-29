using OptionalData
import Base: position

export Ephemeris, load_ephemeris!, state, position, velocity

abstract type Ephemeris end

@OptionalData EPHEMERIS Ephemeris

for (f, n) in zip((:state, :position, :velocity), (6, 3, 3))
    fmut = Symbol(f, "!")
    @eval begin
        function $fmut(arr, ep::Epoch, from::Type{<:CelestialBody},
            to::Type{<:CelestialBody})
            $fmut(arr, get(EPHEMERIS), TDBEpoch(ep), to, from)
        end

        function $f(ep::Epoch, from::Type{<:CelestialBody},
            to::Type{<:CelestialBody})
            arr = zeros($n)
            $fmut(arr, ep, from, to)
        end

        $fmut(arr, ep::Epoch, to::Type{<:CelestialBody}) = $fmut(arr, ep, SSB, to)
        $f(ep::Epoch, to::Type{<:CelestialBody}) = $fmut(zeros($n), ep, SSB, to)

        function $fmut(arr, eph::E, ep::Epoch, to::Type{<:CelestialBody},
            from::Type{<:CelestialBody}) where E<:Ephemeris
            error($fmut, " not implemented for ephemeris $E.")
        end

        function $f(eph::E, ep::Epoch, to::Type{<:CelestialBody},
            from::Type{<:CelestialBody}) where E<:Ephemeris
            error($f, " not implemented for ephemeris $E.")
        end

        $fmut(arr, ep::Epoch, to::Type{SSB}) = zeros($n)
        $fmut(arr, ep::Epoch, from::Type{ParentOfSSB}, to::Type{SSB}) = zeros($n)
        $f(ep::Epoch, to::Type{SSB}) = zeros($n)
        $f(ep::Epoch, from::Type{ParentOfSSB}, to::Type{SSB}) = zeros($n)

        export $f, $fmut
    end
end

function load_ephemeris!(eph::E) where E<:Ephemeris
    push!(EPHEMERIS, eph)
end

function load_ephemeris!(::Type{E}, file) where E<:Ephemeris
    push!(EPHEMERIS, E, file)
end
