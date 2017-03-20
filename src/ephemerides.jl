import Base: position

export Ephemeris

abstract type Ephemeris end

EPHEMERIS = Ref{Ephemeris}()

for f in (:state, :position, :velocity)
    @eval begin
        function $f(ep::Epoch, to::Type{T1}, from::Type{T2}=SSB) where {T1<:CelestialBody, T2<:CelestialBody}
            if !isassigned(EPHEMERIS)
                error("No ephemeris has been loaded. Use 'load_ephemeris'.")
            end
            state(EPHEMERIS[], ep, to, from)
        end

        $f(eph::T, ep, to, from) where T<:Ephemeris = error("'$f' not implemented for ephemeris '$T'.")
        export $f
    end
end

load_ephemeris(eph::Ephemeris) = EPHEMERIS[] = eph
