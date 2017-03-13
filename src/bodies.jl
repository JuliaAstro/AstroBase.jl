export CelestialBody, Planet, NaturalSatellite, MinorBody, Sun

abstract type CelestialBody end
abstract type Planet <: CelestialBody end
abstract type NaturalSatellite <: CelestialBody end
abstract type MinorBody <: CelestialBody end

Base.show{T<:CelestialBody}(io::IO, ::Type{T}) = print(io, Base.datatype_name(T))

const CONSTANTS = (
    :μ,
)

for c in CONSTANTS
    @eval begin
        ($c){T<:CelestialBody}(::Type{T}) = error(
            "No constant '", $(string(c)), "' available for body '", T, "'.")
        export $c
    end
end

struct Sun <: CelestialBody end

μ(::Type{Sun}) = 1.32712440041e11
mean_radius(::Type{Sun}) = 696000.0
id(::Type{Sun}) = 10

const mu = μ
