export CelestialBody, Planet, NaturalSatellite, MinorBody, Barycenter
export Sun, SSB

export naif_id, μ, mu, j2, mean_radius, equatorial_radius, polar_radius
export maximum_elevation, maximum_depression, deviation

abstract type CelestialBody end
abstract type Barycenter <: CelestialBody end
abstract type Planet <: CelestialBody end
abstract type NaturalSatellite <: CelestialBody end
abstract type MinorBody <: CelestialBody end

Base.show{T<:CelestialBody}(io::IO, ::Type{T}) = print(io, Base.datatype_name(T))

const CONSTANTS = (
    :parent,
    :naif_id,
    :μ,
    :j2,
    :mean_radius,
    :equatorial_radius,
    :polar_radius,
    :deviation,
    :maximum_elevation,
    :maximum_depression,
    :ra₀,
    :ra₁,
    :ra₂,
    :dec₀,
    :dec₁,
    :dec₂,
    :w₀,
    :w₁,
    :w₂,
    :a,
    :d,
    :w,
    :θ₀,
    :θ₁,
)

for c in CONSTANTS
    @eval begin
        ($c){T<:CelestialBody}(::Type{T}) = error(
            "No constant '", $(string(c)), "' available for body '", T, "'.")
    end
end

struct SolarSystemBarycenter <: Barycenter end
const SSB = SolarSystemBarycenter

naif_id(::Type{SSB}) = 0
parent(::Type{SSB}) = SSB

struct Sun <: CelestialBody end

μ(::Type{Sun}) = 1.32712440041e11km^3/s^2
mean_radius(::Type{Sun}) = 696000.0km
naif_id(::Type{Sun}) = 10
parent(::Type{Sun}) = SSB

const mu = μ
