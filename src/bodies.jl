import Base: parent
import AstronomicalTime: SEC_PER_CENTURY, SEC_PER_DAY, in_seconds
export CelestialBody, Planet, NaturalSatellite, MinorBody, Barycenter
export Sun, SSB, SolarSystemBarycenter

export naif_id, μ, mu, j2, mean_radius, equatorial_radius, polar_radius
export maximum_elevation, maximum_depression, deviation, parent
export right_ascension, right_ascension_rate, declination, declination_rate
export rotation_angle, rotation_rate

abstract type CelestialBody end
abstract type Barycenter <: CelestialBody end
abstract type Planet <: CelestialBody end
abstract type NaturalSatellite <: CelestialBody end
abstract type MinorBody <: CelestialBody end

Base.show{T<:CelestialBody}(io::IO, ::Type{T}) = print(io, Base.datatype_name(T))

ra₀(::Type{<:CelestialBody}) = 0.0
ra₁(::Type{<:CelestialBody}) = 0.0
ra₂(::Type{<:CelestialBody}) = 0.0
dec₀(::Type{<:CelestialBody}) = 0.0
dec₁(::Type{<:CelestialBody}) = 0.0
dec₂(::Type{<:CelestialBody}) = 0.0
w₀(::Type{<:CelestialBody}) = 0.0
w₁(::Type{<:CelestialBody}) = 0.0
w₂(::Type{<:CelestialBody}) = 0.0
a(::Type{<:CelestialBody}) = [0.0]
d(::Type{<:CelestialBody}) = [0.0]
w(::Type{<:CelestialBody}) = [0.0]
θ₀(::Type{<:CelestialBody}) = [0.0]
θ₁(::Type{<:CelestialBody}) = [0.0]

struct SolarSystemBarycenter <: Barycenter end
const SSB = SolarSystemBarycenter

naif_id(::Type{SSB}) = 0
parent(::Type{SSB}) = SSB

struct Sun <: CelestialBody end

μ(::Type{Sun}) = 1.32712440041e11km^3/s^2
mean_radius(::Type{Sun}) = 696000.0km
naif_id(::Type{Sun}) = 10
parent(::Type{Sun}) = SSB
ra₀(::Type{Sun}) = deg2rad(286.13)
dec₀(::Type{Sun}) = deg2rad(63.87)
w₀(::Type{Sun}) = deg2rad(84.176)
w₁(::Type{Sun}) = deg2rad(84.176)

θ(t::Float64, b::Type{<:CelestialBody}) = θ₀(b) .+ θ₁(b) .* t/SEC_PER_CENTURY

function right_ascension(b::Type{C}, ep) where C<:CelestialBody
    t = in_seconds(ep)
    mod2pi(ra₀(b) + ra₁(b) * t / SEC_PER_CENTURY +
        ra₂(b) * t^2 / SEC_PER_CENTURY^2 +
        sum(a(b) .* sin.(θ(t, b))))
end

function declination(b::Type{C}, ep) where C<:CelestialBody
    t = in_seconds(ep)
    mod2pi(dec₀(b) + dec₁(b) * t / SEC_PER_CENTURY +
        dec₂(b) * t^2 / SEC_PER_CENTURY^2 +
        sum(d(b) .* cos.(θ(t, b))))
end

function rotation_angle(b::Type{C}, ep) where C<:CelestialBody
    t = in_seconds(ep)
    mod2pi(w₀(b) + w₁(b) * t / SEC_PER_DAY +
        w₂(b) * t^2 / SEC_PER_DAY^2 +
        sum(w(b) .* sin.(θ(t, b))))
end

function right_ascension_rate(b::Type{C}, ep) where C<:CelestialBody
    t = in_seconds(ep)
    ra₁(b) / SEC_PER_CENTURY + 2 * ra₂(b) * t / SEC_PER_CENTURY^2 +
        sum(a(b) .* θ₁(b) / SEC_PER_CENTURY .* cos.(θ(t, b)))
end

function declination_rate(b::Type{C}, ep) where C<:CelestialBody
    t = in_seconds(ep)
    dec₁(b) / SEC_PER_CENTURY + 2 * dec₂(b) * t / SEC_PER_CENTURY^2 +
        sum(d(b) .* θ₁(b) / SEC_PER_CENTURY .* sin.(θ(t, b)))
end

function rotation_rate(b::Type{C}, ep) where C<:CelestialBody
    t = in_seconds(ep)
    w₁(b) / SEC_PER_DAY + 2 * w₂(b) * t / SEC_PER_DAY^2 +
        sum(w(b) .* θ₁(b) / SEC_PER_CENTURY .* cos.(θ(t, b)))
end

const mu = μ
