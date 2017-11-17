import Base: parent, show, datatype_name
import AstronomicalTime: SECONDS_PER_CENTURY, SECONDS_PER_DAY, J2000, seconds

export CelestialBody, Planet, NaturalSatellite, MinorBody, Barycenter,
    Sun, SSB, SolarSystemBarycenter

export naif_id, μ, mu, j2, mean_radius, equatorial_radius, polar_radius,
    maximum_elevation, maximum_depression, deviation, parent, show,
    right_ascension, right_ascension_rate, declination, declination_rate,
    rotation_angle, rotation_rate, euler_angles, euler_derivatives

abstract type CelestialBody end
abstract type Barycenter <: CelestialBody end
abstract type Planet <: CelestialBody end
abstract type NaturalSatellite <: CelestialBody end
abstract type MinorBody <: CelestialBody end

α₀(::Type{<:CelestialBody}) = 0.0
α₁(::Type{<:CelestialBody}) = 0.0
α₂(::Type{<:CelestialBody}) = 0.0
δ₀(::Type{<:CelestialBody}) = 0.0
δ₁(::Type{<:CelestialBody}) = 0.0
δ₂(::Type{<:CelestialBody}) = 0.0
w₀(::Type{<:CelestialBody}) = 0.0
w₁(::Type{<:CelestialBody}) = 0.0
w₂(::Type{<:CelestialBody}) = 0.0
α(::Type{<:CelestialBody}) = zeros(1)
δ(::Type{<:CelestialBody}) = zeros(1)
w(::Type{<:CelestialBody}) = zeros(1)
θ₀(::Type{<:CelestialBody}) = zeros(1)
θ₁(::Type{<:CelestialBody}) = zeros(1)

struct SolarSystemBarycenter <: Barycenter end
const SSB = SolarSystemBarycenter
struct ParentOfSSB <: Barycenter end

naif_id(::Type{SSB}) = 0
body(::Type{SSB}) = Sun
parent(::Type{SSB}) = ParentOfSSB
parent(::Type{ParentOfSSB}) = ParentOfSSB

μ(::Type{B}) where B<:Barycenter = μ(body(B))

struct Sun <: CelestialBody end

μ(::Type{Sun}) = 1.32712440041e11
mean_radius(::Type{Sun}) = 696000.0
naif_id(::Type{Sun}) = 10
parent(::Type{Sun}) = SSB
α₀(::Type{Sun}) = deg2rad(286.13)
δ₀(::Type{Sun}) = deg2rad(63.87)
w₀(::Type{Sun}) = deg2rad(84.176)
w₁(::Type{Sun}) = deg2rad(84.176)

θ(t, b::Type{<:CelestialBody}) = θ₀(b) .+ θ₁(b) .* t ./ SECONDS_PER_CENTURY

function right_ascension(b::Type{C}, ep) where C<:CelestialBody
    t = get(seconds(ep, J2000))
    α₀(b) + α₁(b) * t / SECONDS_PER_CENTURY +
        α₂(b) * t^2 / SECONDS_PER_CENTURY^2 +
        sum(α(b) .* sin.(θ(t, b)))
end

function declination(b::Type{C}, ep) where C<:CelestialBody
    t = get(seconds(ep, J2000))
    δ₀(b) + δ₁(b) * t / SECONDS_PER_CENTURY +
        δ₂(b) * t^2 / SECONDS_PER_CENTURY^2 +
        sum(δ(b) .* cos.(θ(t, b)))
end

function rotation_angle(b::Type{C}, ep) where C<:CelestialBody
    t = get(seconds(ep, J2000))
    w₀(b) + w₁(b) * t / SECONDS_PER_DAY +
        w₂(b) * t^2 / SECONDS_PER_DAY^2 +
        sum(w(b) .* sin.(θ(t, b)))
end

function right_ascension_rate(b::Type{C}, ep) where C<:CelestialBody
    t = get(seconds(ep, J2000))
    α₁(b) / SECONDS_PER_CENTURY + 2 * α₂(b) * t / SECONDS_PER_CENTURY^2 +
        sum(α(b) .* θ₁(b) ./ SECONDS_PER_CENTURY .* cos.(θ(t, b)))
end

function declination_rate(b::Type{C}, ep) where C<:CelestialBody
    t = get(seconds(ep, J2000))
    δ₁(b) / SECONDS_PER_CENTURY + 2 * δ₂(b) * t / SECONDS_PER_CENTURY^2 -
        sum(δ(b) .* θ₁(b) ./ SECONDS_PER_CENTURY .* sin.(θ(t, b)))
end

function rotation_rate(b::Type{C}, ep) where C<:CelestialBody
    t = get(seconds(ep, J2000))
    w₁(b) / SECONDS_PER_DAY + 2 * w₂(b) * t / SECONDS_PER_DAY^2 +
        sum(w(b) .* θ₁(b) ./ SECONDS_PER_CENTURY .* cos.(θ(t, b)))
end

function euler_angles(b::Type{C}, ep) where C<:CelestialBody
    right_ascension(b, ep) + π/2, π/2 - declination(b, ep), mod2pi(rotation_angle(b, ep))
end

function euler_derivatives(b::Type{C}, ep) where C<:CelestialBody
    right_ascension_rate(b, ep), -declination_rate(b, ep), rotation_rate(b, ep)
end

const mu = μ
