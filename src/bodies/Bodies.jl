module Bodies

import Base: parent

using AstroTime: SECONDS_PER_DAY, SECONDS_PER_CENTURY, value, j2000, seconds
using ItemGraphs: ItemGraph, SimpleGraph, add_edge!, add_vertex!, items

export CelestialBody,
    NAIFId,
    SolarSystemBarycenter,
    Sun,
    declination,
    declination_rate,
    euler_angles,
    euler_rates,
    from_naifid,
    grav_param,
    naifid,
    parent,
    right_ascension,
    right_ascension_rate,
    rotation_angle,
    rotation_rate,
    ssb,
    sun

const NAIFId = Int
const BODIES = ItemGraph{NAIFId, NAIFId}(SimpleGraph())

abstract type CelestialBody end

function grav_param end
function mean_radius end
function polar_radius end
function equatorial_radius end

Base.show(io::IO, body::CelestialBody) = print(io, string(nameof(typeof(body))))

from_naifid(id::NAIFId) = from_naifid(Val(id))

path_ids(from::CelestialBody, to::CelestialBody) = items(BODIES, naifid(from), naifid(to))

abstract type Barycenter <: CelestialBody end

struct SolarSystemBarycenter <: Barycenter end
const ssb = SolarSystemBarycenter()
Base.show(io::IO, ::SolarSystemBarycenter) = print(io, "Solar System Barycenter")
naifid(::SolarSystemBarycenter) = 0
from_naifid(::Val{0}) = ssb

struct Sun <: CelestialBody end
const sun = Sun()
parent(::Sun) = ssb
naifid(::Sun) = 10
from_naifid(::Val{10}) = sun

alpha0(::CelestialBody) = 0.0
alpha1(::CelestialBody) = 0.0
alpha2(::CelestialBody) = 0.0
delta0(::CelestialBody) = 0.0
delta1(::CelestialBody) = 0.0
delta2(::CelestialBody) = 0.0
omega0(::CelestialBody) = 0.0
omega1(::CelestialBody) = 0.0
omega2(::CelestialBody) = 0.0
alpha(::CelestialBody) = zeros(1)
delta(::CelestialBody) = zeros(1)
omega(::CelestialBody) = zeros(1)
theta0(::CelestialBody) = zeros(1)
theta1(::CelestialBody) = zeros(1)

const Δtc = SECONDS_PER_CENTURY
const Δtd = SECONDS_PER_DAY

theta(b::CelestialBody, t) = theta0(b) .+ theta1(b) .* t / Δtc

function right_ascension(b::CelestialBody, ep)
    t = value(seconds(j2000(ep)))
    α = alpha(b)
    α₀ = alpha0(b)
    α₁ = alpha1(b)
    α₂ = alpha2(b)
    θ = theta(b, t)
    α₀ + α₁ * t / Δtc + α₂ * t^2 / Δtc^2 + sum(α .* sin.(θ))
end

function declination(b::CelestialBody, ep)
    t = value(seconds(j2000(ep)))
    δ = delta(b)
    δ₀ = delta0(b)
    δ₁ = delta1(b)
    δ₂ = delta2(b)
    θ = theta(b, t)
    δ₀ + δ₁ * t / Δtc + δ₂ * t^2 / Δtc^2 + sum(δ .* cos.(θ))
end

function rotation_angle(b::CelestialBody, ep)
    t = value(seconds(j2000(ep)))
    ω = omega(b)
    ω₀ = omega0(b)
    ω₁ = omega1(b)
    ω₂ = omega2(b)
    θ = theta(b, t)
    ω₀ + ω₁ * t / Δtd + ω₂ * t^2 / Δtd^2 + sum(ω .* sin.(θ))
end

function right_ascension_rate(b::CelestialBody, ep)
    t = value(seconds(j2000(ep)))
    α = alpha(b)
    α₁ = alpha1(b)
    α₂ = alpha2(b)
    θ = theta(b, t)
    θ₁ = theta1(b)
    α₁ / Δtc + 2 * α₂ * t / Δtc^2 + sum(α .* θ₁ ./ Δtc .* cos.(θ))
end

function declination_rate(b::CelestialBody, ep)
    t = value(seconds(j2000(ep)))
    δ = delta(b)
    δ₁ = delta1(b)
    δ₂ = delta2(b)
    θ = theta(b, t)
    θ₁ = theta1(b)
    δ₁ / Δtc + 2 * δ₂ * t / Δtc^2 - sum(δ .* θ₁ ./ Δtc .* sin.(θ))
end

function rotation_rate(b::CelestialBody, ep)
    t = value(seconds(j2000(ep)))
    ω = omega(b)
    ω₁ = omega1(b)
    ω₂ = omega2(b)
    θ = theta(b, t)
    θ₁ = theta1(b)
    ω₁ / Δtd + 2 * ω₂ * t / Δtd^2 + sum(ω .* θ₁ ./ Δtc .* cos.(θ))
end

function euler_angles(b::CelestialBody, ep)
    right_ascension(b, ep) + π/2, π/2 - declination(b, ep), mod2pi(rotation_angle(b, ep))
end

function euler_rates(b::CelestialBody, ep)
    right_ascension_rate(b, ep), -declination_rate(b, ep), rotation_rate(b, ep)
end

include("planets.jl")
include("minor.jl")
include("satellites.jl")

const ALL_NAMES = Tuple([["Sun", "SolarSystemBarycenter"]
                         collect(PLANET_NAMES .* "Barycenter");
                         collect(PLANET_NAMES);
                         ["Luna", "Phobos", "Deimos"];
                         collect(JUPITER_SATELLITE_NAMES);
                         collect(SATURN_SATELLITE_NAMES);
                         collect(URANUS_SATELLITE_NAMES);
                         collect(NEPTUNE_SATELLITE_NAMES);
                         collect(MINOR_BODY_NAMES);
                         ["Pluto", "PlutoBarycenter"];
                         collect(PLUTO_SATELLITE_NAMES)])

gm = joinpath(@__DIR__, "..", "..", "gen", "gm.jl")
isfile(gm) || error("`gm.jl` has not been generated, yet. Run `gen/gen_gm.jl`")
include(gm)

pck = joinpath(@__DIR__, "..", "..", "gen", "pck.jl")
isfile(pck) || error("`pck.jl` has not been generated, yet. Run `gen/gen_pck.jl`")
include(pck)

function __init__()
    # Sun and SSB
    add_vertex!(BODIES, 0)
    add_edge!(BODIES, 0, 10)

    # Planets and barycenters
    for i in 1:length(PLANET_NAMES)
        id = 100i + 99
        add_edge!(BODIES, 0, i)
        add_edge!(BODIES, i, id)
    end

    # Earth satellite
    add_edge!(BODIES, 3, 301)

    # Mars satellites
    add_edge!(BODIES, 4, 401)
    add_edge!(BODIES, 4, 402)

    # Jupiter satellites
    for i in 1:length(JUPITER_SATELLITE_NAMES[1:end-1])
        id = 500 + i
        add_edge!(BODIES, 5, id)
    end
    add_edge!(BODIES, 5, 553)

    # Saturn satellites
    for i in 1:length(SATURN_SATELLITE_NAMES)
        id = 600 + i
        add_edge!(BODIES, 6, id)
    end

    # Uranus satellites
    for i in 1:length(URANUS_SATELLITE_NAMES)
        id = 700 + i
        add_edge!(BODIES, 7, id)
    end

    # Neptune satellites
    for i in 1:length(NEPTUNE_SATELLITE_NAMES)
        id = 800 + i
        add_edge!(BODIES, 8, id)
    end

    # Pluto satellites
    for i in 1:length(PLUTO_SATELLITE_NAMES)
        id = 900 + i
        add_edge!(BODIES, 9, id)
    end

    # Minor bodies
    add_edge!(BODIES, 0, 9)
    add_edge!(BODIES, 9, 999)
    add_edge!(BODIES, 0, 2000001)
    add_edge!(BODIES, 0, 2000002)
    add_edge!(BODIES, 0, 2000004)
    add_edge!(BODIES, 0, 2000016)
    add_edge!(BODIES, 0, 2000021)
    add_edge!(BODIES, 0, 2431010)
    add_edge!(BODIES, 0, 2000433)
    add_edge!(BODIES, 0, 2000511)
    add_edge!(BODIES, 0, 9511010)
    add_edge!(BODIES, 0, 2002867)
    add_edge!(BODIES, 0, 2025143)
    add_edge!(BODIES, 0, 1000093)
    add_edge!(BODIES, 0, 1000005)
end

end
