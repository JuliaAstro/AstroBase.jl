#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
module Bodies

import Base: parent

using ItemGraphs: ItemGraph, SimpleGraph, add_edge!, add_vertex!, items

using ..Time: SECONDS_PER_DAY, SECONDS_PER_CENTURY, value, j2000, seconds

export CelestialBody,
    NAIFId,
    SolarSystemBarycenter,
    Sun,
    along_orbit_radius,
    declination,
    declination_rate,
    ellipsoid,
    equatorial_radius,
    euler_angles,
    euler_rates,
    from_naifid,
    grav_param,
    mean_radius,
    naifid,
    parent,
    polar_radius,
    right_ascension,
    right_ascension_rate,
    rotation_angle,
    rotation_rate,
    ssb,
    subplanetary_radius,
    sun

const NAIFId = Int
const BODIES = ItemGraph{NAIFId, NAIFId}(SimpleGraph())

abstract type CelestialBody end

function grav_param end
function mean_radius end
function polar_radius end
function equatorial_radius end
function subplanetary_radius end
function along_orbit_radius end

function ellipsoid(body::CelestialBody)
    subplanetary_radius(body), along_orbit_radius(body), polar_radius(body)
end

Base.show(io::IO, body::CelestialBody) = print(io, string(nameof(typeof(body))))

from_naifid(id::NAIFId) = from_naifid(Val(id))

path_ids(from::CelestialBody, to::CelestialBody) = items(BODIES, naifid(from), naifid(to))

abstract type Barycenter <: CelestialBody end

struct SolarSystemBarycenter <: Barycenter end
const ssb = SolarSystemBarycenter()
Base.show(io::IO, ::SolarSystemBarycenter) = print(io, "Solar System Barycenter")
naifid(::SolarSystemBarycenter) = 0
from_naifid(::Val{0}) = ssb
add_vertex!(BODIES, 0)

struct Sun <: CelestialBody end
const sun = Sun()
parent(::Sun) = ssb
naifid(::Sun) = 10
from_naifid(::Val{10}) = sun
add_edge!(BODIES, 0, 10)

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

end
