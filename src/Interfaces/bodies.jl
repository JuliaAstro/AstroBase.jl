#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

export CelestialBody, Barycenter, Planet, NaturalSatellite, MinorBody, NAIFId
export grav_param, mean_radius, polar_radius, equatorial_radius, subplanetary_radius,
    along_orbit_radius, alpha, alpha0, alpha1, delta, delta0, delta1, omega, omega0,
    omega1, theta0, theta1, ellipsoid


"""
    CelestialBody

Abstract supertype for all celestial bodies and pseudo-bodies.
"""
abstract type CelestialBody end

"""
    Barycenter <: CelestialBody

Abstract supertype for representing the barycenters of the solar system and planetary
systems as pseudo-bodies.
"""
abstract type Barycenter <: CelestialBody end

"""
    Planet <: CelestialBody

Abstract supertype for planets.
"""
abstract type Planet <: CelestialBody end

"""
    NaturalSatellite <: CelestialBody

Abstract supertype for natural satellites (moons).
"""
abstract type NaturalSatellite <: CelestialBody end

"""
    MinorBody <: CelestialBody

Abstract supertype for minor solar system bodies.
"""
abstract type MinorBody <: CelestialBody end

"""
    NAIFId

An integer code for identifying celestial bodies and other objects in space.
"""
const NAIFId = Int

"""
    grav_param([T=Float64,] body::CelestialBody)

Return the gravitational parameter ``\\mu = GM`` for `body` in km^3/s^2.
"""
function grav_param end
function mean_radius end
function polar_radius end
function equatorial_radius end
function subplanetary_radius end
function along_orbit_radius end

function right_ascension_coeffs end
function declination_coeffs end
function rotation_coeffs end
function nutation_precession_coeffs end

const Δtc = SECONDS_PER_CENTURY
const Δtd = SECONDS_PER_DAY

@inline function theta(::Type{T}, body, t) where T
    θ₀, θ₁ = nutation_precession_coeffs(T, body)
    return θ₀ .+ θ₁ .* t ./ Δtc
end

for func in (:right_ascension, :declination, :rotation)
    coeffs = Symbol(func, "_coeffs")
    rate = Symbol(func, "_rate")
    trig = func == :declination ? :cos : :sin
    trigrate = func == :declination ? :sin : :cos
    sign = func == :declination ? -1.0 : 1.0
    dt = func == :rotation ? :SECONDS_PER_DAY : :SECONDS_PER_CENTURY
    func = func == :rotation ? :rotation_angle : func
    @eval begin
        @inline function $func(::Type{T}, ::Type{TT}, body, ep) where {T, TT}
            t = julian_period(TT, ep; unit=seconds)
            c₀, c₁, c₂, c = $coeffs(T, body)
            c_np = isempty(c) ? zero(T) : sum(c .* $trig.(theta(T, body, t)))
            return c₀ + c₁ * t / $dt + c₂ * t^2 / $dt^2 + c_np
        end
        $func(body, ep) = $func(Float64, Float64, body, ep)

        @inline function $rate(::Type{T}, ::Type{TT}, body, ep) where {T, TT}
            t = julian_period(TT, ep; unit=seconds)
            _, c₁, c₂, c = $coeffs(T, body)
            _, θ₁ = nutation_precession_coeffs(T, body)
            c_np = isempty(c) ? zero(T) :
                sum(c .* θ₁ ./ SECONDS_PER_CENTURY .* $trigrate.(theta(T, body, t)))
            return c₁ / $dt + 2c₂ * t / $dt^2 + $sign * c_np
        end
        $rate(body, ep) = $rate(Float64, Float64, body, ep)
    end
end

function euler_angles(b::CelestialBody, ep)
    right_ascension(b, ep) + π/2, π/2 - declination(b, ep), mod2pi(rotation_angle(b, ep))
end

function euler_rates(b::CelestialBody, ep)
    right_ascension_rate(b, ep), -declination_rate(b, ep), rotation_rate(b, ep)
end

function ellipsoid(body::CelestialBody)
    subplanetary_radius(body), along_orbit_radius(body), polar_radius(body)
end

