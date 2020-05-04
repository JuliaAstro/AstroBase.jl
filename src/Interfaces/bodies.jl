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

function alpha end
function alpha0 end
function alpha1 end
function delta end
function delta0 end
function delta1 end
function omega end
function omega0 end
function omega1 end
function theta0 end
function theta1 end

function ellipsoid(body::CelestialBody)
    subplanetary_radius(body), along_orbit_radius(body), polar_radius(body)
end
