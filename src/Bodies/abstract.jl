#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
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

# References

- [NASA NAIF](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html)
"""
const NAIFId = Int

"""
    grav_param([T=Float64,] body::CelestialBody)

Return the gravitational parameter ``\\mu = GM`` of `body` in km^3/s^2.

# Example

```jldoctest
julia> grav_param(earth)
398600.435436096
```

# References

- [NASA NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc)
"""
function grav_param end

"""
    mean_radius([T=Float64,] body::CelestialBody)

Return the mean radius of `body` in km.

# Example

```jldoctest
julia> mean_radius(earth)
6371.008366666666
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function mean_radius end

"""
    polar_radius([T=Float64,] body::CelestialBody)

Return the polar radius of `body` in km.

# Example

```jldoctest
julia> polar_radius(earth)
6356.7519
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function polar_radius end

"""
    equatorial_radius([T=Float64,] body::CelestialBody)

Return the equatorial radius of `body` in km.

!!! note
    The equatorial radius is only available for bodies which have identical subplanetary
    and along-orbit radii.

# Example

```jldoctest
julia> equatorial_radius(earth)
6378.1366
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function equatorial_radius end

"""
    subplanetary_radius([T=Float64,] body::CelestialBody)

Return the subplanetary radius of `body` in km.

# Example

```jldoctest
julia> subplanetary_radius(earth)
6378.1366
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function subplanetary_radius end

"""
    along_orbit_radius([T=Float64,] body::CelestialBody)

Return the along-orbit radius of `body` in km.

# Example

```jldoctest
julia> along_orbit_radius(earth)
6378.1366
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function along_orbit_radius end

"""
    ellipsoid([T=Float64,] body::CelestialBody)

Return the subplanetary, along-orbit, and polar radii of `body` in km which constitute its
tri-axial ellipsoid.

# Example

```jldoctest
julia> ellipsoid(earth)
(6378.1366, 6378.1366, 6356.7519)
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function ellipsoid(::Type{T}, body) where {T}
    subplanetary_radius(T, body), along_orbit_radius(T, body), polar_radius(T, body)
end
ellipsoid(body) = ellipsoid(Float64, body)

"""
    right_ascension([[T=Float64,] PT=Float64,], body, ep)

Return the right ascension of the body-fixed frame of `body` w.r.t. the ICRF at epoch `ep`
in rad.

# Example

```jldoctest
julia> right_ascension(earth, TTEpoch(2000, 1, 1))
1.5314935667739428e-7
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function right_ascension end

"""
    right_ascension_rate([[T=Float64,] PT=Float64,], body, ep)

Return the right ascension rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s.

# Example

```jldoctest
julia> right_ascension_rate(earth, TTEpoch(2000, 1, 1))
-3.545123997161905e-12
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function right_ascension_rate end

"""
    declination([[T=Float64,] PT=Float64,], body, ep)

Return the declination of the body-fixed frame of `body` w.r.t. the ICRF at epoch `ep` in
rad.

# Example

```jldoctest
julia> declination(earth, TTEpoch(2000, 1, 1))
1.5707964598747588
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function declination end

"""
    declination_rate([[T=Float64,] PT=Float64,], body, ep)

Return the declination rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s.

# Example

```jldoctest
julia> declination_rate(earth, TTEpoch(2000, 1, 1))
-3.0805523657085508e-12
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function declination_rate end

"""
    rotation_angle([[T=Float64,] PT=Float64,], body, ep)

Return the rotation angle of the body-fixed frame of `body` w.r.t. the ICRF at epoch `ep` in
rad.

# Example

```jldoctest
julia> rotation_angle(earth, TTEpoch(2000, 1, 1))
0.16849737156984945
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function rotation_angle end

"""
    rotation_rate([[T=Float64,] PT=Float64,], body, ep)

Return the rotation rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s.

# Example

```jldoctest
julia> rotation_rate(earth, TTEpoch(2000, 1, 1))
7.292115373194001e-5
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function rotation_rate end

"""
    euler_angles([[T=Float64,] PT=Float64,], body, ep)

Return the orientation of the body-fixed frame of `body` w.r.t. the ICRF at epoch `ep` as a
set of Euler angles in rad with rotation order *ZXZ*.

# Example

```jldoctest
julia> euler_angles(earth, TTEpoch(2000, 1, 1))
(1.5707964799442533, -1.3307986224120327e-7, 0.16849737156984945)
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function euler_angles(b::CelestialBody, ep)
    right_ascension(b, ep) + π/2, π/2 - declination(b, ep), mod2pi(rotation_angle(b, ep))
end

"""
    euler_rates([[T=Float64,] PT=Float64,], body, ep)

Return the rate of change of orientation of the body-fixed frame of `body` w.r.t. the ICRF
at epoch `ep` as a set of Euler angle rates in rad/s with rotation order *ZXZ*.

# Example

```jldoctest
julia> euler_rates(earth, TTEpoch(2000, 1, 1))
(-3.545123997161905e-12, 3.0805523657085508e-12, 7.292115373194001e-5)
```

# References

- Archinal, Brent Allen, et al. "Report of the IAU working group on cartographic coordinates
    and rotational elements: 2009." *Celestial Mechanics and Dynamical Astronomy* 109.2
    (2011): 101-135.
"""
function euler_rates(b::CelestialBody, ep)
    right_ascension_rate(b, ep), -declination_rate(b, ep), rotation_rate(b, ep)
end

function right_ascension_coeffs end
function declination_coeffs end
function rotation_coeffs end
function nutation_precession_coeffs end

@inline function theta(::Type{T}, body, t) where T
    θ₀, θ₁ = nutation_precession_coeffs(T, body)
    return θ₀ .+ θ₁ .* t ./ SECONDS_PER_CENTURY
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
        @inline function $func(::Type{T}, ::Type{PT}, body, ep) where {T, PT}
            t = julian_period(PT, ep; unit=seconds)
            c₀, c₁, c₂, c = $coeffs(T, body)
            c_np = isempty(c) ? zero(T) : sum(c .* $trig.(theta(T, body, t)))
            return c₀ + c₁ * t / $dt + c₂ * t^2 / $dt^2 + c_np
        end
        $func(body, ep) = $func(Float64, Float64, body, ep)

        @inline function $rate(::Type{T}, ::Type{PT}, body, ep) where {T, PT}
            t = julian_period(PT, ep; unit=seconds)
            _, c₁, c₂, c = $coeffs(T, body)
            _, θ₁ = nutation_precession_coeffs(T, body)
            c_np = isempty(c) ? zero(T) :
                sum(c .* θ₁ ./ SECONDS_PER_CENTURY .* $trigrate.(theta(T, body, t)))
            return c₁ / $dt + 2c₂ * t / $dt^2 + $sign * c_np
        end
        $rate(body, ep) = $rate(Float64, Float64, body, ep)
    end
end

