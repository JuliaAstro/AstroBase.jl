module AstroBase

using Rotations

export tio_locator, sec2rad, rad2sec, J2000, polar_motion, earth_rotation_angle,
  celestial_to_intermediate, greenwich_mean_sidereal_time82, greenwich_mean_sidereal_time00, greenwich_mean_sidereal_time06

const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const SECONDS_PER_DAY = 24.0 * 60.0 * 60.0

"""
    celestial_to_intermediate(x, y, s)

Returns celestial to intermediate-frame-of-date transformation matrix given
the Celestial Intermediate Pole location (`x`, `y` and the CIO locator `s`).

```jldoctest
julia> celestial_to_intermediate(0.2, 0.2, 0.1)
3×3 RotZYZ{Float64}(0.785398, 0.286757, -0.885398):
  0.976728   0.0774803  0.2
 -0.11811    0.972651   0.2
 -0.179034  -0.218968   0.959166
```
"""
function celestial_to_intermediate(x, y, s)
    r2 = x^2 + y^2
    e = r2 > 0.0 ? atan2(y,x) : 0.0
    d = atan(sqrt(r2 / (1.0 - r2)))
    RotZYZ(e, d, -(e + s))
end

"""
    polar_motion(rx, ry, sp)

Form the matrix of polar motion for coordinates of the pole (radians).

# Example
```jldoctest
julia> polar_motion(20, 30, 50)
3×3 RotZYX{Float64}(50.0, -20.0, -30.0):
  0.393785  -0.829946  -0.395124
 -0.10707    0.385514  -0.916469
  0.912945   0.403198   0.0629472
```
"""
function polar_motion(rx, ry, sp)
    RotZYX{Float64}(sp, -rx, -ry)
end

"""
    earth_rotation_angle(jd1, jd2)

Return Earth rotation angle (radians) for a given UT1 2-part Julian Date (jd1, jd2).

# Example
```jldoctest
julia> earth_rotation_angle(2.4578265e6, 0.30434616919175345)
4.912208135094597
```
"""
function earth_rotation_angle(jd1, jd2)
    if jd1 < jd2
        d1 = jd1
        d2 = jd2
    else
        d1 = jd2
        d2 = jd1
    end
    t = d1 + d2 - J2000
    f = mod(d1, 1.0) + mod(d2, 1.0)
    mod2pi(2pi * (f + 0.7790572732640 + 0.00273781191135448 * t))
end

"""
    sec2rad(sec::Real)

Convert an angle in arcseconds to radians.

# Example

```jldoctest
julia> sec2rad(3600 * 30)
0.5235987755982988
```
"""
sec2rad(sec::Real) = deg2rad(sec / 3600)

"""
    rad2sec(rad::Real)
Convert an angle in radians to arcseconds.

# Example

```jldoctest
julia> rad2sec(0.5235987755982988)
107999.99999999999
```
"""
rad2sec(rad::Real) = rad2deg(rad) * 3600

"""
    tio_locator(jd1, jd2)

Returns TIO locator s' position for a given TT 2-part Julian date (jd1, jd2).

# Example

'''jldoctest
julia> AstroBase.tio_locator(2.4578265e6, 0.30434616919175345)
-3.9189245827947945e-11
'''
"""
function tio_locator(jd1, jd2)
    t = (jd1 - J2000 + jd2) / DAYS_PER_CENTURY
    -47e-6 * t * sec2rad(1)
end

function greenwich_mean_sidereal_time82(jd1, jd2)
    A = 24110.54841  -  DAYS_PER_CENTURY / 2.0
    B = 8640184.812866
    C = 0.093104
    D = -6.2e-6

    if jd1 < jd2
        d1 = jd1
        d2 = jd2
    else
        d1 = jd2
        d2 = jd1
    end
    t = (d1 + (d2 - J2000)) / DAYS_PER_CENTURY

    f = SECONDS_PER_DAY * (mod(d1, 1.0) + mod(d2, 1.0))

    mod2pi(sec2rad(((A + (B + (C + D * t) * t) * t) + f)))
end

"""
    greenwich_mean_sidereal_time00(ut1, ut2, tt1, tt2)

Returns Greenwich mean sidereal time(radians) for given two, 2 part Julian dates (TT and UT1).
(consistent with IAU 2000 precession)

# Example

```jldoctest
julia> greenwich_mean_sidereal_time00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
4.9596733720586075
```
"""
function greenwich_mean_sidereal_time00(ut1, ut2, tt1, tt2)
    t = ((tt1 - J2000) + tt2) / DAYS_PER_CENTURY
    mod2pi(earth_rotation_angle(ut1, ut2) + sec2rad(@evalpoly t 0.014506 4612.15739966 1.39667721 -0.00009344 0.00001882))
end

"""
    greenwich_mean_sidereal_time06(ut1, ut2, tt1, tt2)

Returns Greenwich mean sidereal time(radians) for given two, 2 part Julian dates (TT and UT1).
(consistent with IAU 2006 precession)

# Example

```jldoctest
julia> greenwich_mean_sidereal_time06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
4.959673370568533
```
"""
function greenwich_mean_sidereal_time06(ut1, ut2, tt1, tt2)
    t = ((tt1 - J2000) + tt2) / DAYS_PER_CENTURY
    mod2pi(earth_rotation_angle(ut1, ut2) + sec2rad(@evalpoly t 0.014506 4612.156534 1.3915817 -0.00000044 -0.000029956 -0.0000000368 ))
end
end
