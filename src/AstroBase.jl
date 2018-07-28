module AstroBase

using Rotations

export tio_locator, sec2rad, rad2sec, J2000, polar_motion, earth_rotation_angle,
       celestial_to_intermediate, nutation, normalize_pi_angle

const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const U2R = deg2rad(1/1e4 * (1/3600))

include("constants.jl")
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

"""
    normalize_pi_angle(a)

Returns an angle(radians) in the range [-π, π) for a given angle(radians).

# Example

```jldoctest
julia> normalize_pi_angle(21)
2.1504440784612413
```
"""
normalize_pi_angle(a) = rem2pi(a, RoundNearest)

"""
    nutation(jd1, jd2)

Returns nutation in longitude(radians) and obliquity(radians) for a given 2 part Julian date (TT format).

# Example

```jldoctest
julia> nutation(2.4578265e6, 0.30434616919175345)
(-3.7565297299394694e-5, -3.665617105048724e-5
```
"""
function nutation(jd1, jd2)

    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    mean_longitude_moon_minus_mean_longitude_moon_perigee = normalize_pi_angle(
        sec2rad(@evalpoly t 485866.733 715922.633 31.310 0.064 )
      + mod(1325.0 * t, 1.0) * 2pi)

    mean_longitude_sun_minus_mean_longitude_sun_perigee = normalize_pi_angle(
        sec2rad(@evalpoly t 1287099.804 1292581.224 -0.577 - 0.012)
      + mod(99.0 * t, 1.0) * 2pi)

    mean_longitude_moon_minus_mean_longitude_moon_node = normalize_pi_angle(
        sec2rad(@evalpoly t 335778.877 295263.137 -13.257 0.011)
      + mod(1342.0 * t, 1.0) * 2pi)

    mean_elongation_moon_from_sun = normalize_pi_angle(
        sec2rad(@evalpoly t 1072261.307 1105601.328 -6.891 0.019)
      + mod(1236.0 * t, 1.0) * 2pi)

    mean_ascending_node_lunar_orbit_ecliptic_measured_mean_equinox_date = normalize_pi_angle(
        sec2rad(@evalpoly t 450160.280 -482890.539 7.455 0.008)
      + mod(-5.0 * t, 1.0) * 2pi)

    dp = 0.0
    de = 0.0

    for i in reverse(eachindex(multiples_of_arguments_and_coefficients))
        arg = (multiples_of_arguments_and_coefficients[i][1]  * mean_longitude_moon_minus_mean_longitude_moon_perigee
        + multiples_of_arguments_and_coefficients[i][2] * mean_longitude_sun_minus_mean_longitude_sun_perigee
        + multiples_of_arguments_and_coefficients[i][3]  * mean_longitude_moon_minus_mean_longitude_moon_node
        + multiples_of_arguments_and_coefficients[i][4]  * mean_elongation_moon_from_sun
        + multiples_of_arguments_and_coefficients[i][5] * mean_ascending_node_lunar_orbit_ecliptic_measured_mean_equinox_date)

        s = multiples_of_arguments_and_coefficients[i][6] + multiples_of_arguments_and_coefficients[i][7] * t
        c = multiples_of_arguments_and_coefficients[i][8] + multiples_of_arguments_and_coefficients[i][9] * t
        sinarg, cosarg = sincos(arg)
        iszero(s) || (dp += s * sinarg)
        iszero(c) || (de += c * cosarg)
    end

    dp * U2R, de * U2R
end

end
