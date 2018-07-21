module AstroBase

using Rotations

export tio_locator, sec2rad, rad2sec, J2000, polar_motion, earth_rotation_angle,
  celestial_to_intermediate, precession_rate_part_of_nutation, bias_precession_matrix_00

const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0

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
    precession_rate_part_of_nutation(jd1, jd2)

Returns precession corrections for a given 2 part Julian date (TT).

# Example

```jldoctest

```
"""
function precession_rate_part_of_nutation(jd1, jd2)
    PRECESSION = -sec2rad(0.29965)
    OBLIQUITY = -sec2rad(0.02524)

    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY
    PRECESSION * t, OBLIQUITY * t
end

function bias_component_00()
    DPBIAS = sec2rad(-0.041775)
    DEBIAS = sec2rad(-0.0068192)
    DRA0 = sec2rad(-0.0146)

    DPBIAS, DEBIAS, DRA0
end

"""
    bias_precession_matrix_00(jd1, jd2)

Returns a tuple of array, rb (frame bias), rp (precession matrix), rbp (bias precession matrix)
given 2 part Julian date (TT).

# Example

```jldoctest
julia> bias_precession_matrix_00(2.4578265e6, 0.30434616919175345)
(RotZYX(-7.07828e-8, -8.05622e-8, 3.30604e-8), [0.999991 0.00384587 0.00167106; -0.00384587 0.999993 -3.2343e-6; -0.00167106 -3.1924e-6 0.999999], [0.999991 0.00384594 0.00167098; -0.00384594 0.999993 -3.26705e-6; -0.00167098 -3.15946e-6 0.999999])
```
"""
function bias_precession_matrix_00(jd1, jd2)
    EPS0 = sec2rad(84381.448)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    dpsibi, depsbi, dra0 = bias_component_00()

    psia77 = sec2rad((@evalpoly t 0 5038.7784 -1.07259 -0.001147))
    oma77  = EPS0 + sec2rad(((0.05127 + (-0.007726) * t) * t) * t)
    chia   = sec2rad((@evalpoly t 0 10.5526 -2.38064 -0.001125))

    dpsipr, depspr = precession_rate_part_of_nutation(jd1, jd2)
    psia = psia77 + dpsipr
    oma  = oma77 + depspr

    rbw = RotZYX(dra0, dpsibi*sin(EPS0), -depsbi)
    rb = copy(rbw)
    rp = RotXZX(EPS0, -psia, -oma) * RotZ(chia)

    rb, rp, rp * rbw
end
end
