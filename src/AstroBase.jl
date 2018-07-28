module AstroBase

using Reexport

include(joinpath("bodies", "Bodies.jl"))

@reexport using .Bodies

using Rotations

export tio_locator,
    sec2rad,
    rad2sec,
    J2000,
    polar_motion,
    earth_rotation_angle,
    celestial_to_intermediate,
    earth_rotation_angle,
    xy06,
    obliquity_of_ecliptic_06,
    mean_obliquity_of_ecliptic,
    precession_fukushima_williams06,
    nutation_00b

const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const ARCSECONDS_IN_CIRCLE = 1296000.0
const U2R = deg2rad(1/1e7 * (1/3600))
const DPPLAN = -0.135 * deg2rad(1/1e3 * (1/3600))
const DEPLAN =  0.388 * deg2rad(1/1e3 * (1/3600))

include("mfals.jl")
include("nut_const.jl")

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
    e = r2 > 0.0 ? atan(y, x) : 0.0
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
    sec2rad(sec)

Convert an angle in arcseconds to radians.

# Example

```jldoctest
julia> sec2rad(3600 * 30)
0.5235987755982988
```
"""
sec2rad(sec) = deg2rad(sec / 3600)

"""
    rad2sec(rad)
Convert an angle in radians to arcseconds.

# Example

```jldoctest
julia> rad2sec(0.5235987755982988)
107999.99999999999
```
"""
rad2sec(rad) = rad2deg(rad) * 3600

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
    mean_obliquity_of_ecliptic(jd1, jd2)

Returns  obliquity of the ecliptic (radians) for a given Julian 2 part date (TT).

# Example

```jldoctest
julia> mean_obliquity_of_ecliptic(2.4578265e6, 0.30434616919175345)
0.40905376936136706
```
"""
function mean_obliquity_of_ecliptic(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY
    sec2rad(@evalpoly t 84381.448 -46.8150 0.00059 0.001813)
end

"""
    obliquity_of_ecliptic_06(jd1, jd2)

Returns obliquity of ecliptic (radians) for a given 2 part Julain date (TT).

# Example

```jldoctest
julia> obliquity_of_ecliptic_06(2.4578265e6, 0.30434616919175345)
0.409053547482157
```
"""
function obliquity_of_ecliptic_06(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY
    sec2rad(@evalpoly t 84381.406 -46.836769 -0.0001831 0.00200340 -0.000000576 -0.0000000434)
end

"""
    precession_fukushima_williams06(jd1, jd2)

Returns fukushima angles(radians) for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> precession_fukushima_williams06(2.4578265e6, 0.30434616919175345)
(8.616170933989655e-6, 0.4090536093366178, 0.004201176043952816, 0.409053547482157)
```
"""
function precession_fukushima_williams06(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    sec2rad(@evalpoly t -0.052928 10.556378 0.4932044 -0.00031238 -0.000002788 0.0000000260),
    sec2rad(@evalpoly t 84381.412819 -46.811016 0.0511268 0.00053289 -0.000000440 -0.0000000176),
    sec2rad(@evalpoly t -0.041775 5038.481484 1.5584175 -0.00018522 -0.000026452 -0.0000000148),
    obliquity_of_ecliptic_06(jd1, jd2)
end

"""
    mean_anomaly(::Luna, t)

Returns mean anomaly of Moon for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest

julia> mean_anomaly(moon, 23.0)
0.5891752616281019
```
"""
function mean_anomaly(::Luna, t)
    p = @evalpoly t 485868.249036 1717915923.2178 31.8792 0.051635 -0.00024470
    sec2rad(p % ARCSECONDS_IN_CIRCLE)
end

"""
    mean_anomaly(::Sun, t)

Returns mean anomaly of the Sun for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> mean_anomaly(sun, 23.0)
5.857396217361825
```
"""
function mean_anomaly(::Sun, t)
    p = @evalpoly t 1287104.793048 129596581.0481 -0.5532 0.000136 -0.00001149
    sec2rad(p % ARCSECONDS_IN_CIRCLE)
end

"""
    mean_longitude_minus_lan(::Luna, t)

Returnsmean longitude of the Moon minus mean longitude of the ascending node for Julian
centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> mean_longitude_minus_lan(23.0)
3.103138156410118
```
"""
function mean_longitude_minus_lan(::Luna, t)
    p = @evalpoly t 335779.526232 1739527262.8478 -12.7512 -0.001037 0.00000417
    sec2rad(p % ARCSECONDS_IN_CIRCLE)
end

"""
    mean_elongation(::Luna, t)

Returns mean elongation of the Moon from the Sun for Julian centuries since
J2000.0 in TDB.

# Example

```jldoctest
julia> mean_elongation(::Luna, 23.0)
2.8012040574296484
```
"""
function mean_elongation(::Luna, t)
    p = @evalpoly t 1072260.703692 1602961601.2090 -6.3706 0.006593 -0.00003169
    sec2rad(p % ARCSECONDS_IN_CIRCLE)
end

"""
    mean_longitude_ascending_node(::Luna, t)

Return fundamental argument for Julian centuries since J2000.0 in TDB.

# Example
```jldoctest
julia> mean_longitude_ascending_node(luna, 23.0)
4.904897783682109
```
"""
function mean_longitude_ascending_node(::Luna, t)
    p = @evalpoly t 450160.398036 -6962890.5431 7.4722 0.007702 -0.00005939
    sec2rad(p % ARCSECONDS_IN_CIRCLE)
end

"""
    mean_longitude(::Mercury, t)

Returns mean longitude of Mercury for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> mean_longitude(mercury, 23.0)
2.160150897150834
```
"""
mean_longitude(::Mercury, t) = mod2pi(4.402608842 + 2608.7903141574t)

"""
    mean_longitude(::Venus, t)

Returns mean longitude of Venus for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> mean_longitude(venus, 23.0)
0.9030394378238363
```
"""
mean_longitude(::Venus, t) = mod2pi(3.176146697 + 1021.3285546211t)

"""
    mean_longitude(::Earth, t)

Returns mean longitude of Earth for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> mean_longitude(earth, 23.0)
1.501718780251826
```
"""
mean_longitude(::Earth, t) = mod2pi(1.753470314 + 628.3075849991t)

"""
    mean_longitude(::Mars, t)

Returns mean longitude of Mars for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> mean_longitude(mars, 23.0)
5.276431642365657
```
"""
mean_longitude(::Mars, t) = mod2pi(6.203480913 + 334.0612426700t)

"""
    mean_longitude(::Jupiter, t)

Returns mean longitude of Jupiter for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> mean_longitude(jupiter, 23.0)
6.233996285639864
```
"""
mean_longitude(::Jupiter, t) = mod2pi(0.599546497 + 52.9690962641t)

"""
    mean_longitude(::Saturn, t)

Returns mean longitude of Saturn for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> mean_longitude(saturn, 23.0)
1.3735042049922535
```
"""
mean_longitude(::Saturn, t) = mod2pi(0.874016757 + 21.3299104960t)

"""
    mean_longitude(::Uranus, t)

Returns mean longitude of Uranus for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> mean_longitude(uranus, 23.0)
1.5497819750715893
```
"""
mean_longitude(::Uranus, t) = mod2pi(5.481293872 + 7.4781598567t)

"""
    mean_longitude(::Neptune, t)

Returns mean longitude of Neptune for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> mean_longitude(neptune, 23.0)
5.053273953885775
```
"""
mean_longitude(::Neptune, t) = mod2pi(5.311886287 + 3.8133035638t)

"""
    general_precession_in_longitude(t)

Returns general accumulated precession in longitude for Julian centuries since J2000.0 in TDB.

# Example

```jldoctest
julia> general_precession_in_longitude(23.0)
0.56362992539
```
"""
general_precession_in_longitude(t) = @evalpoly t 0.00000538691 0.024381750

"""
    xy06(jd1, jd2)

Returns X, Y coordinates of celestial intermediate pole for a given 2-part Julian date in TT.

# Example

```jldoctest
julia> xy06(2.4578265e6, 0.30440190993249416)
(0.0016558850230577835, -3.986943362456243e-5)
```
"""
function xy06(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY
    # Powers of T.
    pt = [t^i for i = 0:6]

    fa = (mean_anomaly(luna, t),
          mean_anomaly(sun, t),
          mean_longitude_minus_lan(luna, t),
          mean_elongation(luna, t),
          mean_longitude_ascending_node(luna, t),
          mean_longitude(mercury, t),
          mean_longitude(venus, t),
          mean_longitude(earth, t),
          mean_longitude(mars, t),
          mean_longitude(jupiter, t),
          mean_longitude(saturn, t),
          mean_longitude(uranus, t),
          mean_longitude(neptune, t),
          general_precession_in_longitude(t))

    xpr = @evalpoly t x_coeff[1] x_coeff[2] x_coeff[3] x_coeff[4] x_coeff[5] x_coeff[6]
    ypr = @evalpoly t y_coeff[1] y_coeff[2] y_coeff[3] y_coeff[4] y_coeff[5] y_coeff[6]

    xypl = zeros(2)
    xyls = zeros(2)

    ialast = length(amp)
    for ifreq in reverse(eachindex(planetary))
        arg = 0.0
        for i in eachindex(fa)
            m = planetary[ifreq][i]
            arg += float(m) * fa[i]
        end

        sc = sincos(arg)

        ia = pointers_to_amp[ifreq + length(luni_solar)]
        for i in ialast:-1:ia
            j = i - ia + 1
            jxy = jaxy[j]
            jsc = jasc[j]
            jpt = japt[j]
            xypl[jxy] += amp[i] * sc[jsc] * pt[jpt]
        end
        ialast = ia - 1
    end

    for ifreq in reverse(eachindex(luni_solar))
        arg = 0.0
        for i in 1:5
           m = luni_solar[ifreq][i]
           arg += float(m) * fa[i]
        end

        sc = sincos(arg)

        ia = pointers_to_amp[ifreq]
        for i in ialast:-1:ia
            j = i - ia + 1
            jxy = jaxy[j]
            jsc = jasc[j]
            jpt = japt[j]
            xyls[jxy] += amp[i] * sc[jsc] * pt[jpt]
        end
        ialast = ia - 1
    end

    sec2rad((xpr + (xyls[1] + xypl[1]) / 1e6)), sec2rad(ypr + (xyls[2] + xypl[2]) / 1e6)
end

"""
    nutation_00b(jd1, jd2)

Returns luni-solar and planetary nutation for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> nutation_00b(2.4578265e6, 0.30440190993249416)
(-3.7589177912131684e-5, -3.6657431214029895e-5)
```
"""
function nutation_00b(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    el  = sec2rad(mod(485868.249036 + (1717915923.2178) * t, ARCSECONDS_IN_CIRCLE))
    elp = sec2rad(mod(1287104.79305 + (129596581.0481) * t, ARCSECONDS_IN_CIRCLE))
    f   = sec2rad(mod(335779.526232 + (1739527262.8478) * t, ARCSECONDS_IN_CIRCLE))
    d   = sec2rad(mod(1072260.70369 + (1602961601.2090) * t, ARCSECONDS_IN_CIRCLE))
    om  = sec2rad(mod(450160.398036 + (-6962890.5431) * t, ARCSECONDS_IN_CIRCLE))

    dp = 0.0
    de = 0.0

    for i in reverse(eachindex(x_nutation))
        arg = mod(x_nutation[i][1]  * el  +
              x_nutation[i][2] * elp +
              x_nutation[i][3]  * f   +
              x_nutation[i][4]  * d   +
              x_nutation[i][5] * om, 2pi)
        sarg, carg = sincos(arg)

        dp += (x_nutation[i][6] + x_nutation[i][7] * t) * sarg + x_nutation[i][8] * carg
        de += (x_nutation[i][9] + x_nutation[i][10] * t) * carg + x_nutation[i][11] * sarg
    end
    dpsils = dp * U2R
    depsls = de * U2R

    dpsipl = DPPLAN
    depspl = DEPLAN

    dpsils + dpsipl, depsls + depspl

end

end # module
