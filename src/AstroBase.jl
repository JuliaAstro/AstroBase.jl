module AstroBase

using Reexport

include(joinpath("bodies", "Bodies.jl"))

@reexport using .Bodies

using Rotations

export tio_locator, sec2rad, rad2sec, J2000, polar_motion, earth_rotation_angle,
        celestial_to_intermediate, earth_rotation_angle, xy06, mean_anomaly_of_moon, mean_anomaly_of_sun, mean_longitude_of_moon_minus_mean_longitude_of_ascending_node,
         mean_elongation_of_moon_from_sun, mean_longitude_ascending_node_moon, mean_longitude_of_mercury, mean_longitude_of_venus,
         mean_longitude_of_earth, mean_longitude_of_mars, mean_longitude_of_jupiter, mean_longitude_of_saturn,
         mean_longitude_of_uranus, mean_longitude_of_neptune, general_precession_in_longitude

const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const TURNAS = 1296000.0

include("mfals.jl")

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
    mean_anomaly_of_moon(t)

Returns mean anomaly of Moon for a given Julian century.

# Example

```jldoctest
julia> mean_anomaly_of_moon(23.0)
0.5891752616281019
```
"""
mean_anomaly_of_moon(t)  =  sec2rad(mod((@evalpoly t 485868.249036 1717915923.2178 31.8792 0.051635 -0.00024470) , TURNAS))

"""
    mean_anomaly_of_sun(t)

Returns mean anomaly of Sun for a given Julian century.

# Example

```jldoctest
julia> mean_anomaly_of_sun(23.0)
5.857396217361825
```
"""
mean_anomaly_of_sun(t) =  sec2rad(mod((@evalpoly t 1287104.793048 129596581.0481 -0.5532 0.000136 -0.00001149),  TURNAS))

"""
    mean_longitude_of_moon_minus_mean_longitude_of_ascending_node(t)

Returns mean longitude of the Moon for a given Julian century.

# Example

```jldoctest
julia> mean_longitude_of_moon_minus_mean_longitude_of_ascending_node(23.0)
3.103138156410118
```
"""
mean_longitude_of_moon_minus_mean_longitude_of_ascending_node(t)  =  sec2rad(mod((@evalpoly t 335779.526232 1739527262.8478 -12.7512 -0.001037 0.00000417), TURNAS))

"""
    mean_elongation_of_moon_from_sun(t)

Returns mean elongation of the Moon from the Sun for given Julian century.

# Example

```jldoctest
julia> mean_elongation_of_moon_from_sun(23.0)
2.8012040574296484
```
"""
mean_elongation_of_moon_from_sun(t)  =  sec2rad(mod((@evalpoly t 1072260.703692 1602961601.2090 -6.3706 0.006593 -0.00003169), TURNAS))

"""
    mean_longitude_ascending_node_moon(t)

Return fundamental argument for a given Julian century.

# Example
```jldoctest
julia> mean_longitude_ascending_node_moon(23.0)
4.904897783682109
```
"""
mean_longitude_ascending_node_moon(t) =  sec2rad(mod((@evalpoly t 450160.398036 -6962890.5431 7.4722 0.007702 -0.00005939), TURNAS ))

"""
    mean_longitude_of_mercury(t)

Returns mean longitude of Mercury for a given Julian century.

# Example

```jldoctest
julia> mean_longitude_of_mercury(23.0)
2.160150897150834
```
"""
mean_longitude_of_mercury(t) =  mod2pi(4.402608842 + 2608.7903141574t)

"""
    mean_longitude_of_venus(t)

Returns mean longitude of Venus for a given Julian century.

# Example

```jldoctest
julia> mean_longitude_of_venus(23.0)
0.9030394378238363
```
"""
mean_longitude_of_venus(t) =  mod2pi(3.176146697 + 1021.3285546211t)

"""
    mean_longitude_of_earth(t)

Returns mean longitude of Earth for a given Julian century.

# Example

```jldoctest
julia> mean_longitude_of_earth(23.0)
1.501718780251826
```
"""
mean_longitude_of_earth(t)  =  mod2pi(1.753470314 + 628.3075849991t)

"""
    mean_longitude_of_mars(t)

Returns mean longitude of Mars for a given Julian century.

# Example

```jldoctest
julia> mean_longitude_of_mars(23.0)
5.276431642365657
```
"""
mean_longitude_of_mars(t) =  mod2pi(6.203480913 + 334.0612426700t)

"""
    mean_longitude_of_jupiter(t)

Returns mean longitude of Jupiter for a given Julian century.

# Example

```jldoctest
julia> mean_longitude_of_jupiter(23.0)
6.233996285639864
```
"""
mean_longitude_of_jupiter(t) =  mod2pi(0.599546497 + 52.9690962641t)

"""
    mean_longitude_of_saturn(t)

Returns mean longitude of Saturn for a given Julian century.

# Example

```jldoctest
julia> mean_longitude_of_saturn(23.0)
1.3735042049922535
```
"""
mean_longitude_of_saturn(t) =  mod2pi(0.874016757 + 21.3299104960t)

"""
    mean_longitude_of_uranus(t)

Returns mean longitude of Uranus for a given Julian century.

# Example

```jldoctest
julia> mean_longitude_of_uranus(23.0)
1.5497819750715893
```
"""
mean_longitude_of_uranus(t) =  mod2pi(5.481293872 + 7.4781598567t)

"""
    mean_longitude_of_neptune(t)

Returns mean longitude of Neptune for a given Julian century.

# Example

```jldoctest
julia> mean_longitude_of_neptune(23.0)
5.053273953885775
```
"""
mean_longitude_of_neptune(t) =  mod2pi(5.311886287 + 3.8133035638t)

"""
    general_precession_in_longitude(t)

Returns general accumulated precession in longitude for a given Julian century.

# Example

```jldoctest
julia> general_precession_in_longitude(23.0)
0.56362992539
```
"""
general_precession_in_longitude(t) =  (0.024381750 + 0.00000538691t) * t

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
    NFLS = length(fundamental_argument_multipliers)
    NFPL = length(fundamental_argument_multipliers_pl)
    NA = length(amp)
    MAXPT= 5
    pt = Vector{Float64}(MAXPT + 1)
    fa = Vector{Float64}(14)
    xypr, xypl, xyls, sc = zeros(2), zeros(2), zeros(2), zeros(2)

    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY
    # Powers of T.
    w = 1.0
    for i in 1:MAXPT+1
        pt[i] = w
        w *= t
    end

    fa[1]  = mean_anomaly_of_moon(t)
    fa[2]  = mean_anomaly_of_sun(t)
    fa[3]  = mean_longitude_of_moon_minus_mean_longitude_of_ascending_node(t)
    fa[4]  = mean_elongation_of_moon_from_sun(t)
    fa[5]  = mean_longitude_ascending_node_moon(t)
    fa[6]  = mean_longitude_of_mercury(t)
    fa[7]  = mean_longitude_of_venus(t)
    fa[8]  = mean_longitude_of_earth(t)
    fa[9]  = mean_longitude_of_mars(t)
    fa[10] = mean_longitude_of_jupiter(t)
    fa[11] = mean_longitude_of_saturn(t)
    fa[12] = mean_longitude_of_uranus(t)
    fa[13] = mean_longitude_of_neptune(t)
    fa[14] = general_precession_in_longitude(t)

    for i in 1:2
           xypr[i] = @evalpoly t xy_polynomial_coefficients[i][1] xy_polynomial_coefficients[i][2] xy_polynomial_coefficients[i][3] xy_polynomial_coefficients[i][4] xy_polynomial_coefficients[i][5] xy_polynomial_coefficients[i][6]
    end

    ialast = NA
    for ifreq in NFPL:-1:1
        arg = 0.0
        for i in range(1,14)
           m = fundamental_argument_multipliers_pl[ifreq][i]
           if (m != 0)
               arg += float(m) * fa[i]
           end
        end

        sc[2], sc[1] = reim(cis(arg))

        ia = pointers_to_amp[ifreq + NFLS]
        for i in (ialast + 1):-1:(ia + 1)
               j = i - ia
               jxy = jaxy[j]
               jsc = jasc[j]
               jpt = japt[j]
               xypl[jxy] += amp[i-1] * sc[jsc] * pt[jpt]
        end
        ialast = ia - 1
    end

    for ifreq in NFLS:-1:1
        arg = 0.0
        for i in 1:5
           m = fundamental_argument_multipliers[ifreq][i]
           if (m != 0)
               arg += float(m) * fa[i]
           end
        end

        sc[2], sc[1] = reim(cis(arg))

        ia = pointers_to_amp[ifreq]
        for i in (ialast + 1):-1:(ia + 1)
               j = i - ia
               jxy = jaxy[j]
               jsc = jasc[j]
               jpt = japt[j]
               xyls[jxy] += amp[i-1] * sc[jsc] * pt[jpt]
        end
        ialast = ia - 1
    end

    sec2rad((xypr[1] + (xyls[1] + xypl[1]) / 1e6)), sec2rad(xypr[2] + (xyls[2] + xypl[2]) / 1e6)
end
end # module
