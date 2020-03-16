module EarthAttitude

using ..Bodies
using ..Util: sec2rad
using ReferenceFrameRotations: angle_to_dcm, angleaxis_to_dcm, compose_rotation

export
    J2000,
    bias_precession_matrix_00,
    celestial_to_intermediate,
    earth_rotation_angle,
    equation_of_equinoxes_00,
    equation_of_equinoxes_complementary_terms,
    equation_of_origins,
    fukushima_williams_matrix,
    general_precession_in_longitude,
    greenwich_mean_sidereal_time00,
    greenwich_mean_sidereal_time06,
    greenwich_mean_sidereal_time82,
    mean_elongation,
    mean_longitude,
    mean_longitude_ascending_node,
    mean_longitude_minus_lan,
    mean_obliquity_of_ecliptic,
    numat,
    nutation,
    nutation_00a,
    nutation_00b,
    obliquity_of_ecliptic_06,
    polar_motion,
    precession_fukushima_williams06,
    precession_rate_part_of_nutation,
    s00,
    s06,
    tio_locator,
    xy06,
    nutation_a06,
    greenwich_apparent_sidereal_time06,
    equation_of_equinoxes_94,
    nutation_matrix80,
    precession_nutation00,
    precession_nutation_a00,
    precession_nutation_b00,
    precession_nutation_matrix_a00,
    precession_nutation_matrix_b00,
    precession_nutation_matrix_a06,
    s00a,
    s00b,
    nutation_matrix_day_a00,
    nutation_matrix_day_b00,
    num06a,
    nutation_matrix_day,
    xys00a,
    xys00b,
    s06a,
    equation_of_equinoxes_a00,
    equation_of_equinoxes_b00,
    greenwich_mean_sidereal_time_a06,
    equation_of_equinoxes_a06,
    greenwich_apparent_sidereal_time94,
    greenwich_apparent_sidereal_time_a00,
    greenwich_apparent_sidereal_time_b00,
    greenwich_apparent_sidereal_time_a06,
    celestial_to_intermediate_frame_of_date,
    celestial_to_intermediate_matrix,
    celestial_to_intermediate_matrix_a00,
    celestial_to_intermediate_matrix_b00,
    celestial_to_terrestrial_matrix,
    celestial_to_terrestrial_a00,
    celestial_to_terrestrial_b00,
    c2tpe

include("iau_models.jl")
include("nutation.jl")

const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const ARCSECONDS_IN_CIRCLE = 1296000.0
const SECONDS_PER_DAY = 24.0 * 60.0 * 60.0

include(joinpath("constants", "mfals.jl"))
include(joinpath("constants", "nut_const.jl"))
include(joinpath("constants", "cio_locator.jl"))
include(joinpath("constants", "EE00.jl"))
include(joinpath("constants", "NUTATION80.jl"))

"""
    celestial_to_intermediate(x, y, s)

Returns celestial to intermediate-frame-of-date transformation matrix given
the Celestial Intermediate Pole location (`x`, `y` and the CIO locator `s`).

# Example

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
    angle_to_dcm(e + s, -d, -e, :ZYZ)
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
    angle_to_dcm(ry, rx, -sp, :XYZ)
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
    tio_locator(jd1, jd2)

Returns TIO locator s' position for a given TT 2-part Julian date (jd1, jd2).

# Example

```jldoctest
julia> AstroBase.tio_locator(2.4578265e6, 0.30434616919175345)
-3.9189245827947945e-11
```
"""
function tio_locator(jd1, jd2)
    t = (jd1 - J2000 + jd2) / DAYS_PER_CENTURY
    sec2rad(-47e-6 * t)
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

# Example

Returns mean anomaly of the Sun for Julian centuries since J2000.0 in TDB.
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

Returns mean longitude of the Moon minus mean longitude of the ascending node for Julian
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
        arg = mod(x_nutation[i][1] * el
                  + x_nutation[i][2] * elp
                  + x_nutation[i][3] * f
                  + x_nutation[i][4] * d
                  + x_nutation[i][5] * om, 2pi)
        sarg, carg = sincos(arg)

        dp += (x_nutation[i][6] + x_nutation[i][7] * t) * sarg + x_nutation[i][8] * carg
        de += (x_nutation[i][9] + x_nutation[i][10] * t) * carg + x_nutation[i][11] * sarg
    end

    sec2rad(-0.135e-3) + sec2rad(1e-7dp), sec2rad(0.388e-3) + sec2rad(1e-7de)
end

 """
    fukushima_williams_matrix(gamb, phib, psi, eps)

Returns  obliquity of the ecliptic (radians) for a given Julian 2 part date (TT).

# Example

```jldoctest
julia> fukushima_williams_matrix(0.2,0.3,0.5,0.6)
3×3 RotMatrix{Float64}:
  0.951082   0.21718   0.219716
 -0.274534   0.920305  0.278692
 -0.14168   -0.325378  0.93491
```
"""
function fukushima_williams_matrix(gamb, phib, psi, eps)
    compose_rotation(angleaxis_to_dcm(eps, [1.0, 0.0, 0.0]), angle_to_dcm(psi, -phib, -gamb, :ZXZ))
end

"""
    numat(epsa, dpsi, deps)

Returns nutation matrix for a given epsa(mean obliquity), dpsi and deps nutation.

# Example

```jldoctest
julia> numat(0.7, 1.4, 1.3)
3×3 RotXZX{Float64}(0.7, -1.4, -2.0):
  0.169967  -0.410092   0.896067
 -0.753714   0.531687   0.386296
 -0.634844  -0.741035  -0.218722
```
"""
function numat(epsa, dpsi, deps)
    angle_to_dcm(epsa + deps, dpsi, -epsa, :XZX)
end

"""
    greenwich_mean_sidereal_time82(jd1, jd2)

Returns Greenwich mean sidereal time(radians) for given 2 part Julian dates (UT1).
(consistent with IAU 1982 model)

# Example

```jldoctest
julia> greenwich_mean_sidereal_time82(2.4578265e6, 0.30434616919175345)
4.916054244834956
```
"""
function greenwich_mean_sidereal_time82(jd1, jd2)
    A = 24110.54841  -  SECONDS_PER_DAY / 2.0
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
    mod2pi(7.272205216643039903848712e-5 * (@evalpoly t A + f B C D))
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
    mod2pi(earth_rotation_angle(ut1, ut2) +
           sec2rad(@evalpoly t 0.014506 4612.15739966 1.39667721 -0.00009344 0.00001882))
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
    mod2pi(earth_rotation_angle(ut1, ut2) +
           sec2rad(@evalpoly t 0.014506 4612.156534 1.3915817 -0.00000044 -0.000029956 -0.0000000368))
end

"""
    precession_rate_part_of_nutation(jd1, jd2)

Returns precession corrections for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> precession_rate_part_of_nutation(2400000.5, 53736)
-0.8716465172668347629e-7, -0.7342018386722813087e-8
```
"""
function precession_rate_part_of_nutation(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY
    -sec2rad(0.29965) * t, -sec2rad(0.02524) * t
end

"""
    bias_precession_matrix_00(jd1, jd2)

Returns a tuple of array, rb (frame bias), rp (precession matrix), rbp (bias precession matrix)
given 2 part Julian date (TT).

# Example

julia> bias_precession_matrix_00(2.4578265e6, 0.30434616919175345)
(RotZYX(-7.07828e-8, -8.05622e-8, 3.30604e-8), [0.999991 0.00384587 0.00167106; -0.00384587 0.999993 -3.2343e-6; -0.00167106 -3.1924e-6 0.999999], [0.999991 0.00384594 0.00167098; -0.00384594 0.999993 -3.26705e-6; -0.00167098 -3.15946e-6 0.999999])
```
"""
function bias_precession_matrix_00(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    dpbias = sec2rad(-0.041775)
    debias = sec2rad(-0.0068192)
    dra0 = sec2rad(-0.0146)
    eps0 = sec2rad(84381.448)

    psia77 = sec2rad((@evalpoly t 0 5038.7784 -1.07259 -0.001147))
    oma77  = eps0 + sec2rad(@evalpoly t 0 0 0.05127 -0.007726)
    chia   = sec2rad((@evalpoly t 0 10.5526 -2.38064 -0.001125))

    dpsipr, depspr = precession_rate_part_of_nutation(jd1, jd2)
    psia = psia77 + dpsipr
    oma  = oma77 + depspr

    rbw = angle_to_dcm(debias, -(dpbias * sin(eps0)), -dra0, :XYZ)
    rp = compose_rotation(angleaxis_to_dcm(-chia, [0.0, 0.0, 1.0]),
                          angle_to_dcm(oma, psia, -eps0, :XZX))

    rbw, rp, compose_rotation(rp, rbw)
end

"""
    equation_of_origins(rnpb, s)

Returns the equation of origins(radians) for given nutation-bias-precession matrix and the CIO locator.

 # Example

 ```jldoctest
equation_of_origins(rand(3,3), 0.2)
1.7738370040531068
 ```
"""
function equation_of_origins(rnpb, s)
    x = rnpb[1, 3]
    ax = x / (1.0 + rnpb[3, 3])
    xs = 1.0 - ax * x
    ys = -ax * rnpb[2, 3]
    zs = -x
    p = rnpb[1, 1] * xs + rnpb[2, 1] * ys + rnpb[3, 1] * zs
    q = rnpb[1, 2] * xs + rnpb[2, 2] * ys + rnpb[3, 2] * zs
    p != 0 || q != 0 ? s - atan(q, p) : s
end

"""
    nutation(jd1, jd2)

Returns nutation in longitude(radians) and obliquity(radians) for a given 2 part Julian date (TT format).

# Example

julia> nutation(2.4578265e6, 0.30434616919175345)
(-3.7565297299394694e-5, -3.665617105048724e-5
```
"""
function nutation(jd1, jd2)

    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    mean_longitude_moon_minus_mean_longitude_moon_perigee = rem2pi(
        sec2rad(@evalpoly t 485866.733 715922.633 31.310 0.064 )
      + mod(1325.0 * t, 1.0) * 2pi, RoundNearest)

    mean_longitude_sun_minus_mean_longitude_sun_perigee = rem2pi(
        sec2rad(@evalpoly t 1287099.804 1292581.224 -0.577 - 0.012)
      + mod(99.0 * t, 1.0) * 2pi, RoundNearest)

    mean_longitude_moon_minus_mean_longitude_moon_node = rem2pi(
        sec2rad(@evalpoly t 335778.877 295263.137 -13.257 0.011)
      + mod(1342.0 * t, 1.0) * 2pi, RoundNearest)

    mean_elongation_moon_from_sun = rem2pi(
        sec2rad(@evalpoly t 1072261.307 1105601.328 -6.891 0.019)
      + mod(1236.0 * t, 1.0) * 2pi, RoundNearest)

    mean_ascending_node_lunar_orbit_ecliptic_measured_mean_equinox_date = rem2pi(
        sec2rad(@evalpoly t 450160.280 -482890.539 7.455 0.008)
      + mod(-5.0 * t, 1.0) * 2pi, RoundNearest)

    dp = 0.0
    de = 0.0

    for i in reverse(eachindex(multiples_of_coefficients))
        arg = (multiples_of_coefficients[i][1]  * mean_longitude_moon_minus_mean_longitude_moon_perigee
        + multiples_of_coefficients[i][2] * mean_longitude_sun_minus_mean_longitude_sun_perigee
        + multiples_of_coefficients[i][3]  * mean_longitude_moon_minus_mean_longitude_moon_node
        + multiples_of_coefficients[i][4]  * mean_elongation_moon_from_sun
        + multiples_of_coefficients[i][5] * mean_ascending_node_lunar_orbit_ecliptic_measured_mean_equinox_date)

        s = multiples_of_arguments[i][1] + multiples_of_arguments[i][2] * t
        c = multiples_of_arguments[i][3] + multiples_of_arguments[i][4] * t
        sinarg, cosarg = sincos(arg)
        dp += s * sinarg
        de += c * cosarg
    end

    sec2rad(1e-4dp), sec2rad(1e-4de)
end

"""
    equation_of_equinoxes_complementary_terms(jd1, jd2)

Returns complementary terms for a given 2 part Julian date (TT).

# Example

julia> equation_of_equinoxes_complementary_terms(2.4578265e6, 0.30434616919175345)
5.706799075288604e-9
```
"""
function equation_of_equinoxes_complementary_terms(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    fa = (mean_anomaly(luna, t),
          mean_anomaly(sun, t),
          mean_longitude_minus_lan(luna, t),
          mean_elongation(luna, t),
          mean_longitude_ascending_node(luna, t),
          mean_longitude(venus, t),
          mean_longitude(earth, t),
          general_precession_in_longitude(t))

    s0 = 0.0
    s1 = 0.0

    for i in reverse(eachindex(e0_coefficent))
        a = 0.0
        for j in 1:8
            a += e0_coefficent[i][j] * fa[j]
        end
        s0 += e0_arg[i][1] * sin(a) + e0_arg[i][2] * cos(a)
    end

    for i in reverse(eachindex(e1_coefficent))
        a = 0.0;
        for j in 1:8
            a += e1_coefficent[i][j] * fa[j]
        end
        s1 += e1_arg[i][1] * sin(a) + e1_arg[i][2] * cos(a)
    end

    sec2rad(s0 + s1 * t)
end

"""
    equation_of_equinoxes_00(jd1, jd2, epsa, dpsi)

Return equation of equinoxes for given 2 part Julian date (TT), mean obliquity and nutation in longitude.

# Example

```jldoctest
julia> equation_of_equinoxes_00(2.4578265e6, 0.30440190993249416, 1.5, 1.7)
0.12025324854189404
```
"""
function equation_of_equinoxes_00(jd1, jd2, epsa, dpsi)
    dpsi * cos(epsa) + equation_of_equinoxes_complementary_terms(jd1, jd2)
end

function cio_locator(coeffs, terms0, terms1, terms2, terms3, terms4, jd1, jd2, x, y)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    fa = (
        mean_anomaly(luna, t),
        mean_anomaly(sun, t),
        mean_longitude_minus_lan(luna, t),
        mean_elongation(luna, t),
        mean_longitude_ascending_node(luna,t),
        mean_longitude(venus, t),
        mean_longitude(earth, t),
        general_precession_in_longitude(t),
    )

    w = collect(coeffs)

    for ct in reverse(terms0)
        a = 0.0
        for i in reverse(eachindex(fa))
            a += ct.coeffs[i] * fa[i]
        end
        s, c = sincos(a)
        w[1] += ct.sc * s + ct.cc * c
    end

    for ct in reverse(terms1)
        a = 0.0
        for i in reverse(eachindex(fa))
            a += ct.coeffs[i] * fa[i]
        end
        s, c = sincos(a)
        w[2] += ct.sc * s + ct.cc * c
    end

    for ct in reverse(terms2)
        a = 0.0
        for i in reverse(eachindex(fa))
            a += ct.coeffs[i] * fa[i]
        end
        s, c = sincos(a)
        w[3] += ct.sc * s + ct.cc * c
    end

    for ct in reverse(terms3)
        a = 0.0
        for i in reverse(eachindex(fa))
            a += ct.coeffs[i] * fa[i]
        end
        s, c = sincos(a)
        w[4] += ct.sc * s + ct.cc * c
    end

    for ct in reverse(terms4)
        a = 0.0
        for i in reverse(eachindex(fa))
            a += ct.coeffs[i] * fa[i]
        end
        s, c = sincos(a)
        w[5] += ct.sc * s + ct.cc * c
    end

    sec2rad((@evalpoly t w[1] w[2] w[3] w[4] w[5] w[6])) - x * y / 2.0
end

const s00(jd1, jd2, x, y) = cio_locator(CIO_COEFFS_00,
                                        CIO_TERMS_0TH_00,
                                        CIO_TERMS_1ST_00,
                                        CIO_TERMS_2ND_00,
                                        CIO_TERMS_3RD_00,
                                        CIO_TERMS_4TH_00,
                                        jd1, jd2, x, y)

const s06(jd1, jd2, x, y) = cio_locator(CIO_COEFFS_06,
                                        CIO_TERMS_0TH_06,
                                        CIO_TERMS_1ST_06,
                                        CIO_TERMS_2ND_06,
                                        CIO_TERMS_3RD_06,
                                        CIO_TERMS_4TH_06,
                                        jd1, jd2, x, y)

"""
    s00(jd1, jd2, x, y)

Returns Celestial Intermediate Origin(CIO) for a given 2-part Julian date (jd1, jd2)
with CIP coordinates (x, y). Compatible with IAU-2000 precession-nutation.

# Example

```jldoctest
julia> AstroBase.s00(2.4578265e6, 0.30434616919175345, 20, 50)
-500.00000000383193
```
"""
s00

"""
    s06(jd1, jd2, x, y)

Returns Celestial Intermediate Origin(CIO) for a given 2-part Julian date (jd1, jd2)
with CIP coordinates (x, y)

 # Example

julia> AstroBase.s06(2.4578265e6, 0.30434616919175345, 20, 50)
-500.00000000383193
```
"""
s06

"""
    nutation_00a(jd1, jd2)

Returns luni-solar and planetary nutations for a given 2 part Julian date(TT).

# Example

```jldoctest
julia> nutation_00a(2.4578265e6, 0.30440190993249416)
(-3.7589903391357206e-5, -3.6659049617818334e-5)
```
"""
function nutation_00a(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    el = mean_anomaly(luna, t)
    f  = mean_longitude_minus_lan(luna, t)
    om = mean_longitude_ascending_node(luna, t)

    apa  = general_precession_in_longitude(t)
    alme = mean_longitude(mercury, t)
    alve = mean_longitude(venus, t)
    alea = mean_longitude(earth, t)
    alma = mean_longitude(mars, t)
    alju = mean_longitude(jupiter, t)
    alsa = mean_longitude(saturn, t)
    alur = mean_longitude(uranus, t)

    elp = sec2rad(mod((@evalpoly t 1287104.79305 129596581.0481 -0.5532 0.000136 -0.00001149),
                      ARCSECONDS_IN_CIRCLE))
    d   = sec2rad(mod((@evalpoly t 1072260.70369 1602961601.2090 -6.3706 0.006593 -0.00003169),
                      ARCSECONDS_IN_CIRCLE))

    al   = mod2pi(2.35555598 + 8328.6914269554 * t)
    af   = mod2pi(1.627905234 + 8433.466158131 * t)
    ad   = mod2pi(5.198466741 + 7771.3771468121 * t)
    aom  = mod2pi(2.18243920 - 33.757045 * t)
    alne = mod2pi(5.321159000 + 3.8127774000 * t)

    dpls = 0.0
    dels = 0.0

    for i in reverse(eachindex(xls_nutation))
        arg = mod2pi(xls_nutation[i][1]  * el +
                    xls_nutation[i][2] * elp +
                    xls_nutation[i][3]  * f +
                    xls_nutation[i][4]  * d +
                    xls_nutation[i][5] * om)
        sarg, carg = sin(arg), cos(arg)
        dpls += (xls_nutation[i][6] + xls_nutation[i][7] * t) * sarg + xls_nutation[i][8] * carg
        dels += (xls_nutation[i][9] + xls_nutation[i][10] * t) * carg + xls_nutation[i][11] * sarg
    end

    dppl = 0.0
    depl = 0.0

    for i in reverse(eachindex(xpl_nutation))
        arg = mod2pi(xpl_nutation[i][1] * al +
             xpl_nutation[i][2] * af   +
             xpl_nutation[i][3] * ad   +
             xpl_nutation[i][4] * aom  +
             xpl_nutation[i][5] * alme +
             xpl_nutation[i][6] * alve +
             xpl_nutation[i][7] * alea +
             xpl_nutation[i][8] * alma +
             xpl_nutation[i][9] * alju +
             xpl_nutation[i][10] * alsa +
             xpl_nutation[i][11] * alur +
             xpl_nutation[i][12] * alne +
             xpl_nutation[i][13] * apa)
        sarg, carg = sin(arg), cos(arg)

        dppl += xpl_nutation[i][14] * sarg + xpl_nutation[i][15] * carg
        depl += xpl_nutation[i][16] * sarg + xpl_nutation[i][17] * carg
    end

    sec2rad(1e-7dpls) + sec2rad(1e-7dppl), sec2rad(1e-7dels) + sec2rad(1e-7depl)
end

"""
    nutation_a06(jd1, jd2)

Returns nutation, luni-solar and planetary for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> nutation_a06(2.4578265e6, 0.30434616919175345)
(-3.758986477837638e-5, -3.665903094375903e-5)
```
"""
function nutation_a06(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY
    dp, de = nutation_00a(jd1, jd2)

    dp + dp * (0.4697e-6 - 2.7774e-6 * t), de + de * (-2.7774e-6 * t)
end

"""
    greenwich_mean_sidereal_time06(uta, utb, tta, ttb, rnpb)

Returns Greenwich mean sidereal time(radians) for given two, 2 part Julian dates (TT and UT1).

# Example

```jldoctest
julia> gst06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, rand(3,3))
3.738379033673736
```
"""
function greenwich_apparent_sidereal_time06(uta, utb, tta, ttb, rnpb)
    x, y = rnpb[1, 3], rnpb[2, 3]
    s = s06(tta, ttb, x, y)
    era = earth_rotation_angle(uta, utb)
    eors = equation_of_origins(rnpb,s)
    mod2pi(era - eors)
end

"""
    equation_of_equinoxes_94(jd1, jd2)

Returns equation of equinoxes for a given 2 part Julian date (TT). (IAU 1994 model)

# Example

```jldoctest
julia> equation_of_equinoxes_94(2.4578265e6, 0.30434616919175345)
-3.446039061100115e-5
```
"""
function equation_of_equinoxes_94(jd1, jd2)
    t = ((jd1 - J2000) + jd2) / DAYS_PER_CENTURY

    om = rem2pi(sec2rad(@evalpoly t 450160.280 -482890.539 7.455 0.008)
      + mod(-5.0 * t, 1.0) * 2pi, RoundNearest)

    dpsi, deps = nutation(jd1, jd2)
    eps0 = mean_obliquity_of_ecliptic(jd1, jd2)

    dpsi*cos(eps0) + sec2rad((0.00264 * sin(om) + 0.000063 * sin(om + om)))
end

"""
    nutation_matrix80(jd1, jd2)

Returns nutation matrix for a given 2 part Julian date (TDB).

# Example

```jldoctest
julia> nutation_matrix80(2.4578265e6, 0.30434616919175345)
3×3 Rotations.RotXZX{Float64}(0.409054, 3.75653e-5, -0.409017):
 1.0         -3.44666e-5  -1.494e-5
 3.44661e-5   1.0         -3.66564e-5
 1.49413e-5   3.66559e-5   1.0
 ```
"""
function nutation_matrix80(jd1, jd2)
    dpsi, deps = nutation(jd1, jd2)
    epsa = mean_obliquity_of_ecliptic(jd1, jd2)
    numat(epsa, dpsi, deps)
end

"""
    precession_nutation00(jd1, jd2, dpsi, deps)

Returns
    epsa
    rb
    rp
    rbp
    rn
    rbpn
 for a given 2 part Julian date (TT) and nutations.

# Example

```jldoctest
julia> pn00(2.4578265e6, 0.30434616919175345, 0.2, 0.2)
(0.40905374831590824, [1.0 7.07828e-8 -8.05622e-8; -7.07828e-8 1.0 -3.30604e-8; 8.05622e-8 3.30604e-8 1.0], [0.999991 0.00384587 0.00167106; -0.00384587 0.999993 -3.24146e-6; -0.00167106 -3.18524e-6 0.999999], [0.999991 0.00384594 0.00167098; -0.00384594 0.999993 -3.27421e-6; -0.00167098 -3.1523e-6 0.999999], [0.980067 0.162947 0.113657; -0.182279 0.965066 0.188206; -0.079019 -0.205172 0.975531], [0.979225 0.166314 0.11601; -0.186046 0.964433 0.187765; -0.080656 -0.205447 0.975339])
```
"""
function precession_nutation00(jd1, jd2, dpsi, deps)
    dpsipr, depspr = precession_rate_part_of_nutation(jd1, jd2)
    epsa = mean_obliquity_of_ecliptic(jd1, jd2) + depspr
    rb, rp, rbpw = bias_precession_matrix_00(jd1, jd2)
    rnw = numat(epsa, dpsi, deps)
    epsa, rb, rp, rbpw, rnw, rbpw * rnw
end

"""
    precession_nutation_a00(jd1, jd2, dpsi, deps)

Returns
    dpsi, deps
    epsa
    rb
    rp
    rbp
    rn
    rbpn
for a given 2 part Julian date (TT) and nutations.

# Example

```jldoctest
julia> precession_nutation_a00(2.4578265e6, 0.30434616919175345)
(-3.758986507815969e-5, -3.6659048454874034e-5, (0.40905374831590824, [1.0 7.07828e-8 -8.05622e-8; -7.07828e-8 1.0 -3.30604e-8; 8.05622e-8 3.30604e-8 1.0], [0.999991 0.00384587 0.00167106; -0.00384587 0.999993 -3.24146e-6; -0.00167106 -3.18524e-6 0.999999], [0.999991 0.00384594 0.00167098; -0.00384594 0.999993 -3.27421e-6; -0.00167098 -3.1523e-6 0.999999], [1.0 -3.44892e-5 -1.49498e-5; 3.44886e-5 1.0 -3.66593e-5; 1.4951e-5 3.66588e-5 1.0], [0.999991 0.00381152 0.00165589; -0.00381145 0.999993 -3.98758e-5; -0.00165603 3.35641e-5 0.999999]))
```
"""
function precession_nutation_a00(jd1, jd2)
    dpsi, deps = nutation_00a(jd1, jd2)
    dpsi, deps, precession_nutation00(jd1, jd2, dpsi, deps)
end

"""
    precession_nutation_b00(jd1, jd2, dpsi, deps)

Returns
    dpsi, deps
    epsa
    rb
    rp
    rbp
    rn
    rbpn
for a given 2 part Julian date (TT) and nutations.

# Example

```jldoctest
julia> precession_nutation_b00(2.4578265e6, 0.30434616919175345)
(-3.758913944197044e-5, -3.665743003046017e-5, (0.40905374831590824, [1.0 7.07828e-8 -8.05622e-8; -7.07828e-8 1.0 -3.30604e-8; 8.05622e-8 3.30604e-8 1.0], [0.999991 0.00384587 0.00167106; -0.00384587 0.999993 -3.24146e-6; -0.00167106 -3.18524e-6 0.999999], [0.999991 0.00384594 0.00167098; -0.00384594 0.999993 -3.27421e-6; -0.00167098 -3.1523e-6 0.999999], [1.0 -3.44885e-5 -1.49495e-5; 3.44879e-5 1.0 -3.66577e-5; 1.49508e-5 3.66572e-5 1.0], [0.999991 0.00381152 0.00165589; -0.00381146 0.999993 -3.98741e-5; -0.00165603 3.35625e-5 0.999999]))
```
"""
function precession_nutation_b00(jd1, jd2)
    dpsi, deps = nutation_00b(jd1, jd2)
    dpsi, deps, precession_nutation00(jd1, jd2, dpsi, deps)
end


"""
    precession_nutation_matrix_a00(jd1, jd2)

Returns classical NBP matrix for a given 2 part Julian date(TT)

# Example

```jldoctest
julia> precession_nutation_matrix_a00(2.4578265e6, 0.30434616919175345)
3×3 Rotations.RotMatrix{3,Float64,9}:
  0.999991    0.00381152   0.00165589
 -0.00381145  0.999993    -3.98758e-5
 -0.00165603  3.35641e-5   0.999999
```
"""
function precession_nutation_matrix_a00(jd1, jd2)
    precession_nutation_a00(jd1, jd2)[3][6]
end

"""
    precession_nutation_matrix_b00(jd1, jd2)

Returns classical NBP matrix for a given 2 part Julian date(TT)

# Example

```jldoctest
julia> precession_nutation_matrix_b00(2.4578265e6, 0.30434616919175345)
3×3 Rotations.RotMatrix{3,Float64,9}:
  0.999991    0.00381152   0.00165589
 -0.00381146  0.999993    -3.98741e-5
 -0.00165603  3.35625e-5   0.999999
```
"""
function precession_nutation_matrix_b00(jd1, jd2)
    precession_nutation_b00(jd1, jd2)[3][6]
end

"""
    precession_nutation_matrix_a06(jd1, jd2)

Returns classical NBP matrix for a given 2 part Julian date(TT)

# Example

```jldoctest
julia> precession_nutation_matrix_a06(2.4578265e6, 0.30434616919175345)
3×3 Rotations.RotMatrix{3,Float64,9}:
  0.999991    0.00381151   0.00165589
 -0.00381145  0.999993    -3.98694e-5
 -0.00165602  3.35578e-5   0.999999
```
"""
function precession_nutation_matrix_a06(jd1, jd2)
    gamb, phib, psib, epsa = precession_fukushima_williams06(jd1, jd2)
    dp, de = nutation_a06(jd1, jd2)
    fukushima_williams_matrix(gamb, phib, psib + dp, epsa + de)
end
"""
    s00a(jd1, jd2)

Returns the CIO locator (radians) for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> s00a(2.4578265e6, 0.30434616919175345)
2.9183401976144617e-8
```
"""
function s00a(jd1, jd2)
    rbpn = precession_nutation_matrix_a00(jd1, jd2)
    s00(jd1, jd2, rbpn[1, 3], rbpn[2, 3])
end

"""
    s00b(jd1, jd2)

Returns the CIO locator (radians) for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> s00b(2.4578265e6, 0.30434616919175345)
2.9182068810695328e-8
```
"""
function s00b(jd1, jd2)
    rbpn = precession_nutation_matrix_b00(jd1, jd2)
    s00(jd1, jd2, rbpn[1,3], rbpn[2,3])
end

"""
    nutation_matrix_day_a00(jd1, jd2)

Returns nutation matrix for a given 2 part Julian date(TT).

# Example

```jldoctest
julia> precession_nutation_a00(2.4578265e6, 0.30434616919175345)[3][5]
3×3 Rotations.RotXZX{Float64}(0.409054, 3.75891e-5, -0.409017):
 1.0         -3.44885e-5  -1.49495e-5
 3.44879e-5   1.0         -3.66577e-5
 1.49508e-5   3.66572e-5   1.0
```
"""
function nutation_matrix_day_a00(jd1, jd2)
    precession_nutation_a00(jd1, jd2)[3][5]
end

"""
    nutation_matrix_day_b00(jd1, jd2)

Returns nutation matrix for a given 2 part Julian date(TT).

# Example

```jldoctest
julia> precession_nutation_b00(2.4578265e6, 0.30434616919175345)[3][5]
3×3 Rotations.RotXZX{Float64}(0.409054, 3.75891e-5, -0.409017):
 1.0         -3.44885e-5  -1.49495e-5
 3.44879e-5   1.0         -3.66577e-5
 1.49508e-5   3.66572e-5   1.0
```
"""
function nutation_matrix_day_b00(jd1, jd2)
    precession_nutation_b00(jd1, jd2)[3][5]
end

"""
    nutation_matrix_day(jd1,jd2)

Returns nutation matrix for a given 2 part Julian date(TT).

# Example

```jldoctest
julia> precession_nutation_b00(2.4578265e6, 0.30434616919175345)[3][5]
3×3 Rotations.RotXZX{Float64}(0.409054, 3.75891e-5, -0.409017):
 1.0         -3.44885e-5  -1.49495e-5
 3.44879e-5   1.0         -3.66577e-5
 1.49508e-5   3.66572e-5   1.0
```
"""
function nutation_matrix_day(jd1, jd2)
    eps = mean_obliquity_of_ecliptic(jd1, jd2)
    dp, de = nutation_a06(jd1, jd2)

    numat(eps, dp, de)
end

"""
    xys00a(jd1, jd2)

Returns Celestial intermediate pole (x, y) and the CIO locator for a given 2 part Julian date. (TT)

# Example

```jldoctest
julia> xys00a(2.4578265e6, 0.30434616919175345)
(0.001655885227051451, -3.987617878202287e-5, 2.9183401976144617e-8)
```
"""
function xys00a(jd1, jd2)
    rbpn = precession_nutation_matrix_a00(jd1, jd2)
    rbpn[1,3], rbpn[2,3], s00(jd1, jd2, rbpn[1,3], rbpn[2,3])
end

"""
xys00a(jd1, jd2)

Returns Celestial intermediate pole (x, y) and the CIO locator for a given 2 part Julian date. (TT)

# Example

```jldoctest
julia> xys00b(2.4578265e6, 0.30434616919175345)
(0.001655885521808701, -3.987456146931851e-5, 2.9182068810695328e-8)
```
"""
function xys00b(jd1, jd2)
    rbpn = precession_nutation_matrix_b00(jd1, jd2)
    rbpn[1,3], rbpn[2,3], s00(jd1, jd2, rbpn[1,3], rbpn[2,3])
end

"""
    s06a(jd1, jd2)

Returns the CIO locator for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> s06a(2.4578265e6, 0.30434616919175345)
2.917767075719633e-8
```
"""
function s06a(jd1, jd2)
    rnpb = precession_nutation_matrix_a06(jd1, jd2)
    s06(jd1, jd2, rnpb[1,3], rnpb[2,3])
end

"""
    equation_of_equinoxes_a00(jd1, jd2)

Returns equation of equinoxes for a given 2 part Julian date(TT).

# Example

```jldoctest
julia> equation_of_equinoxes_a00(2.4578265e6, 0.30434616919175345)
-3.4482238309105874e-5
```
"""
function equation_of_equinoxes_a00(jd1, jd2)
    dpsipr, depspr = precession_rate_part_of_nutation(jd1, jd2)
    epsa = mean_obliquity_of_ecliptic(jd1, jd2) + depspr
    dpsi, deps = nutation_00a(jd1, jd2)
    dpsi * cos(epsa) + equation_of_equinoxes_complementary_terms(jd1, jd2)
end

"""
    equation_of_equinoxes_b00(jd1, jd2)

Returns equation of equinoxes for a given 2 part Julian date(TT).

# Example

```jldoctest
julia> equation_of_equinoxes_b00(2.4578265e6, 0.30434616919175345)
-3.4482238309105874e-5
```
"""
function equation_of_equinoxes_b00(jd1, jd2)
    dpsipr, depspr = precession_rate_part_of_nutation(jd1, jd2)
    epsa = mean_obliquity_of_ecliptic(jd1, jd2) + depspr
    dpsi, deps = nutation_00b(jd1, jd2)
    dpsi * cos(epsa) + equation_of_equinoxes_complementary_terms(jd1, jd2)
end



"""
    greenwich_mean_sidereal_time_a06(uta, utb, tta, ttb, rnpb)

Returns greenwich apparent sidereal time for given two, 2 part Julian dates (UT1, TT).

# Example

```jldoctest

```
"""
function greenwich_mean_sidereal_time_a06(uta, utb, tta, ttb)
    rnpb = precession_nutation_matrix_a06(tta, ttb)
    greenwich_apparent_sidereal_time06(uta, utb, tta, ttb, rnpb)
end

"""
    equation_of_equinoxes_a06(jd1, jd2)

Returns equation of equinoxes for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> equation_of_equinoxes_a06(2.4578265e6, 0.30434616919175345)
-3.448290610741367e-5
```
"""
function equation_of_equinoxes_a06(jd1, jd2)
    gst06a = greenwich_mean_sidereal_time_a06(0.0, 0.0, jd1, jd2)
    gmst06 = greenwich_mean_sidereal_time06(0.0, 0.0, jd1, jd2)

    rem2pi(gst06a - gmst06, RoundNearest)
end

"""
    greenwich_apparent_sidereal_time94(uta, utb)

Returns greenwich appart sidereal time for given 2 part Julian date (UT1).

# Example

```jldoctest
julia> greenwich_apparent_sidereal_time94(2.4578265e6, 0.30434616919175345)
4.9160197844443445
```
"""
function greenwich_apparent_sidereal_time94(uta, utb)
    gmst82 = greenwich_mean_sidereal_time82(uta, utb)
    eqeq94 = equation_of_equinoxes_94(uta, utb)
    mod2pi(gmst82 + eqeq94)
end

"""
    greenwich_apparent_sidereal_time_a00(uta, utb, tta, ttb)

Returns Greenwich apparent sidereal time for given two, 2 part Julian date. (UT1, TT)

# Example

```jldoctest
julia> greenwich_apparent_sidereal_time_a00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
4.959632360790386
```
"""
function greenwich_apparent_sidereal_time_a00(uta, utb, tta, ttb)
    gmst00 = greenwich_mean_sidereal_time00(uta, utb, tta, ttb)
    ee00a = equation_of_equinoxes_a00(tta, ttb)
    mod2pi(gmst00 + ee00a)
end

"""
    greenwich_apparent_sidereal_time_b00(uta, utb, tta, ttb)

Returns Greenwich apparent sidereal time for given two, 2 part Julian date. (UT1, TT)

# Example

```jldoctest
julia> greenwich_apparent_sidereal_time_b00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
4.95963236232993
```
"""
function greenwich_apparent_sidereal_time_b00(uta, utb)
    gmst00 = greenwich_mean_sidereal_time00(uta, utb, uta, utb)
    ee00b = equation_of_equinoxes_b00(uta, utb)
    mod2pi(gmst00 + ee00b)
end

"""
    greenwich_apparent_sidereal_time_a06(uta, utb, tta, ttb)

Returns Greenwich apparent sidereal time for given two, 2 part Julian dates. (UT1, TT)

# Example

```jldoctest
julia> greenwich_apparent_sidereal_time_a06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851)
4.959632359298287
```
"""
function greenwich_apparent_sidereal_time_a06(uta, utb, tta, ttb)
    rnpb = precession_nutation_matrix_a06(tta, ttb)
    greenwich_apparent_sidereal_time06(uta, utb, tta, ttb, rnpb)
end

"""
    celestial_to_intermediate_frame_of_date(jd1, jd2, x, y)

Returns celestial to intermediate-frame-of-date matrix for a given 2 part Julian Date.

# Example

```jldoctest
julia> celestial_to_intermediate_frame_of_date(2.4578265e6, 0.30434616919175345, 0.2, 0.3)
3×3 Rotations.RotZYZ{Float64}(0.982794, 0.368863, -0.952794):
  0.977932    -0.0604048  0.2
 -0.00243134   0.953936   0.3
 -0.208909    -0.293866   0.932738
```
"""
function celestial_to_intermediate_frame_of_date(jd1, jd2, x, y)
    celestial_to_intermediate(x, y, s00(jd1, jd2, x, y))
end

"""
    celestial_to_intermediate_matrix(jd1, jd2, rbpn)

Returns celestial-to-intermediate matrix for given 2 part Julian date and celestial-to-true matrix.

# Example

```jldoctest
julia> celestial_to_intermediate_matrix(2.4578265e6, 0.30434616919175345, rand(3,3))
3×3 Rotations.RotZYZ{Float64}(0.921674, 1.10597, -0.729284):
  0.732866  -0.413437  0.54035
 -0.13652    0.688686  0.712092
 -0.666536  -0.595636  0.448272
```
"""
function celestial_to_intermediate_matrix(jd1, jd2, rbpn)
    x, y = rbpn[1,3], rbpn[2,3]
    celestial_to_intermediate_frame_of_date(jd1, jd2, x, y)
end

"""
    celestial_to_intermediate_matrix_a00(jd1, jd2)

Returns celestial-to-intermediate matrix for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> celestial_to_intermediate_matrix_a00(2.4578265e6, 0.30434616919175345)
3×3 Rotations.RotZYZ{Float64}(-0.0240766, 0.00165637, 0.0240765):
  0.999999    6.21979e-8   0.00165589
  3.83181e-9  1.0         -3.98758e-5
 -0.00165589  3.98757e-5   0.999999
```
"""
function celestial_to_intermediate_matrix_a00(jd1, jd2)
    rbpn = precession_nutation_matrix_a00(jd1, jd2)
    celestial_to_intermediate_matrix(jd1, jd2, rbpn)
end

"""
    celestial_to_intermediate_matrix_b00(jd1, jd2)

Returns celestial-to-intermediate matrix for a given 2 part Julian date (TT).

# Example

```jldoctest
julia> celestial_to_intermediate_matrix_b00(2.4578265e6, 0.30434616919175345)
3×3 Rotations.RotZYZ{Float64}(-0.0240766, 0.00165637, 0.0240765):
0.999999    6.21979e-8   0.00165589
3.83181e-9  1.0         -3.98758e-5
-0.00165589  3.98757e-5   0.999999
```
"""
function celestial_to_intermediate_matrix_b00(jd1, jd2)
    rbpn = precession_nutation_matrix_b00(jd1, jd2)
    celestial_to_intermediate_matrix(jd1, jd2, rbpn)
end

"""
    celestial_to_terrestrial_matrix(tta, ttb, uta, utb, x, y, xp, yp)

Returns celestial to terrestrial matrix for given two, 2 part Julian dates (UT1, TT) and CIP coordinates.

# Example

```jldoctest
julia> celestial_to_terrestrial_matrix(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3, 0.5, 0.7)
3×3 Rotations.RotMatrix{3,Float64,9}:
  0.279858   0.75876   0.588186
 -0.66989   -0.284531  0.685777
  0.687697  -0.58594   0.428658
```
"""
function celestial_to_terrestrial_matrix(tta, ttb, uta, utb, x, y, xp, yp)
    rc2i = celestial_to_intermediate_frame_of_date(tta, ttb, x, y)
    era = earth_rotation_angle(uta, utb)
    sp = tio_locator(tta, ttb)

    rpom = polar_motion(xp, yp, sp)
    compose_rotation(rpom, angleaxis_to_dcm(-era, [0.0, 0.0, 1.0]), rc2i)
end

"""
    celestial_to_terrestrial_a00(tta, ttb, uta, utb, xp, yp)

Returns celestial to terrestrial matrix given two, 2 part Julian dates (UT1, TT) and coordinates of the pole.

# Example

```jldoctest
julia> celestial_to_terrestrial_a00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3)
3×3 Rotations.RotMatrix{3,Float64,9}:
  0.235088    0.934797  0.266247
 -0.966879    0.196899  0.162409
  0.0993959  -0.295609  0.950124
```
"""
function celestial_to_terrestrial_a00(tta, ttb, uta, utb, xp, yp)
    rc2i = celestial_to_intermediate_matrix_a00(tta, ttb)
    era = earth_rotation_angle(uta, utb)
    sp = tio_locator(tta, ttb)

    rpom = polar_motion(xp, yp, sp)
    compose_rotation(rpom, angleaxis_to_dcm(-era, [0.0, 0.0, 1.0]), rc2i)
end

"""
    celestial_to_terrestrial_b00(tta, ttb, uta, utb, xp, yp)

Returns celestial to terrestrial matrix given two, 2 part Julian dates (UT1, TT) and coordinates of the pole.

# Example

```jldoctest
julia> celestial_to_terrestrial_b00(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3)
3×3 Rotations.RotMatrix{3,Float64,9}:
  0.235088    0.934797  0.266247
 -0.966879    0.196899  0.162409
  0.0993959  -0.295609  0.950124
```
"""
function celestial_to_terrestrial_b00(tta, ttb, uta, utb, xp, yp)
    rc2i = celestial_to_intermediate_matrix_b00(tta, ttb)
    era = earth_rotation_angle(uta, utb)

    rpom = polar_motion(xp, yp, 0.0)
    compose_rotation(rpom, angleaxis_to_dcm(-era, [0.0, 0.0, 1.0]), rc2i)
end

"""
    c2tpe(tta, ttb, uta, utb, dpsi, deps, xp, yp)

Returns the celestial to terrestrial matrix given two, 2 part Julian date (UT1, TT), the nutation and coordinates of the pole(radians).

# Example

```jldoctest
julia> c2tpe(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, 0.1, 0.3, 0.5, 0.7)
3×3 Rotations.RotMatrix{3,Float64,9}:
  0.253143   0.777796  0.575284
 -0.676266  -0.28296   0.680146
  0.691797  -0.561219  0.454367
```
"""
function c2tpe(tta, ttb, uta, utb, dpsi, deps, xp, yp)
    epsa, rb, rp, rbp, rn, rbpn = precession_nutation00(tta, ttb, dpsi, deps)
    gmst = greenwich_mean_sidereal_time00(uta, utb, tta, ttb)
    ee = equation_of_equinoxes_00(tta, ttb, epsa, dpsi)
    sp = tio_locator(tta, ttb)

    rpom = polar_motion(xp, yp, sp)
    compose_rotation(rpom, angleaxis_to_dcm(-(gmst + ee), [0.0, 0.0, 1.0]), rbpn)
end

end
