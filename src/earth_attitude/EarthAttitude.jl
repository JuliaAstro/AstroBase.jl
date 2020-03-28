module EarthAttitude

using AstroTime: days
using ..Bodies
using ..Util: sec2rad
using ReferenceFrameRotations: angle_to_dcm, angleaxis_to_dcm, compose_rotation

export
    J2000,
    celestial_to_intermediate,
    earth_rotation_angle,
    equation_of_equinoxes_00,
    equation_of_equinoxes_complementary_terms,
    equation_of_origins,
    general_precession_in_longitude,
    greenwich_mean_sidereal_time00,
    greenwich_mean_sidereal_time06,
    greenwich_mean_sidereal_time82,
    mean_elongation,
    mean_longitude,
    mean_longitude_ascending_node,
    mean_longitude_minus_lan,
    polar_motion,
    s00,
    s06,
    tio_locator,
    xy06,
    greenwich_apparent_sidereal_time06,
    equation_of_equinoxes_94,
    s00a,
    s00b,
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
include("obliquity.jl")
include("nutation.jl")
include("precession.jl")
include("precession_nutation.jl")
include("icrs.jl")

const J2000 = 2451545.0
const DAYS_PER_CENTURY = 36525.0
const ARCSECONDS_IN_CIRCLE = 1296000.0
const SECONDS_PER_DAY = 24.0 * 60.0 * 60.0

include(joinpath("constants", "mfals.jl"))
include(joinpath("constants", "cio_locator.jl"))
include(joinpath("constants", "EE00.jl"))

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
    angle_to_dcm(e, d, -(e+s), :ZYZ)
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
function polar_motion(xp, yp, sp)
    angle_to_dcm(sp, -xp, -yp, :ZYX)
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
    equation_of_origins(rnpb, s)

Returns the equation of origins(radians) for given nutation-bias-precession matrix and the CIO locator.

 # Example

 ```jldoctest
equation_of_origins(rand(3,3), 0.2)
1.7738370040531068
 ```
"""
function equation_of_origins(rnpb, s)
    x = rnpb[3,1]
    ax = x / (1.0 + rnpb[3,3])
    xs = 1.0 - ax * x
    ys = -ax * rnpb[3,2]
    zs = -x
    p = rnpb[1,1] * xs + rnpb[1,2] * ys + rnpb[1,3] * zs
    q = rnpb[2,1] * xs + rnpb[2,2] * ys + rnpb[2,3] * zs
    p != 0 || q != 0 ? s - atan(q, p) : s
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
    greenwich_mean_sidereal_time06(uta, utb, tta, ttb, rnpb)

Returns Greenwich mean sidereal time(radians) for given two, 2 part Julian dates (TT and UT1).

# Example

```jldoctest
julia> gst06(2.4579405e6, 0.0, 2.4579405e6, -0.0007966009351851851, rand(3,3))
3.738379033673736
```
"""
function greenwich_apparent_sidereal_time06(uta, utb, tta, ttb, rnpb)
    x, y = cip_coords(rnpb)
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

    dpsi, deps = nutation(iau1980, TTEpoch(jd1 * days, jd2 * days, origin=:julian))
    eps0 = obliquity(iau1980, TTEpoch(jd1 * days, jd2 * days, origin=:julian))

    dpsi*cos(eps0) + sec2rad((0.00264 * sin(om) + 0.000063 * sin(om + om)))
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
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    rbpn = precession_nutation_matrix(iau2000a, ep)
    x, y = cip_coords(rbpn)
    s00(jd1, jd2, x, y)
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
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    rbpn = precession_nutation_matrix(iau2000b, ep)
    x, y = cip_coords(rbpn)
    s00(jd1, jd2, x, y)
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
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    rbpn = precession_nutation_matrix(iau2000a, ep)
    x, y = cip_coords(rbpn)
    return x, y, s00(jd1, jd2, x, y)
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
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    rbpn = precession_nutation_matrix(iau2000b, ep)
    x, y = cip_coords(rbpn)
    return x, y, s00(jd1, jd2, x, y)
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
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    rbpn = precession_nutation_matrix(iau2006a, ep)
    s06(jd1, jd2, cip_coords(rbpn)...)
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
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    dpsipr, depspr = precession(iau2000, ep)
    epsa = obliquity(iau1980, ep) + depspr
    dpsi, deps = nutation(iau2000a, ep)
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
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    dpsipr, depspr = precession(iau2000, ep)
    epsa = obliquity(iau1980, ep) + depspr
    dpsi, deps = nutation(iau2000b, ep)
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
    ep = TTEpoch(tta * days, ttb * days, origin=:julian)
    rnpb = precession_nutation_matrix(iau2006a, ep)
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
    ep = TTEpoch(tta * days, ttb * days, origin=:julian)
    rnpb = precession_nutation_matrix(iau2006a, ep)
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
    x, y = cip_coords(rbpn)
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
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    rbpn = precession_nutation_matrix(iau2000a, ep)
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
    ep = TTEpoch(jd1 * days, jd2 * days, origin=:julian)
    rbpn = precession_nutation_matrix(iau2000b, ep)
    celestial_to_intermediate_matrix(jd1, jd2, rbpn)
end

function c2tcio(rc2i, era, rpom)
    rera = angleaxis_to_dcm(era, [0.0, 0.0, 1.0])
    return compose_rotation(rc2i, rera, rpom)
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
    return c2tcio(rc2i, era, rpom)
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
    return c2tcio(rc2i, era, rpom)
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
    return c2tcio(rc2i, era, rpom)
end

function c2teqx(rbpn, gst, rpom)
    rgst = angleaxis_to_dcm(gst, [0.0, 0.0, 1.0])
    return compose_rotation(rbpn, rgst, rpom)
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
    tt = TTEpoch(tta * days, ttb * days, origin=:julian)
    epsa, rb, rp, rbp, rn, rbpn = precession_nutation(iau2000, tt, dpsi, deps)
    gmst = greenwich_mean_sidereal_time00(uta, utb, tta, ttb)
    ee = equation_of_equinoxes_00(tta, ttb, epsa, dpsi)
    sp = tio_locator(tta, ttb)

    rpom = polar_motion(xp, yp, sp)
    return c2teqx(rbpn, gmst + ee, rpom)
end

end
