module AstroBase

using Reexport
using Rotations

include("util.jl")
include(joinpath("bodies", "Bodies.jl"))
include("EarthAttitude.jl")

@reexport using .Bodies
@reexport using .EarthAttitude

export nutation_matrix80,
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
    rnw = numat(epsa,dpsi, deps)
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
    precession_nutation_matrix_a00

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
    precession_nutation_matrix_a00

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
approx error
"""
function s00a(jd1, jd2)
    rbpn = precession_nutation_matrix_a00(jd1, jd2)
    s00(jd1, jd2, rbpn[1, 3], rbpn[2, 3])
end

"""
approx error
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
approx error
"""
function xys00a(jd1, jd2)
    rbpn = precession_nutation_matrix_a00(jd1, jd2)
    rbpn[1,3], rbpn[2,3], s00(jd1, jd2, rbpn[1,3], rbpn[2,3])
end

"""
approx error
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
    rc2i * RotZ(era) * rpom
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
    rc2i * RotZ(era) * rpom
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
    rc2i * RotZ(era) * rpom
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
    rbpn * RotZ(gmst + ee) * rpom
end
end # module
